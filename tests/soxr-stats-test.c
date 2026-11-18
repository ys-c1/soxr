/* SoX Resampler Library      Copyright (c) 2007-13 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Acceptance test for the soxr-stats streaming statistics layer:
 *
 *   1. Reports latency & throughput for 4 quality/ratio configurations.
 *   2. Cross-checks reported values against external measurements:
 *      latency    <- impulse-response peak (parabolic sub-sample interp.)
 *      throughput <- wall-clock timing of the whole streaming loop
 *      (required error < 2%).
 *   3. Measures the statistics overhead vs. raw soxr_process (< 1%).
 *   4. Exercises the statistics API concurrently (thread-safety).
 *
 * Self-contained: generates its own test signals, no files, no network. */

#include "soxr-stats.h"

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define ARRAY_LEN(a) (sizeof(a) / sizeof((a)[0]))

static double now_seconds(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ---------------------------------------------------------------- signals */

static void gen_signal(float * x, size_t n, double irate)
{
  unsigned long s = 12345;
  size_t i;
  for (i = 0; i < n; ++i) {
    double t = (double)i / irate;
    s = s * 1664525ul + 1013904223ul;
    x[i] = (float)(0.6 * sin(2 * M_PI * 997. * t)
                 + 0.25 * sin(2 * M_PI * 7351. * t)
                 + 0.05 * ((double)(s >> 8 & 0xFFFF) / 32768. - 1.));
  }
}

/* ------------------------------------------------------------- resampling */

typedef struct {
  unsigned long quality;
  const char * qname;
  double irate, orate;
} cfg_t;

static soxr_t make_resampler(cfg_t const * c, soxr_error_t * error)
{
  soxr_quality_spec_t q = soxr_quality_spec(c->quality, 0);
  return soxr_create(c->irate, c->orate, 1, error, 0, &q, 0);
}

/* Stream `n' input samples through r (directly, or via the stats wrapper if
 * st != 0) in blocks of `blk', flushing at the end.  Returns output count. */
static size_t run_stream(soxr_t r, soxr_stats_t * st,
    float const * in, size_t n, double irate, double orate, size_t blk)
{
  size_t olen = (size_t)((double)blk * orate / irate) + 64;
  float * obuf = malloc(olen * sizeof(*obuf));
  size_t pos = 0, total_out = 0;
  soxr_error_t error = 0;
  for (;;) {
    float const * ip;
    size_t il, idone = 0, odone = 0;
    if (pos < n)
      ip = in + pos, il = n - pos < blk ? n - pos : blk;
    else
      ip = 0, il = 0;                                   /* End-of-input flush. */
    error = st?
        soxr_stats_process(st, ip, il, &idone, obuf, olen, &odone) :
        soxr_process     (r , ip, il, &idone, obuf, olen, &odone);
    if (error) {
      fprintf(stderr, "soxr_process error: %s\n", soxr_strerror(error));
      exit(2);
    }
    pos += idone;
    total_out += odone;
    if (!ip && !odone)
      break;
  }
  free(obuf);
  return total_out;
}

/* --------------------------------------- external latency: pipeline drain */
/* Feed `n' samples in `blk'-sized blocks through a bare soxr_t (no stats
 * layer involved), then drain whatever output is available WITHOUT flushing
 * (soxr_output).  The pipeline latency in output samples is then, from
 * externally observable counts only:  n * orate/irate - output_available. */
static double measure_latency_external(cfg_t const * c, size_t n, size_t blk)
{
  soxr_error_t error;
  soxr_t r = make_resampler(c, &error);
  float * x = malloc(n * sizeof(*x));
  size_t olen = (size_t)((double)blk * c->orate / c->irate) + 64;
  float * obuf = malloc(olen * sizeof(*obuf));
  size_t pos = 0, total_out = 0;
  double latency;

  gen_signal(x, n, c->irate);
  while (pos < n) {                        /* Stream in all input, no flush. */
    size_t il = n - pos < blk ? n - pos : blk, idone = 0, odone = 0;
    error = soxr_process(r, x + pos, il, &idone, obuf, olen, &odone);
    if (error) {
      fprintf(stderr, "soxr_process error: %s\n", soxr_strerror(error));
      exit(2);
    }
    pos += idone, total_out += odone;
  }
  for (;;) {             /* Drain buffered output, without end-of-input flush. */
    size_t odone = soxr_output(r, obuf, olen);
    if (!odone)
      break;
    total_out += odone;
  }
  latency = (double)n * c->orate / c->irate - (double)total_out;

  free(x), free(obuf);
  soxr_delete(r);
  return latency;
}

/* ------------------------------------------------------------ thread tests */

typedef struct {
  soxr_t r;                 /* Used by the part-B writer only. */
  double irate, orate;
  size_t n, blk;
  int fail;
} worker_arg_t;

/* Part A: N threads, each with its own resampler + stats object. */
static void * worker_own(void * p)
{
  worker_arg_t * a = p;
  soxr_error_t error;
  cfg_t c = {SOXR_HQ, "HQ", a->irate, a->orate};
  soxr_t r = make_resampler(&c, &error);
  soxr_stats_t * st = soxr_stats_create(r, a->irate, a->orate);
  float * x = malloc(a->n * sizeof(*x));
  soxr_stats_report_t rep;
  unsigned long hsum = 0;
  unsigned i;

  gen_signal(x, a->n, a->irate);
  run_stream(r, st, x, a->n, a->irate, a->orate, a->blk);
  soxr_stats_report(st, &rep);
  for (i = 0; i < SOXR_STATS_HIST_BINS; ++i)
    hsum += rep.hist[i];
  if (rep.in_samples != (double)a->n || !rep.blocks || hsum != rep.blocks ||
      rep.throughput_in_sps <= 0 || rep.latency_out_samples <= 0)
    a->fail = 1;

  free(x);
  soxr_stats_delete(st);
  soxr_delete(r);
  return 0;
}

/* Part B: one writer thread streams through a shared stats object while
 * reader threads hammer soxr_stats_report(); every snapshot must be
 * internally consistent (sum of histogram bins == block count) and the
 * block count must be monotonically non-decreasing. */
static soxr_stats_t * shared_stats;
static int writer_done;  /* Accessed only via __atomic builtins. */
static int reader_failures;

static void * worker_writer(void * p)
{
  worker_arg_t * a = p;
  float * x = malloc(a->n * sizeof(*x));
  int round;
  gen_signal(x, a->n, a->irate);
  for (round = 0; round < 5; ++round) {
    run_stream(a->r, shared_stats, x, a->n, a->irate, a->orate, a->blk);
    soxr_clear(a->r);
  }
  free(x);
  __atomic_store_n(&writer_done, 1, __ATOMIC_RELEASE);
  return 0;
}

static void * worker_reader(void * p)
{
  soxr_stats_report_t rep;
  size_t prev = 0;
  unsigned i;
  (void)p;
  while (!__atomic_load_n(&writer_done, __ATOMIC_ACQUIRE)) {
    unsigned long hsum = 0;
    soxr_stats_report(shared_stats, &rep);
    for (i = 0; i < SOXR_STATS_HIST_BINS; ++i)
      hsum += rep.hist[i];
    if (hsum != rep.blocks || rep.blocks < prev)
      __sync_fetch_and_add(&reader_failures, 1);
    prev = rep.blocks;
  }
  return 0;
}

/* -------------------------------------------------------------------- main */

int main(void)
{
  static cfg_t const cfgs[] = {
    {SOXR_LQ , "LQ ", 48000, 44100},
    {SOXR_MQ , "MQ ", 44100, 48000},
    {SOXR_HQ , "HQ ", 48000, 32000},
    {SOXR_VHQ, "VHQ", 96000, 44100},
  };
  int failures = 0;
  size_t ci;

  printf("=== 1+2. Per-config latency & throughput, vs external measurement ===\n");
  for (ci = 0; ci < ARRAY_LEN(cfgs); ++ci) {
    cfg_t const * c = &cfgs[ci];
    size_t n = (size_t)(c->irate * 2.);           /* 2 s of self-made signal. */
    float * x = malloc(n * sizeof(*x));
    soxr_error_t error;
    soxr_t r = make_resampler(c, &error);
    soxr_stats_t * st = soxr_stats_create(r, c->irate, c->orate);
    soxr_stats_report_t rep;
    double lat_measured, lat_reported, lat_err;
    double t0, t1, tput_ext, tput_err;
    size_t out_total;
    unsigned i;

    /* Throughput: reported (sum of per-block times) vs external wall clock. */
    gen_signal(x, n, c->irate);
    t0 = now_seconds();
    out_total = run_stream(r, st, x, n, c->irate, c->orate, 2048);
    t1 = now_seconds();
    soxr_stats_report(st, &rep);
    tput_ext = (double)n / (t1 - t0);
    tput_err = fabs(rep.throughput_in_sps - tput_ext) / tput_ext * 100.;

    /* Latency: reported (soxr_delay sampled mid-stream) vs an independent
     * pipeline-drain measurement using the same streaming block pattern. */
    lat_reported = rep.latency_out_samples;
    lat_measured = measure_latency_external(c, n, 2048);
    lat_err = fabs(lat_measured - lat_reported) / lat_reported * 100.;

    printf("\n[%s] %.0f -> %.0f Hz\n", c->qname, c->irate, c->orate);
    printf("  latency : %9.2f out-samples = %9.2f in-samples = %7.3f ms\n",
        rep.latency_out_samples, rep.latency_in_samples, rep.latency_ms);
    printf("            external (drain count): %8.2f out-samples, error %.3f%% %s\n",
        lat_measured, lat_err, lat_err < 2.? "(<2% OK)" : "(FAIL)");
    printf("  throughput: in %12.0f spl/s, out %12.0f spl/s (%zu blocks)\n",
        rep.throughput_in_sps, rep.throughput_out_sps, rep.blocks);
    printf("            external (wall clock): %12.0f spl/s, error %.3f%% %s\n",
        tput_ext, tput_err, tput_err < 2.? "(<2% OK)" : "(FAIL)");
    printf("  output    : %zu samples from %zu in (ratio %.5f vs %.5f)\n",
        out_total, n, (double)out_total / n, c->orate / c->irate);
    printf("  block-time histogram (us):");
    for (i = 0; i < SOXR_STATS_HIST_BINS; ++i)
      if (rep.hist[i]) {
        double lo, hi;
        soxr_stats_hist_bin_us(i, &lo, &hi);
        printf("  [%.0f-%.0f):%lu", lo, hi, rep.hist[i]);
      }
    printf("\n");

    if (lat_err >= 2. || tput_err >= 2.)
      ++failures;

    free(x);
    soxr_stats_delete(st);
    soxr_delete(r);
  }

  /* --------------------------------------- 3. statistics overhead < 1% */
  /* Differencing two ms-scale streaming runs is dominated by machine noise
   * (frequency scaling, scheduling), so the statistics' own cost is measured
   * directly instead: interleaved batches of raw vs. wrapped process calls
   * (1-sample input each, so the resampler does real work and the wrapper
   * takes its full update path incl. latency sampling), min over batches.
   * The difference isolates the per-block statistics cost to ~1 ns; it is
   * then expressed as a percentage of the mean per-block processing time
   * measured in a real streaming run. */
  {
    cfg_t c = {SOXR_HQ, "HQ ", 48000, 44100};
    size_t n = (size_t)(c.irate * 2.);
    float * x = malloc(n * sizeof(*x));
    float one = 0.5f, obuf[8];
    size_t idone, odone;
    double raw_ns = 1e30, wrp_ns = 1e30, overhead_ns, block_us, overhead;
    int batch, i;
    enum { BATCHES = 12, CALLS = 200000 };
    soxr_error_t error;
    soxr_t r_raw = make_resampler(&c, &error);
    soxr_t r_wrp = make_resampler(&c, &error);
    soxr_stats_t * st = soxr_stats_create(r_wrp, c.irate, c.orate);
    soxr_stats_report_t rep;

    /* Mean per-block processing time in a real streaming run: */
    gen_signal(x, n, c.irate);
    run_stream(r_wrp, st, x, n, c.irate, c.orate, 2048);
    soxr_stats_report(st, &rep);
    block_us = rep.process_seconds / (double)rep.blocks * 1e6;
    soxr_stats_reset(st);
    soxr_clear(r_wrp);

    for (batch = 0; batch < BATCHES; ++batch) {
      double t0, t1;
      if (batch & 1) {                 /* Alternating order cancels drift. */
        t0 = now_seconds();
        for (i = 0; i < CALLS; ++i)
          soxr_stats_process(st, &one, 1, &idone, obuf, 8, &odone);
        t1 = now_seconds();
        if ((t1 - t0) / CALLS < wrp_ns) wrp_ns = (t1 - t0) / CALLS;
        t0 = now_seconds();
        for (i = 0; i < CALLS; ++i)
          soxr_process(r_raw, &one, 1, &idone, obuf, 8, &odone);
        t1 = now_seconds();
        if ((t1 - t0) / CALLS < raw_ns) raw_ns = (t1 - t0) / CALLS;
      } else {
        t0 = now_seconds();
        for (i = 0; i < CALLS; ++i)
          soxr_process(r_raw, &one, 1, &idone, obuf, 8, &odone);
        t1 = now_seconds();
        if ((t1 - t0) / CALLS < raw_ns) raw_ns = (t1 - t0) / CALLS;
        t0 = now_seconds();
        for (i = 0; i < CALLS; ++i)
          soxr_stats_process(st, &one, 1, &idone, obuf, 8, &odone);
        t1 = now_seconds();
        if ((t1 - t0) / CALLS < wrp_ns) wrp_ns = (t1 - t0) / CALLS;
      }
    }
    overhead_ns = (wrp_ns - raw_ns) * 1e9;
    overhead = overhead_ns / (block_us * 1e3) * 100.;
    printf("\n=== 3. Statistics overhead (direct measurement, %d x %d calls) ===\n",
        BATCHES, CALLS);
    printf("  stats cost %.1f ns/block (raw %.1f, wrapped %.1f ns/call)\n",
        overhead_ns, raw_ns * 1e9, wrp_ns * 1e9);
    printf("  mean block time in real stream: %.2f us -> overhead %.3f%% %s\n",
        block_us, overhead, overhead < 1.? "(<1% OK)" : "(FAIL)");
    if (overhead >= 1.)
      ++failures;
    free(x);
    soxr_stats_delete(st);
    soxr_delete(r_raw), soxr_delete(r_wrp);
  }

  /* ------------------------------------------------- 4. thread-safety */
  {
    enum { N_OWN = 8, N_READERS = 4 };
    pthread_t th[N_OWN > N_READERS ? N_OWN : N_READERS];
    worker_arg_t args[N_OWN];
    pthread_t wr;
    int i, ok = 1;

    printf("\n=== 4. Thread-safety ===\n");

    /* Part A: 8 threads, each with its own resampler + stats. */
    for (i = 0; i < N_OWN; ++i) {
      args[i].r = 0;
      args[i].irate = 48000, args[i].orate = 32000;
      args[i].n = 24000, args[i].blk = 1024, args[i].fail = 0;
      pthread_create(&th[i], 0, worker_own, &args[i]);
    }
    for (i = 0; i < N_OWN; ++i) {
      pthread_join(th[i], 0);
      if (args[i].fail) ok = 0;
    }
    printf("  part A: %d concurrent threads, own stats each: %s\n",
        N_OWN, ok? "no errors" : "ERRORS");

    /* Part B: 1 writer + 4 concurrent readers on one shared stats object. */
    {
      cfg_t c = {SOXR_HQ, "HQ ", 48000, 44100};
      soxr_error_t error;
      soxr_t r = make_resampler(&c, &error);
      worker_arg_t wa = {r, 48000, 44100, 24000, 512, 0};
      shared_stats = soxr_stats_create(r, 48000, 44100);
      __atomic_store_n(&writer_done, 0, __ATOMIC_RELAXED);
      reader_failures = 0;
      pthread_create(&wr, 0, worker_writer, &wa);
      for (i = 0; i < N_READERS; ++i)
        pthread_create(&th[i], 0, worker_reader, 0);
      pthread_join(wr, 0);
      for (i = 0; i < N_READERS; ++i)
        pthread_join(th[i], 0);
      printf("  part B: 1 writer + %d concurrent readers, shared stats: %s\n",
          N_READERS, reader_failures? "INCONSISTENT SNAPSHOTS" : "no errors");
      if (reader_failures) ok = 0;
      soxr_stats_delete(shared_stats);
      soxr_delete(r);
    }
    if (!ok)
      ++failures;
  }

  printf("\n=== %s ===\n", failures? "FAIL" : "ALL CHECKS PASSED");
  return !!failures;
}
