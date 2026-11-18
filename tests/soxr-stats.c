/* SoX Resampler Library      Copyright (c) 2007-13 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Implementation of the streaming statistics layer declared in soxr-stats.h.
 * Only soxr's public API is used.  Per-block statistics overhead is two
 * timestamp reads, one short mutex hold and a few register updates; on x86
 * the time-stamp counter is used directly (~7ns/read) to keep the overhead
 * far below 1% of typical per-block processing time, with a portable
 * clock_gettime(CLOCK_MONOTONIC) fallback elsewhere. */

#include "soxr-stats.h"

#include <math.h>           /* ldexp (display helper only, not the hot path) */
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
  #define SOXR_STATS_USE_TSC 1
  #include <x86intrin.h>
#endif

struct soxr_stats {
  soxr_t resampler;          /* Not owned. */
  double irate, orate;
  pthread_mutex_t mutex;
  /* Everything below is guarded by `mutex': */
  size_t blocks;
  double in_samples, out_samples;
  double process_seconds;
  double latency_out;         /* soxr_delay() sampled while consuming input. */
  int latency_set;
  unsigned long hist[SOXR_STATS_HIST_BINS];
};

#if SOXR_STATS_USE_TSC
static double tsc_seconds;  /* Seconds per TSC cycle, calibrated once. */
static pthread_once_t tsc_once = PTHREAD_ONCE_INIT;

static void tsc_calibrate(void)
{
  struct timespec ts0, ts1, req = {0, 4000000};  /* 4 ms calibration window. */
  unsigned long long c0, c1;
  clock_gettime(CLOCK_MONOTONIC, &ts0);
  c0 = __rdtsc();
  nanosleep(&req, 0);
  c1 = __rdtsc();
  clock_gettime(CLOCK_MONOTONIC, &ts1);
  tsc_seconds = ((double)(ts1.tv_sec - ts0.tv_sec) +
      1e-9 * (double)(ts1.tv_nsec - ts0.tv_nsec)) / (double)(c1 - c0);
}

static double now_seconds(void)
{
  return (double)__rdtsc() * tsc_seconds;
}
#else
static double now_seconds(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}
#endif

static unsigned hist_bin(double seconds)
{
  /* Bin b (b>0) covers [2^(b-1), 2^b) us; bin 0 covers <1us.  Integer-only
   * (no libm in the hot path): floor(log2(us)) = 63 - clzll(us). */
  unsigned long long us = (unsigned long long)(seconds * 1e6);
  unsigned b;
  if (!us)
    return 0;
#if defined(__GNUC__) || defined(__clang__)
  b = 63u - (unsigned)__builtin_clzll(us) + 1u;
#else
  {
    unsigned long long v = us;
    for (b = 1; v >>= 1;)
      ++b;
  }
#endif
  return b < SOXR_STATS_HIST_BINS ? b : SOXR_STATS_HIST_BINS - 1;
}

soxr_stats_t * soxr_stats_create(soxr_t resampler, double irate, double orate)
{
  soxr_stats_t * s = malloc(sizeof(*s));
  if (!s)
    return 0;
#if SOXR_STATS_USE_TSC
  pthread_once(&tsc_once, tsc_calibrate);  /* One-off, outside any timing. */
#endif
  s->resampler = resampler;
  s->irate = irate;
  s->orate = orate;
  pthread_mutex_init(&s->mutex, 0);
  s->blocks = 0;
  s->in_samples = s->out_samples = s->process_seconds = 0;
  s->latency_out = 0, s->latency_set = 0;
  memset(s->hist, 0, sizeof(s->hist));
  return s;
}

void soxr_stats_delete(soxr_stats_t * s)
{
  if (s) {
    pthread_mutex_destroy(&s->mutex);
    free(s);
  }
}

soxr_error_t soxr_stats_process(soxr_stats_t * s,
    void const * in, size_t ilen, size_t * idone,
    void       * out, size_t olen, size_t * odone)
{
  size_t idone1 = 0, odone1 = 0, id, od;
  double t0 = now_seconds();
  soxr_error_t error = soxr_process(s->resampler,
      in, ilen, idone? idone : &idone1,
      out, olen, odone? odone : &odone1);
  double dt = now_seconds() - t0;

  id = idone? *idone : idone1;
  od = odone? *odone : odone1;
  {
    unsigned bin = hist_bin(dt);          /* Pure computation: before lock. */
    pthread_mutex_lock(&s->mutex);
    s->blocks += 1;
    s->in_samples  += (double)id;
    s->out_samples += (double)od;
    s->process_seconds += dt;
    s->hist[bin] += 1;
    if (id) {
      /* soxr_delay() is the instantaneous pipeline contents, which decays to
       * 0 as the stream is flushed; while input is being consumed it reads
       * the steady-state input->output delay, the figure of interest. */
      s->latency_out = soxr_delay(s->resampler);
      s->latency_set = 1;
    }
    pthread_mutex_unlock(&s->mutex);
  }
  return error;
}

void soxr_stats_reset(soxr_stats_t * s)
{
  pthread_mutex_lock(&s->mutex);
  s->blocks = 0;
  s->in_samples = s->out_samples = s->process_seconds = 0;
  s->latency_out = 0, s->latency_set = 0;
  memset(s->hist, 0, sizeof(s->hist));
  pthread_mutex_unlock(&s->mutex);
}

void soxr_stats_report(soxr_stats_t * s, soxr_stats_report_t * r)
{
  pthread_mutex_lock(&s->mutex);
  r->blocks          = s->blocks;
  r->in_samples      = s->in_samples;
  r->out_samples     = s->out_samples;
  r->process_seconds = s->process_seconds;
  r->latency_out_samples = s->latency_set? s->latency_out :
      soxr_delay(s->resampler);
  memcpy(r->hist, s->hist, sizeof(r->hist));
  pthread_mutex_unlock(&s->mutex);
  r->latency_in_samples  = r->latency_out_samples * s->irate / s->orate;
  r->latency_ms          = 1e3 * r->latency_out_samples / s->orate;
  r->throughput_in_sps   = r->process_seconds > 0?
      r->in_samples  / r->process_seconds : 0;
  r->throughput_out_sps  = r->process_seconds > 0?
      r->out_samples / r->process_seconds : 0;
}

void soxr_stats_hist_bin_us(unsigned bin, double * lo_us, double * hi_us)
{
  if (bin == 0) {
    *lo_us = 0, *hi_us = 1;
  } else {
    *lo_us = ldexp(1., (int)bin - 1);
    *hi_us = ldexp(1., (int)bin);
  }
}
