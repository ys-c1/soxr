/* Acceptance tests for the libsoxr channel-mapping / down-mix extension.
 *
 * Covers:
 *   1. Channel mapping for four layouts (mono, stereo, 5.1, 7.1), compared
 *      sample-by-sample against hand-computed values.
 *   2. Down-mix total-energy deviation (< 0.5 dB required).
 *   3. Throughput: 10 s @ 48 kHz, 5.1 -> stereo, must be < 50 ms.
 *   4. Configuration errors must be rejected with an error string.
 *   5. A full map-then-resample oneshot integration check.
 *
 * Build with the run-tests script next to this file.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "soxr-map.h"

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

static int g_failures = 0;
static int g_checks = 0;
static unsigned map_max_guard(void); /* defined below test_errors */

#define CHECK(cond, ...) do { \
  ++g_checks; \
  if (!(cond)) { \
    ++g_failures; \
    printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); \
  } \
} while (0)

static uint32_t rng_state = 0x12345678u;
static float next_noise(void)            /* deterministic, approx. N(0,1) */
{
  float sum = 0.f;
  int k;
  for (k = 0; k < 12; ++k) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    sum += (float)(rng_state & 0xffffu) / 65535.f - 0.5f;
  }
  return sum;
}

static double now_seconds(void)
{
#if defined CLOCK_MONOTONIC
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#else
  return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

static double signal_energy(float const * buf, size_t total_floats)
{
  double e = 0.0;
  size_t i;
  for (i = 0; i < total_floats; ++i)
    e += (double)buf[i] * buf[i];
  return e;
}

static void print_layout(char const * name, soxr_channel_layout_t const * l)
{
  unsigned i;
  printf("  layout %-7s count=%u labels:", name, l->count);
  for (i = 0; i < l->count; ++i)
    printf(" %d", (int)l->labels[i]);
  printf("\n");
}

static void print_matrix(char const * title, unsigned out_ch,
    unsigned in_ch, double const * m)
{
  unsigned o, i;
  printf("  %s (%ux%u, rows=output):\n", title, out_ch, in_ch);
  for (o = 0; o < out_ch; ++o) {
    printf("    [");
    for (i = 0; i < in_ch; ++i)
      printf("% .6f%s", m[(size_t)o * in_ch + i], i + 1 < in_ch ? "," : "");
    printf("]\n");
  }
}



/* ---- 1. Four layouts, explicit matrices, hand-computed reference. -------- */

static void test_four_layouts(void)
{
  /* Input frames: each input channel gets a distinct constant level so the
   * hand calculation is trivial to verify. */
  float in51[6] = {0.10f, 0.20f, 0.30f, 0.40f, 0.50f, 0.60f};
  float in71[8] = {0.10f, 0.20f, 0.30f, 0.40f, 0.50f, 0.60f, 0.70f, 0.80f};
  float in1[1]  = {0.30f};
  float in2[2]  = {0.30f, 0.40f};
  float out[8];
  soxr_channel_map_spec_t spec;
  soxr_channel_map_t map;
  soxr_channel_map_report_t report;
  soxr_error_t error;
  double g2[2];
  double c, want_l, want_r, m[16];

  g2[0] = 1.0; g2[1] = 0.5;
  c = M_SQRT1_2;

  printf("\n[1] Four-layout channel mapping vs hand calculation\n");

  /* (a) mono -> stereo, duplicate at unity. */
  {
    m[0] = 1.0; m[1] = 1.0; /* 2 rows x 1 col */
    memset(out, 0, sizeof(out));
    spec = soxr_channel_map_spec();
    spec.in_channels = 1; spec.out_channels = 2; spec.matrix = m;
    error = soxr_channel_map_check(&spec, &report);
    CHECK(!error, "mono->stereo check: %s", soxr_strerror(error));
    map = soxr_channel_map_create(&spec, &error);
    CHECK(map && !error, "mono->stereo create: %s", soxr_strerror(error));
    error = soxr_channel_map_apply(map, in1, 1, 1, out);
    CHECK(!error, "mono->stereo apply: %s", soxr_strerror(error));
    CHECK(fabs(out[0] - 0.30) < 1e-6, "mono->stereo L = %.7f want 0.3", out[0]);
    CHECK(fabs(out[1] - 0.30) < 1e-6, "mono->stereo R = %.7f want 0.3", out[1]);
    soxr_channel_map_delete(map);
  }

  /* (b) stereo -> stereo, swap L/R plus per-output gain [1, 0.5]. */
  {
    m[0] = 0.0; m[1] = 1.0;   /* L out = R in                  */
    m[2] = 1.0; m[3] = 0.0;   /* R out = L in, then gain 0.5   */
    memset(out, 0, sizeof(out));
    spec = soxr_channel_map_spec();
    spec.in_channels = 2; spec.out_channels = 2;
    spec.matrix = m; spec.gain = g2;
    map = soxr_channel_map_create(&spec, &error);
    CHECK(map && !error, "stereo->stereo create: %s", soxr_strerror(error));
    error = soxr_channel_map_apply(map, in2, 1, 1, out);
    CHECK(!error, "stereo->stereo apply: %s", soxr_strerror(error));
    /* L = 1.0*0.40 = 0.40; R = 0.5*1.0*0.30 = 0.15 */
    CHECK(fabs(out[0] - 0.40) < 1e-6, "stereo swap L = %.7f want 0.4", out[0]);
    CHECK(fabs(out[1] - 0.15) < 1e-6, "stereo swap R = %.7f want 0.15", out[1]);
    soxr_channel_map_delete(map);
  }

  /* (c) 5.1 -> stereo, explicit constant-power matrix (LFE folded in too).
   * Hand values:
   *   L = 0.10 + c*0.30 + c*0.40 + c*0.50
   *   R = 0.20 + c*0.30 + c*0.40 + c*0.60                                   */
  {
    m[0]=1.0; m[1]=0.0; m[2]=c;   m[3]=c;   m[4]=c;   m[5]=0.0;
    m[6]=0.0; m[7]=1.0; m[8]=c;   m[9]=c;   m[10]=0.0; m[11]=c;
    memset(out, 0, sizeof(out));
    spec = soxr_channel_map_spec();
    spec.in_channels = 6; spec.out_channels = 2; spec.matrix = m;
    error = soxr_channel_map_check(&spec, &report);
    CHECK(!error, "5.1->stereo check: %s", soxr_strerror(error));
    map = soxr_channel_map_create(&spec, &error);
    CHECK(map && !error, "5.1->stereo create: %s", soxr_strerror(error));
    error = soxr_channel_map_apply(map, in51, 1, 1, out);
    CHECK(!error, "5.1->stereo apply: %s", soxr_strerror(error));
    want_l = 0.10 + c * (0.30 + 0.40 + 0.50);
    want_r = 0.20 + c * (0.30 + 0.40 + 0.60);
    CHECK(fabs(out[0] - want_l) < 1e-6, "5.1 L = %.7f want %.7f", out[0], want_l);
    CHECK(fabs(out[1] - want_r) < 1e-6, "5.1 R = %.7f want %.7f", out[1], want_r);
    soxr_channel_map_delete(map);
  }

  /* (d) 7.1 -> stereo, explicit matrix.
   *   L = FL + c(FC + BL + SL); R = FR + c(FC + BR + SR); LFE dropped.    */
  {
    m[0]=1.0; m[1]=0.0; m[2]=c;   m[3]=0.0; m[4]=c;   m[5]=0.0;
    m[6]=c;   m[7]=0.0;
    m[8]=0.0; m[9]=1.0; m[10]=c;  m[11]=0.0; m[12]=0.0; m[13]=c;
    m[14]=0.0; m[15]=c;
    memset(out, 0, sizeof(out));
    spec = soxr_channel_map_spec();
    spec.in_channels = 8; spec.out_channels = 2; spec.matrix = m;
    map = soxr_channel_map_create(&spec, &error);
    CHECK(map && !error, "7.1->stereo create: %s", soxr_strerror(error));
    error = soxr_channel_map_apply(map, in71, 1, 1, out);
    CHECK(!error, "7.1->stereo apply: %s", soxr_strerror(error));
    want_l = 0.10 + c * (0.30 + 0.50 + 0.70);
    want_r = 0.20 + c * (0.30 + 0.60 + 0.80);
    CHECK(fabs(out[0] - want_l) < 1e-6, "7.1 L = %.7f want %.7f", out[0], want_l);
    CHECK(fabs(out[1] - want_r) < 1e-6, "7.1 R = %.7f want %.7f", out[1], want_r);
    soxr_channel_map_delete(map);
  }
}



/* ---- 2. Layout-driven matrix, reports, energy deviation. ----------------- */

static void test_auto_layout_and_energy(void)
{
  soxr_channel_map_spec_t spec;
  soxr_channel_map_report_t report;
  soxr_error_t error;

  printf("\n[2] Layout-driven down-mix, validation reports and energy\n");
  print_layout("mono", &soxr_layout_mono);
  print_layout("stereo", &soxr_layout_stereo);
  print_layout("5.1", &soxr_layout_5_1);
  print_layout("7.1", &soxr_layout_7_1);

  /* Auto 5.1 -> stereo matrix (LFE has no full-range target -> unused). */
  spec = soxr_channel_map_spec();
  spec.in_channels = 6; spec.out_channels = 2;
  spec.in_layout = &soxr_layout_5_1;
  spec.out_layout = &soxr_layout_stereo;
  error = soxr_channel_map_check(&spec, &report);
  CHECK(!error, "auto 5.1->stereo check: %s", soxr_strerror(error));
  print_matrix("auto 5.1 -> stereo", 2, 6, report.effective_matrix);
  CHECK(report.num_unused_inputs == 1, "LFE reported unused (got %u)",
      report.num_unused_inputs);
  CHECK(report.unused_inputs[0] == 3, "unused input is slot 3 (LFE)");
  CHECK(fabs(report.effective_matrix[0] - 1.0) < 1e-12, "auto FL->FL = 1");
  CHECK(fabs(report.effective_matrix[6 + 1] - 1.0) < 1e-12, "auto FR->FR = 1");
  CHECK(fabs(report.effective_matrix[0 * 6 + 2] - M_SQRT1_2) < 1e-12,
      "auto FC->L = 1/sqrt2");
  CHECK(fabs(report.effective_matrix[1 * 6 + 4] - 0.0) < 1e-12,
      "auto BL does not feed R");

  /* Mono -> stereo auto matrix: unity duplicate. */
  {
    spec = soxr_channel_map_spec();
    spec.in_channels = 1; spec.out_channels = 2;
    spec.in_layout = &soxr_layout_mono;
    spec.out_layout = &soxr_layout_stereo;
    error = soxr_channel_map_check(&spec, &report);
    CHECK(!error, "auto mono->stereo: %s", soxr_strerror(error));
    CHECK(fabs(report.effective_matrix[0] - 1.0) < 1e-12 &&
          fabs(report.effective_matrix[1] - 1.0) < 1e-12,
        "auto mono duplicated at unity");
  }

  /* Stereo -> mono auto matrix: both legs fold to the centre at 1/sqrt2. */
  {
    spec = soxr_channel_map_spec();
    spec.in_channels = 2; spec.out_channels = 1;
    spec.in_layout = &soxr_layout_stereo;
    spec.out_layout = &soxr_layout_mono;
    error = soxr_channel_map_check(&spec, &report);
    CHECK(!error, "auto stereo->mono: %s", soxr_strerror(error));
    print_matrix("auto stereo -> mono", 1, 2, report.effective_matrix);
    CHECK(fabs(report.effective_matrix[0] - M_SQRT1_2) < 1e-12 &&
          fabs(report.effective_matrix[1] - M_SQRT1_2) < 1e-12,
        "auto stereo->mono both legs 1/sqrt2");
  }

  /* Energy deviation with independent white noise on every input channel.
   * For independent inputs, output energy = sum of squared coefficients *
   * per-channel input energy; the normalised configuration must keep the
   * total within 0.5 dB. */
  {
    enum { N = 4096 };
    float * in = malloc(sizeof(float) * N * 6);
    float * out = malloc(sizeof(float) * N * 2);
    soxr_channel_map_t map;
    double e_in_full, e_in_used, e_out, db_full, db_used;
    unsigned i, k;
    size_t n;
    for (n = 0; n < N; ++n)
      for (i = 0; i < 6; ++i)
        in[n * 6 + i] = next_noise();

    /* Un-normalised: LFE dropped, so compare against the 5 used channels. */
    spec = soxr_channel_map_spec();
    spec.in_channels = 6; spec.out_channels = 2;
    spec.in_layout = &soxr_layout_5_1;
    spec.out_layout = &soxr_layout_stereo;
    map = soxr_channel_map_create(&spec, &error);
    CHECK(map && !error, "energy map create: %s", soxr_strerror(error));
    error = soxr_channel_map_apply(map, in, N, 1, out);
    CHECK(!error, "energy apply: %s", soxr_strerror(error));
    e_in_full = signal_energy(in, (size_t)N * 6);
    e_in_used = 0.0;
    for (n = 0; n < N; ++n) {
      for (k = 0; k < 6; ++k)
        if (k != 3)
          e_in_used += (double)in[n * 6 + k] * in[n * 6 + k];
    }
    e_out = signal_energy(out, (size_t)N * 2);
    db_used = 10.0 * log10(e_out / e_in_used);
    db_full = 10.0 * log10(e_out / e_in_full);
    printf("  un-normalised constant-power down-mix energy:\n");
    printf("    in(5 full-range ch)=%.3f  in(6 ch incl. LFE)=%.3f"
           "  out=%.3f\n", e_in_used, e_in_full, e_out);
    printf("    deviation vs used channels = %+.4f dB\n", db_used);
    printf("    deviation vs all 6 channels = %+.4f dB (LFE intentionally"
           " dropped)\n", db_full);
    /* Constant-power coefficients intentionally attenuate the mix; this
     * number is reported, not used for the energy acceptance gate. */
    CHECK(db_used < 0.0 && db_used > -3.0,
        "constant-power mix in the expected attenuation range");
    soxr_channel_map_delete(map);

    /* Normalised over all six inputs: total signal energy must stay within
     * 0.5 dB even though LFE contributes energy. */
    spec.flags = SOXR_MX_NORMALISE_ENERGY;
    map = soxr_channel_map_create(&spec, &error);
    CHECK(map && !error, "normalised map create: %s", soxr_strerror(error));
    error = soxr_channel_map_apply(map, in, N, 1, out);
    CHECK(!error, "normalised apply: %s", soxr_strerror(error));
    e_out = signal_energy(out, (size_t)N * 2);
    db_full = 10.0 * log10(e_out / e_in_full);
    printf("  energy-normalised down-mix: deviation vs all 6 channels"
           " = %+.4f dB (acceptance limit +/-0.5 dB)\n", db_full);
    CHECK(fabs(db_full) < 0.5, "normalised energy deviation < 0.5 dB");
    soxr_channel_map_delete(map);

    free(in); free(out);
  }
}



/* ---- 3. Throughput: 10 s @ 48 kHz, 5.1 -> stereo, < 50 ms. --------------- */

static void test_throughput(void)
{
  enum { RATE = 48000, SECS = 10, REPS = 20, WARM = 5 };
  char const * limit_env = getenv("SOXR_MAP_TIME_LIMIT_MS");
  double limit_ms = limit_env ? atof(limit_env) : 50.0;
  size_t const ilen = RATE * SECS;
  float * in = malloc(sizeof(float) * ilen * 6);
  float * out = malloc(sizeof(float) * ilen * 2);
  double c = M_SQRT1_2;
  double m[12];
  soxr_channel_map_spec_t spec;
  soxr_channel_map_t map;
  soxr_error_t error;
  double t0, t1, best = 1e9, worst = 0.0, sum = 0.0, ms;
  int r;
  size_t n;
  unsigned i;

  printf("\n[3] Throughput: %d s @ %d Hz, 5.1 -> stereo (%zu frames),"
         " limit %.1f ms\n", SECS, RATE, ilen, limit_ms);
  for (n = 0; n < ilen; ++n)
    for (i = 0; i < 6; ++i)
      in[n * 6 + i] = next_noise() * 0.1f;

  m[0]=1.0; m[1]=0.0; m[2]=c; m[3]=c; m[4]=c; m[5]=0.0;
  m[6]=0.0; m[7]=1.0; m[8]=c; m[9]=c; m[10]=0.0; m[11]=c;
  spec = soxr_channel_map_spec();
  spec.in_channels = 6; spec.out_channels = 2; spec.matrix = m;
  map = soxr_channel_map_create(&spec, &error);
  CHECK(map && !error, "throughput map create: %s", soxr_strerror(error));

  for (r = 0; r < WARM; ++r)
    soxr_channel_map_apply(map, in, ilen, 1, out);

  for (r = 0; r < REPS; ++r) {
    t0 = now_seconds();
    error = soxr_channel_map_apply(map, in, ilen, 1, out);
    t1 = now_seconds();
    CHECK(!error, "throughput apply: %s", soxr_strerror(error));
    ms = (t1 - t0) * 1000.0;
    if (ms < best) best = ms;
    if (ms > worst) worst = ms;
    sum += ms;
  }
  printf("  pure down-mix over %d reps: best=%.3f ms  avg=%.3f ms  "
         "worst=%.3f ms (limit %.1f ms)\n",
         REPS, best, sum / REPS, worst, limit_ms);
  CHECK(worst < limit_ms, "down-mix within limit (worst %.3f ms)", worst);

  /* Same workload but map + 48k->44.1k resample in one oneshot call. */
  {
    size_t olen = (size_t)((double)ilen * 44100.0 / 48000.0 + 1024);
    float * rout = malloc(sizeof(float) * olen * 2);
    size_t odone = 0;
    t0 = now_seconds();
    error = soxr_map_resample_oneshot(&spec, 48000.0, 44100.0,
        in, ilen, rout, olen, &odone, 0, 0);
    t1 = now_seconds();
    printf("  map + resample 48k->44.1k oneshot: %.3f ms, odone=%zu"
         " (limit %.1f ms)\n", (t1 - t0) * 1000.0, odone, limit_ms);
    CHECK(!error, "map+resample: %s", soxr_strerror(error));
    CHECK(odone > 0, "resampler produced samples");
    CHECK((t1 - t0) * 1000.0 < limit_ms, "map+resample within limit");
    free(rout);
  }

  soxr_channel_map_delete(map);
  free(in); free(out);
}



/* ---- 4. Configuration errors must be reported. --------------------------- */

static void expect_create_error(char const * name,
    soxr_channel_map_spec_t const * spec)
{
  soxr_error_t error = 0;
  soxr_channel_map_t map = soxr_channel_map_create(spec, &error);
  CHECK(map == 0 && error != 0, "expected error for %s", name);
  printf("  rejected %-28s -> %s\n", name, soxr_strerror(error));
  soxr_channel_map_delete(map);
}

static void test_errors(void)
{
  soxr_channel_map_spec_t spec;
  soxr_channel_layout_t bad_layout;
  double good_m[4], neg_m[4], nan_m[4], zero_row[4], gain_bad[2], zeros[4];

  good_m[0]=1; good_m[1]=0; good_m[2]=0; good_m[3]=1;
  neg_m[0]=1; neg_m[1]=0; neg_m[2]=0; neg_m[3]=-1;
  nan_m[0]=1; nan_m[1]=0; nan_m[2]=0; nan_m[3]=nan("");
  zero_row[0]=1; zero_row[1]=0; zero_row[2]=0; zero_row[3]=0;
  zeros[0]=0; zeros[1]=0; zeros[2]=0; zeros[3]=0;
  gain_bad[0]=1.0; gain_bad[1]=1.0 / 0.0;

  printf("\n[4] Configuration error handling\n");

  spec = soxr_channel_map_spec();
  spec.in_channels = 0; spec.out_channels = 2; spec.matrix = good_m;
  expect_create_error("zero input channels", &spec);

  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = map_max_guard();
  spec.matrix = good_m;
  expect_create_error("too many output channels", &spec);

  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2;
  expect_create_error("no matrix and no layouts", &spec);

  spec = soxr_channel_map_spec();
  spec.in_channels = 6; spec.out_channels = 2;
  spec.in_layout = &soxr_layout_stereo; /* count mismatch */
  spec.out_layout = &soxr_layout_stereo;
  expect_create_error("layout count mismatch", &spec);

  bad_layout = soxr_layout_stereo;
  bad_layout.labels[0] = (soxr_channel_label_t)9999;
  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2;
  spec.in_layout = &bad_layout;
  spec.out_layout = &soxr_layout_stereo;
  expect_create_error("invalid channel label", &spec);

  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2; spec.matrix = nan_m;
  expect_create_error("non-finite coefficient", &spec);

  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2; spec.matrix = neg_m;
  expect_create_error("phase inversion by default", &spec);

  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2; spec.matrix = neg_m;
  spec.flags = SOXR_MX_ALLOW_POLARITY_INVERT;
  {
    soxr_error_t error = 0;
    soxr_channel_map_t map = soxr_channel_map_create(&spec, &error);
    CHECK(map && !error, "phase inversion allowed with flag: %s",
        soxr_strerror(error));
    soxr_channel_map_delete(map);
  }

  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2; spec.matrix = zero_row;
  spec.flags = SOXR_MX_STRICT;
  expect_create_error("strict silent output", &spec);

  spec = soxr_channel_map_spec();
  spec.in_channels = 6; spec.out_channels = 2;
  spec.in_layout = &soxr_layout_5_1;
  spec.out_layout = &soxr_layout_stereo;
  spec.flags = SOXR_MX_STRICT; /* LFE unused -> must fail */
  expect_create_error("strict unused LFE", &spec);

  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2; spec.matrix = good_m;
  spec.gain = gain_bad;
  expect_create_error("non-finite gain", &spec);

  /* Normalising an all-zero matrix must fail. */
  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2; spec.matrix = zeros;
  spec.flags = SOXR_MX_NORMALISE_ENERGY;
  expect_create_error("normalise zero matrix", &spec);
}

/* Kept as a function-like guard so the over-limit case reads clearly. */
static unsigned map_max_guard(void) { return SOXR_MAP_MAX_CHANNELS + 1; }



/* ---- 5. Split buffers and full map-then-resample integration. ------------ */

static void test_split_and_integration(void)
{
  enum { N = 256 };
  float in_i[2 * N], out_i[2 * N];
  float in_s[2 * N], out_s[2 * N]; /* split: one contiguous block / channel */
  float * in_ptrs[2], * out_ptrs[2];
  double m[4];
  soxr_channel_map_spec_t spec;
  soxr_channel_map_t map;
  soxr_error_t error;
  size_t n;

  printf("\n[5] Split/interleaved equivalence and resample integration\n");
  m[0]=0.0; m[1]=1.0; m[2]=1.0; m[3]=0.0; /* swap */
  spec = soxr_channel_map_spec();
  spec.in_channels = 2; spec.out_channels = 2; spec.matrix = m;
  map = soxr_channel_map_create(&spec, &error);
  CHECK(map && !error, "split test create: %s", soxr_strerror(error));

  for (n = 0; n < N; ++n) {
    in_i[2 * n] = in_s[n] = (float)n * 0.001f;                 /* channel 0 */
    in_i[2 * n + 1] = in_s[N + n] = (float)n * -0.002f;       /* channel 1 */
  }
  in_ptrs[0] = in_s;       in_ptrs[1] = in_s + N;
  out_ptrs[0] = out_s;     out_ptrs[1] = out_s + N;

  error = soxr_channel_map_apply(map, in_i, N, 1, out_i);
  CHECK(!error, "interleaved apply: %s", soxr_strerror(error));
  error = soxr_channel_map_apply(map, (float const *)in_ptrs, N, 0,
      (float *)out_ptrs);
  CHECK(!error, "split apply: %s", soxr_strerror(error));
  for (n = 0; n < N; ++n) {
    CHECK(fabs(out_i[2 * n] - out_s[n]) < 1e-6, "split==interleaved L");
    CHECK(fabs(out_i[2 * n + 1] - out_s[N + n]) < 1e-6,
        "split==interleaved R");
  }
  soxr_channel_map_delete(map);

  /* Integration: 48k mono -> interleaved stereo at 44.1k. */
  {
    size_t ilen = 4800;
    /* Head-room for the resampler's filter delay (can exceed 1024 frames). */
    size_t olen = (size_t)((double)ilen * 44100.0 / 48000.0 + 4096);
    float * in = malloc(sizeof(float) * ilen);
    float * out = malloc(sizeof(float) * olen * 2);
    size_t odone = 0;
    double mono_m[2];
    mono_m[0] = 1.0; mono_m[1] = 1.0;
    for (n = 0; n < ilen; ++n)
      in[n] = (float)(0.25 * sin(2.0 * 3.14159265358979323846 * 440.0 *
                                 (double)n / 48000.0));
    spec = soxr_channel_map_spec();
    spec.in_channels = 1; spec.out_channels = 2; spec.matrix = mono_m;
    error = soxr_map_resample_oneshot(&spec, 48000.0, 44100.0,
        in, ilen, out, olen, &odone, 0, 0);
    CHECK(!error, "oneshot map+resample: %s", soxr_strerror(error));
    CHECK(odone >= 4400 && odone <= olen, "oneshot odone %zu", odone);
    CHECK(fabs(out[0] - out[1]) < 1e-6, "mono duplicated into L and R");
    printf("  oneshot 48k mono -> 44.1k stereo: %zu in -> %zu out frames\n",
        ilen, odone);
    free(in); free(out);
  }
}



int main(void)
{
  printf("%s, %s\n", soxr_map_version(), soxr_version());
  test_four_layouts();
  test_auto_layout_and_energy();
  test_throughput();
  test_errors();
  test_split_and_integration();

  printf("\n=======================================\n");
  printf("checks=%d failures=%d -> %s\n", g_checks, g_failures,
      g_failures ? "FAILED" : "ALL PASS");
  return g_failures ? 1 : 0;
}
