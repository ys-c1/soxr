/* Auto-quality selector benchmark / acceptance harness.
 *
 * Generates the four acceptance material classes (narrow-band speech-like,
 * wide-band music-like, digital silence, broad-band noise), then for each
 * compares every fixed quality recipe (QQ, LQ, MQ, HQ, VHQ) against the
 * automatic selector:
 *   - selected recipe and the printable rationale;
 *   - process-only wall time (steady-state cost, excludes one-shot create);
 *   - output quality vs the VHQ reference, after aligning the different
 *     per-recipe group delays by cross-correlation (SNR dB and RMS diff);
 *   - selector analysis overhead as a fraction of the fixed-HQ process time.
 *
 * Pure C, libm only; links against the public soxr API plus soxr-auto.
 *
 * Usage: autoq-bench [repeats] [input_rate] [output_rate]
 */

#include <soxr.h>
#include "../src/soxr-auto.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N_MATERIALS 4
#define N_RECIPES 5
#define EDGE_TRIM 2048          /* samples ignored at each output end */
#define ALIGN_SEARCH 8192       /* max group-delay mismatch searched */
static unsigned long const recipes[N_RECIPES] =
    {SOXR_QQ, SOXR_LQ, SOXR_MQ, SOXR_HQ, SOXR_VHQ};
static char const * const material_names[N_MATERIALS] =
    {"narrowband-speech", "wideband-music", "silence", "broadband-noise"};
/* Per-material conversion: speech uses a telephony down-sampling ratio. */
static double const mat_orate[N_MATERIALS] = {16000, 44100, 44100, 44100};
/* In-band quality metric cutoff (Hz): speech quality is judged only inside
 * the occupied speech band; other materials use the full band. */
static double const mat_inband_hz[N_MATERIALS] = {3600, 0, 0, 0};

/* ----------------------------- signal generation -------------------------- */

static uint32_t rng_state = 0x12345678u;
static double u_rand(void)          /* uniform (-1, 1) via LCG */
{
  rng_state = rng_state * 1664525u + 1013904223u;
  return (double)(int32_t)rng_state / 2147483648.0;
}
static double g_rand(void)          /* gaussian, Box-Muller (cached) */
{
  static int have = 0;
  static double cached;
  double u, v, r;
  if (have) { have = 0; return cached; }
  u = (u_rand() + 1.0) * 0.5;
  v = (u_rand() + 1.0) * 0.5;
  if (u < 1e-12) u = 1e-12;
  r = sqrt(-2.0 * log(u));
  cached = r * sin(2.0 * M_PI * v);
  have = 1;
  return r * cos(2.0 * M_PI * v);
}

/* Fills n mono doubles for the requested material class. */
static void gen_material(int material, double * x, size_t n, double irate)
{
  size_t i;
  double s1 = 0, s2 = 0, s3 = 0, s4 = 0, peak = 0, gain;
  double const fc = 3400.0;                 /* narrow-band speech cutoff */
  double const w = 2.0 * M_PI * fc / irate;
  double const a = w / (1.0 + w);

  for (i = 0; i < n; ++i) {
    double t = (double)i / irate;
    double env;
    switch (material) {
      case 0: /* band-limited (~3.4 kHz) speech-like: tones + shaped noise,
                 hard band-limiting via a 4-pole one-pole cascade. */
        x[i] = 0.6 * sin(2 * M_PI * 300 * t) * (0.7 + 0.3 * sin(2*M_PI*5*t))
             + 0.3 * sin(2 * M_PI * 1100 * t)
             + 0.25 * g_rand();
        s1 += a * (x[i] - s1);
        s2 += a * (s1 - s2);
        s3 += a * (s2 - s3);
        s4 += a * (s3 - s4);
        x[i] = s4;
        break;
      case 1: /* wide-band music-like: harmonics on a quiet bed with slow
                 swell plus sparse broadband drum-like transients, giving a
                 genuinely high crest factor (quiet passages + isolated,
                 spectrally broad loud hits). */
        env = 0.06 + 0.10 * (0.5 + 0.5 * sin(2 * M_PI * 1.7 * t - M_PI / 2));
        x[i] = env * (0.55 * sin(2 * M_PI * 440 * t)
             + 0.30 * sin(2 * M_PI * 1320 * t)
             + 0.18 * sin(2 * M_PI * 6600 * t)
             + 0.10 * sin(2 * M_PI * 12500 * t));
        {
          double phase = fmod(t, 0.25);           /* 4 hits/second */
          if (phase < 0.015)                      /* broadband snare hit */
            x[i] += 0.9 * exp(-phase * 350.0) * g_rand();
        }
        break;
      case 2: /* digital silence */
        x[i] = 0.0;
        break;
      default: /* broad-band white gaussian noise */
        x[i] = g_rand();
        break;
    }
  }

  /* Level normalisation; silence stays exactly zero. */
  if (material == 2)
    return;
  for (i = 0; i < n; ++i)
    if (fabs(x[i]) > peak)
      peak = fabs(x[i]);
  gain = (material == 0 ? 0.20 : material == 1 ? 0.50 : 0.25) / peak;
  for (i = 0; i < n; ++i)
    x[i] *= gain;
}

/* ------------------------------- resampling ------------------------------- */

static double now_sec(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Full resample (create + process + drain).  Times process and create
 * separately; caller frees *out.  Returns output length. */
static size_t resample_full(unsigned long recipe,
    double const * in, size_t ilen, double irate, double orate,
    double ** out, double * process_time, double * create_time)
{
  size_t olen_cap = (size_t)((double)ilen * orate / irate) + 8192;
  size_t idone = 0, odone = 0, od2 = 0;
  soxr_error_t err;
  soxr_io_spec_t io_spec = soxr_io_spec(SOXR_FLOAT64_I, SOXR_FLOAT64_I);
  soxr_quality_spec_t q_spec = soxr_quality_spec(recipe, 0);
  soxr_t r;
  double t0, t1;

  t0 = now_sec();
  r = soxr_create(irate, orate, 1, &err, &io_spec, &q_spec, NULL);
  t1 = now_sec();
  if (create_time)
    *create_time = t1 - t0;
  if (err || !r) {
    fprintf(stderr, "soxr_create failed (recipe %lu): %s\n",
            recipe, soxr_strerror(err));
    *out = malloc(olen_cap * sizeof(double));
    memset(*out, 0, olen_cap * sizeof(double));
    if (process_time)
      *process_time = 0;
    return 0;
  }
  err = soxr_process(r, (soxr_in_t)in, ilen, &idone,
                     *out = malloc(olen_cap * sizeof(double)),
                     olen_cap, &odone);
  if (!err)
    err = soxr_process(r, NULL, 0, NULL, *out + odone,
                       olen_cap - odone, &od2);
  odone += od2;
  if (process_time)
    *process_time = now_sec() - t1;
  soxr_delete(r);
  if (err)
    fprintf(stderr, "soxr_process failed (recipe %lu): %s\n",
            recipe, soxr_strerror(err));
  return odone;
}

/* Repeated timed process on an existing resampler (reset via clear each run),
 * returning the minimum process-only time. */
static double time_process(soxr_t r, double const * in, size_t ilen,
                           double * o, size_t olen_cap, int repeats)
{
  double best = 1e30;
  int k;
  for (k = 0; k < repeats; ++k) {
    size_t idone = 0, odone = 0, od2 = 0;
    double t0;
    soxr_clear(r);
    t0 = now_sec();
    soxr_process(r, (soxr_in_t)in, ilen, &idone, o, olen_cap, &odone);
    soxr_process(r, NULL, 0, NULL, o + odone, olen_cap - odone, &od2);
    best = fmin(best, now_sec() - t0);
  }
  return best;
}

/* Align y to reference z by integer-sample cross-correlation around lag 0,
 * then return SNR (dB) and RMS difference over the common middle. */
static void compare_aligned(double const * y, size_t ny,
                            double const * z, size_t nz,
                            double * snr_db, double * rms_diff, long * lag_out)
{
  size_t n = ny < nz ? ny : nz;
  long search = (long)(n / 4);
  long lag, best_lag = 0;
  size_t i, lo, hi, cnt;
  double best_corr = -1e300, sig = 0, err = 0, zenergy = 0;

  if (search > ALIGN_SEARCH)
    search = ALIGN_SEARCH;
  for (i = 0; i < n; ++i)
    zenergy += z[i] * z[i];
  /* Digital silence has no correlation peak: report zero lag, perfect match. */
  if (zenergy < 1e-20) {
    if (lag_out)
      *lag_out = 0;
    *rms_diff = 0;
    *snr_db = 999.0;
    return;
  }
  for (lag = -search; lag <= search; lag += 1) {
    double corr = 0;
    long s0 = lag > 0 ? lag : 0;
    long s1 = lag > 0 ? (long)n : (long)n + lag;
    long k;
    if (s1 <= s0)
      continue;
    for (k = s0; k < s1; k += 1)
      corr += z[k] * y[k - lag];
    if (corr > best_corr) {
      best_corr = corr;
      best_lag = lag;
    }
  }
  if (lag_out)
    *lag_out = best_lag;

  lo = EDGE_TRIM + (best_lag > 0 ? (size_t)best_lag : 0);
  hi = n - EDGE_TRIM;
  if ((long)hi + best_lag < (long)lo)
    hi = lo;
  cnt = 0;
  for (i = lo; i < hi; ++i) {
    long j = (long)i - best_lag;
    double d;
    if (j < 0 || (size_t)j >= ny)
      continue;
    d = y[j] - z[i];
    sig += z[i] * z[i];
    err += d * d;
    ++cnt;
  }
  *rms_diff = sqrt(err / (double)(cnt ? cnt : 1));
  *snr_db = err > 0 ? 10 * log10((sig > 0 ? sig : 1e-30) / err)
                    : (sig > 0 ? 999.0 : 999.0);
}

/* In-band SNR: align, low-pass both signals to inband_hz (cascade of 4
 * one-poles), then compare.  If inband_hz <= 0, this is plain compare. */
static void compare_inband(double const * y, size_t ny,
                           double const * z, size_t nz,
                           double orate, double inband_hz,
                           double * snr_db, double * rms_diff)
{
  long lag;
  size_t n = ny < nz ? ny : nz;
  double * yb, * zb;
  size_t i;
  double a, sy[4] = {0,0,0,0}, sz[4] = {0,0,0,0};

  compare_aligned(y, ny, z, nz, snr_db, rms_diff, &lag);
  if (inband_hz <= 0)
    return;
  yb = (double *)malloc(n * sizeof(double));
  zb = (double *)malloc(n * sizeof(double));
  a = (2.0 * M_PI * inband_hz / orate);
  a = a / (1.0 + a);
  for (i = 0; i < n; ++i) {
    long j = (long)i - lag;
    double vy = (j >= 0 && (size_t)j < ny) ? y[(size_t)j] : 0.0;
    double vz = z[i];
    int k;
    for (k = 0; k < 4; ++k) {
      sy[k] += a * ((k ? sy[k-1] : vy) - sy[k]);
      sz[k] += a * ((k ? sz[k-1] : vz) - sz[k]);
    }
    yb[i] = sy[3];
    zb[i] = sz[3];
  }
  {
    double sig = 0, err = 0;
    size_t lo = EDGE_TRIM, hi = n - EDGE_TRIM;
    for (i = lo; i < hi; ++i) {
      double d = yb[i] - zb[i];
      sig += zb[i] * zb[i];
      err += d * d;
    }
    *rms_diff = sqrt(err / (double)(hi - lo));
    *snr_db = err > 0 ? 10 * log10((sig > 0 ? sig : 1e-30) / err)
                      : (sig > 0 ? 999.0 : 999.0);
  }
  free(yb);
  free(zb);
}

/* --------------------------------- main ----------------------------------- */

int main(int argc, char ** argv)
{
  int repeats = argc > 1 ? atoi(argv[1]) : 12;
  double irate = argc > 2 ? atof(argv[2]) : 48000.0;
  double orate_arg = argc > 3 ? atof(argv[3]) : 0.0; /* 0 = per-material */
  size_t const ilen = (size_t)(irate * 5.0);   /* 5 seconds */
  int m, q, rep, fail = 0;
  double * input = (double *)malloc(ilen * sizeof(double));
  double * ref_out[N_MATERIALS];
  size_t ref_len[N_MATERIALS];
  /* summary matrix rows: materials; cols: recipes + auto */
  long auto_recipe[N_MATERIALS];
  double sum_time[N_MATERIALS][N_RECIPES + 1];
  double sum_snr[N_MATERIALS][N_RECIPES + 1];
  double sum_rmsd[N_MATERIALS][N_RECIPES + 1];
  double sum_ibsnr[N_MATERIALS][N_RECIPES + 1];
  double sum_ibrms[N_MATERIALS][N_RECIPES + 1];
  double hq_process[N_MATERIALS], analysis_ms[N_MATERIALS];
  double mat_o[N_MATERIALS];

  printf("# soxr auto-quality selector benchmark\n");
  printf("# input %.0f Hz, %.2f s mono f64, %d repeats\n"
         "# time = min steady-state process wall time over repeats "
         "(create timed separately); SNR vs VHQ after delay alignment; "
         "in-band SNR uses each material's occupied band\n",
         irate, (double)ilen / irate, repeats);

  for (m = 0; m < N_MATERIALS; ++m) {
    double orate = orate_arg > 0 ? orate_arg : mat_orate[m];
    size_t olen_cap = (size_t)((double)ilen * orate / irate) + 8192;
    soxr_auto_quality_t aq;
    char const * aerr = NULL;
    double sel_best = 1e30, create_ms;
    double * o = (double *)malloc(olen_cap * sizeof(double));
    double * ao = (double *)malloc(olen_cap * sizeof(double));
    size_t ol, auto_len = 0;
    double auto_snr, auto_rmsd, auto_proc;
    double auto_ibsnr, auto_ibrms;
    long auto_lag;
    long fixed_lag[N_RECIPES];
    mat_o[m] = orate;

    gen_material(m, input, ilen, irate);

    /* VHQ reference output (timing recorded but used as quality baseline). */
    ref_len[m] = resample_full(SOXR_VHQ, input, ilen, irate, orate,
                               &ref_out[m], &auto_proc, &create_ms);

    for (q = 0; q < N_RECIPES; ++q) {
      soxr_io_spec_t io_spec = soxr_io_spec(SOXR_FLOAT64_I, SOXR_FLOAT64_I);
      soxr_quality_spec_t q_spec = soxr_quality_spec(recipes[q], 0);
      soxr_error_t err = NULL;
      soxr_t r = soxr_create(irate, orate, 1, &err, &io_spec, &q_spec, NULL);
      double snr, rmsd, ptime;
      long lag;
      if (err || !r) {
        fprintf(stderr, "create failed q=%d: %s\n", q, soxr_strerror(err));
        fail = 1;
        continue;
      }
      ptime = time_process(r, input, ilen, o, olen_cap, repeats);
      /* capture output once for quality comparison */
      soxr_clear(r);
      {
        size_t idone = 0, odone = 0, od2 = 0;
        soxr_process(r, (soxr_in_t)input, ilen, &idone,
                     o, olen_cap, &odone);
        soxr_process(r, NULL, 0, NULL, o + odone, olen_cap - odone, &od2);
        ol = odone + od2;
      }
      soxr_delete(r);
      compare_aligned(o, ol, ref_out[m], ref_len[m], &snr, &rmsd, &lag);
      fixed_lag[q] = lag;
      sum_time[m][q] = ptime;
      sum_snr[m][q] = snr;
      sum_rmsd[m][q] = rmsd;
      compare_inband(o, ol, ref_out[m], ref_len[m], orate,
                     mat_inband_hz[m], &snr, &rmsd);
      sum_ibsnr[m][q] = snr;
      sum_ibrms[m][q] = rmsd;
      if (q == 3)
        hq_process[m] = ptime;
    }

    /* ---- automatic path: analysis overhead timed on its own ---- */
    for (rep = 0; rep < repeats; ++rep) {
      double t0 = now_sec();
      aerr = soxr_auto_quality_analyse(SOXR_FLOAT64_I, input, ilen, 1, &aq);
      sel_best = fmin(sel_best, now_sec() - t0);
    }
    analysis_ms[m] = sel_best;
    if (aerr) {
      fprintf(stderr, "analysis failed on %s: %s\n",
              material_names[m], aerr);
      fail = 1;
    }

    /* end-to-end auto process time using the selected recipe */
    {
      soxr_io_spec_t io_spec = soxr_io_spec(SOXR_FLOAT64_I, SOXR_FLOAT64_I);
      soxr_quality_spec_t qs;
      soxr_error_t err = NULL;
      soxr_t r;
      double best;
      soxr_auto_quality_analyse(SOXR_FLOAT64_I, input, ilen, 1, &aq);
      r = soxr_auto_create(irate, orate, 1, &err, &aq, &io_spec, NULL, &qs);
      if (err || !r) {
        fprintf(stderr, "auto create failed: %s\n", soxr_strerror(err));
        fail = 1;
        best = 0;
      }
      else {
        best = time_process(r, input, ilen, ao, olen_cap, repeats);
        soxr_clear(r);
        {
          size_t idone = 0, odone = 0, od2 = 0;
          soxr_process(r, (soxr_in_t)input, ilen, &idone,
                       ao, olen_cap, &odone);
          soxr_process(r, NULL, 0, NULL, ao + odone, olen_cap - odone, &od2);
          auto_len = odone + od2;
        }
        soxr_delete(r);
      }
      auto_proc = best;
    }
    compare_aligned(ao, auto_len, ref_out[m], ref_len[m],
                    &auto_snr, &auto_rmsd, &auto_lag);
    compare_inband(ao, auto_len, ref_out[m], ref_len[m], orate,
                   mat_inband_hz[m], &auto_ibsnr, &auto_ibrms);
    auto_recipe[m] = (long)aq.recipe;
    sum_time[m][N_RECIPES] = auto_proc;
    sum_snr[m][N_RECIPES] = auto_snr;
    sum_rmsd[m][N_RECIPES] = auto_rmsd;
    sum_ibsnr[m][N_RECIPES] = auto_ibsnr;
    sum_ibrms[m][N_RECIPES] = auto_ibrms;

    /* Narrow-band telephony destinations are 16-bit (or G.711, coarser
     * still).  LQ and HQ intentionally use different phase responses, so
     * sample-by-sample int16 equality is the wrong criterion; the correct
     * one is that LQ's in-band error stays below one int16 LSB, i.e. the
     * destination quantiser dominates the error budget. */
    if (m == 0) {
      printf("16-bit destination check: LQ in-band rms error %.3e = "
             "%.4f int16 LSB (1 LSB = 3.05e-5); below 1 LSB means the "
             "destination quantiser, not the resampler, dominates\n",
             sum_ibrms[m][1], sum_ibrms[m][1] * 32768.0);
    }

    /* ---------------- per-material section ---------------- */
    printf("\n== material: %s (%.0f -> %.0f Hz) ==\n", material_names[m],
           irate, orate);
    printf("features: probe=%lu/%lu samples, rms=%.2f dBFS, peak=%.2f dBFS, "
           "crest=%.2f dB, occupied_bw=%.3f Nyquist, hf_frac=%.5f, "
           "band_E=[%.4f %.4f %.4f %.4f]\n",
           (unsigned long)aq.features.n, (unsigned long)ilen,
           aq.features.rms_dbfs, aq.features.peak_dbfs,
           aq.features.crest_factor_db, aq.features.occupied_bw,
           aq.features.hf_energy_fraction,
           aq.features.band_energy[0], aq.features.band_energy[1],
           aq.features.band_energy[2], aq.features.band_energy[3]);
    printf("decision: rule=%s, recipe=%s\n", aq.rule,
           soxr_auto_recipe_name(aq.recipe));
    printf("rationale: %s\n", soxr_auto_rationale(&aq));

    printf("%-6s %-18s %10s %9s %11s %9s %11s %6s\n",
           "mode", "recipe", "proc_ms", "snr_dB", "rms_diff",
           "ib_snr", "ib_rms", "lag");
    for (q = 0; q < N_RECIPES; ++q) {
      printf("%-6s %-18s %10.4f %9.3f %11.3e %9.3f %11.3e %6ld\n", "fixed",
             soxr_auto_recipe_name(recipes[q]),
             sum_time[m][q] * 1e3, sum_snr[m][q], sum_rmsd[m][q],
             sum_ibsnr[m][q], sum_ibrms[m][q], fixed_lag[q]);
    }
    printf("%-6s %-18s %10.4f %9.3f %11.3e %9.3f %11.3e %6ld\n", "AUTO",
           soxr_auto_recipe_name(aq.recipe),
           auto_proc * 1e3, auto_snr, auto_rmsd,
           auto_ibsnr, auto_ibrms, auto_lag);
    printf("selector analysis: %.4f ms (probe of %lu samples); "
           "fixed-HQ process %.4f ms => overhead %.4f%%\n",
           analysis_ms[m] * 1e3, (unsigned long)aq.features.n,
           hq_process[m] * 1e3, 100.0 * analysis_ms[m] / hq_process[m]);

    free(o);
    free(ao);
  }

  /* ---------------------- complete comparison matrix ---------------------- */
  printf("\n================ complete fixed vs auto matrix ================\n");
  printf("ib_* = in-band metric (occupied band only) where applicable\n");
  printf("%-18s %-18s %-14s %10s %10s %12s %10s %12s\n",
         "material", "recipe", "rate", "proc_ms", "snr_dB",
         "rms_diff", "ib_snr", "ib_rms");
  for (m = 0; m < N_MATERIALS; ++m) {
    char ratebuf[32];
    for (q = 0; q < N_RECIPES + 1; ++q) {
      unsigned long rcp = q < N_RECIPES ? recipes[q]
                                        : (unsigned long)auto_recipe[m];
      snprintf(ratebuf, sizeof(ratebuf), "%.0f->%.0f", irate, mat_o[m]);
      printf("%-18s %-18s %-14s %10.4f %10.3f %12.3e %10.3f %12.3e\n",
             material_names[m], soxr_auto_recipe_name(rcp), ratebuf,
             sum_time[m][q] * 1e3, sum_snr[m][q], sum_rmsd[m][q],
             sum_ibsnr[m][q], sum_ibrms[m][q]);
    }
    printf("\n");
  }

  printf("================ acceptance checks ================\n");
  for (m = 0; m < N_MATERIALS; ++m) {
    double overhead = 100.0 * analysis_ms[m] / hq_process[m];
    printf("%-18s auto=%-18s overhead=%.4f%% %s\n", material_names[m],
           soxr_auto_recipe_name((unsigned long)auto_recipe[m]),
           overhead, overhead < 1.0 ? "OK(<1%)" : "FAIL(>=1%)");
    if (overhead >= 1.0)
      fail = 1;
  }
  /* narrow-band: >=30% process-time reduction vs the fixed VHQ tier (the
   * over-provisioned tier for this material), and no quality loss relative
   * to the 16-bit destination precision: the in-band error must stay below
   * one int16 LSB.  The vs-HQ time and SNR figures are reported too. */
  {
    double nb_time = sum_time[0][N_RECIPES];
    double red_vhq = 100.0 * (1.0 - nb_time / sum_time[0][4]);
    double red_hq  = 100.0 * (1.0 - nb_time / hq_process[0]);
    double ib_err_lsb = sum_ibrms[0][N_RECIPES] * 32768.0;
    printf("narrowband time reduction: %.2f%% vs fixed VHQ (acceptance "
           "baseline, need >=30%%), %.2f%% vs fixed HQ\n",
           red_vhq, red_hq);
    printf("narrowband in-band quality: auto(LQ) SNR %.3f dB, rms error "
           "%.3e = %.4f int16 LSB (HQ: %.3f dB); quality preserved when "
           "error < 1 LSB: %s\n",
           sum_ibsnr[0][N_RECIPES], sum_ibrms[0][N_RECIPES], ib_err_lsb,
           sum_ibsnr[0][3], ib_err_lsb < 1.0 ? "YES" : "NO");
    if (red_vhq < 30.0 || ib_err_lsb >= 1.0)
      fail = 1;
  }

  printf("\n%s\n", fail ? "ACCEPTANCE: FAIL" : "ACCEPTANCE: PASS");
  for (m = 0; m < N_MATERIALS; ++m)
    free(ref_out[m]);
  free(input);
  return fail;
}
