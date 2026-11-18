/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */



/* Auto-quality selector implementation.
 *
 * Design constraints (see soxr-auto.h):
 *   - sits purely on top of the public soxr API; the resampling kernel is
 *     untouched and no external dependency is added;
 *   - analysis runs over a bounded PROBE window (first N samples) so its cost
 *     is independent of stream length: 3 band edges, each a cascade of two
 *     first-order IIR low-pass stages (-12 dB/octave), partition the spectrum
 *     into 4 bands, alongside running RMS and peak accumulators.  No FFT, no
 *     heap allocation for <= 8 channels;
 *   - every decision threshold is named below and echoed in the rationale.
 */

#include "soxr-auto.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ----------------------------- Tuning knobs ------------------------------- */

/* Number of samples (per channel) inspected; cost is O(PROBE).
 * Overridable at compile time for latency- or accuracy-sensitive builds. */
#ifndef SOXR_AUTO_PROBE_SAMPLES
#define SOXR_AUTO_PROBE_SAMPLES 512
#endif

/* RMS below this is treated as digital/near silence. */
#define SILENCE_RMS_DBFS      (-70.0)

/* Per-band energy fraction at/above which a band counts as occupied.
 * Filter transition leakage stays below this; genuine band content does not. */
#define BAND_ACTIVE_FRACTION  (0.03)

/* Energy above 0.5 Nyquist above this fraction + low crest => noise-like. */
#define NOISE_HF_FRACTION     (0.30)
#define NOISE_MAX_CREST_DB    (14.0)

/* Highest active band edge at/below this Nyquist fraction => narrow-band. */
#define NARROWBAND_LIMIT      (0.25)

/* Crest factor at/above this => high dynamic-range material. */
#define HIGH_DYNAMIC_CREST_DB (15.0)

/* Band edges as fractions of input Nyquist; 4 bands sit between them. */
static double const band_edge[3] = {0.125, 0.25, 0.50};
static double const band_top[4]  = {0.125, 0.25, 0.50, 1.0};

/* Samples skipped at the probe start to let filter transients settle. */
#define TRANSIENT_SKIP        16



char const * soxr_auto_recipe_name(unsigned long recipe)
{
  switch (recipe) {
    case SOXR_QQ:  return "QQ (quick cubic)";
    case SOXR_LQ:  return "LQ (16-bit, fast)";
    case SOXR_MQ:  return "MQ (16-bit, medium)";
    case SOXR_HQ:  return "HQ (20-bit)";
    case SOXR_VHQ: return "VHQ (28-bit)";
    default:       return "unknown";
  }
}



/* Read one interleaved sample as float in [-1, 1]. */
static float read_sample(soxr_datatype_t t, void const * base, size_t idx)
{
  switch (t) {
    case SOXR_FLOAT32_I: return ((float const *)base)[idx];
    case SOXR_FLOAT64_I: return (float)((double const *)base)[idx];
    case SOXR_INT32_I:   return (float)(((int32_t const *)base)[idx]
                                * (1.0 / 2147483648.0));
    case SOXR_INT16_I:   return (float)(((int16_t const *)base)[idx]
                                * (1.0 / 32768.0));
    default:             return 0.f;
  }
}



unsigned long soxr_auto_quality_select(
    soxr_auto_features_t const * f,
    char * rule, size_t rule_size,
    char * reason, size_t reason_size)
{
  unsigned long recipe;
  char const * why;
  char const * id;

  if (f->rms_dbfs <= SILENCE_RMS_DBFS) {
    recipe = SOXR_QQ;
    id = "silence";
    why = "near-silence: nothing above the noise floor, QQ is bit-exact on "
          "zeros and avoids filter setup cost";
  }
  else if (f->hf_energy_fraction >= NOISE_HF_FRACTION &&
           f->crest_factor_db <= NOISE_MAX_CREST_DB) {
    recipe = SOXR_MQ;
    id = "broadband-noise";
    why = "broad-band noise-like input: energy spread to Nyquist with low "
          "crest factor, 16-bit MQ matches the achievable fidelity";
  }
  else if (f->occupied_bw <= NARROWBAND_LIMIT) {
    recipe = SOXR_LQ;
    id = "narrowband";
    why = "narrow-band content: no significant energy above 0.25 Nyquist; "
          "LQ's 16-bit error floor is below the in-band signal level and its "
          "wide transition band only attenuates unused spectrum";
  }
  else if (f->crest_factor_db >= HIGH_DYNAMIC_CREST_DB) {
    recipe = SOXR_VHQ;
    id = "wideband-high-dr";
    why = "wide-band high dynamic-range material: 28-bit VHQ preserves "
          "quiet passages against the resampling noise floor";
  }
  else {
    recipe = SOXR_HQ;
    id = "wideband-general";
    why = "general wide-band material: the 20-bit HQ default matches both "
          "the occupied bandwidth and the dynamic range";
  }

  if (rule && rule_size) {
    strncpy(rule, id, rule_size - 1);
    rule[rule_size - 1] = '\0';
  }
  if (reason && reason_size) {
    snprintf(reason, reason_size,
        "select %s: rms=%.1f dBFS, crest=%.1f dB, occupied_bw=%.3f "
        "Nyquist, hf_frac=%.4f; rule=%s (%s).",
        soxr_auto_recipe_name(recipe), f->rms_dbfs, f->crest_factor_db,
        f->occupied_bw, f->hf_energy_fraction, id, why);
  }
  return recipe;
}



char const * soxr_auto_rationale(soxr_auto_quality_t const * a)
{
  return a && a->reason[0] ? a->reason : "no analysis available";
}



char const * soxr_auto_quality_analyse(
    soxr_datatype_t itype,
    soxr_cbuf_t data, size_t len, unsigned num_channels,
    soxr_auto_quality_t * result)
{
  /* Per channel: 3 edges x 2 one-pole stages, level and band accumulators. */
  #define MAX_STACK_CHANS 8
  #define N_EDGE 3
  #define N_STAGE 2
  float st_stack[MAX_STACK_CHANS][N_EDGE][N_STAGE];
  double e2_stack[MAX_STACK_CHANS];        /* sum of x^2 */
  double pk_stack[MAX_STACK_CHANS];        /* max |x| */
  double band_stack[MAX_STACK_CHANS][4];   /* band energy */
  float *st_c[MAX_STACK_CHANS];
  double *e2_c, *pk_c, *band_c[MAX_STACK_CHANS];
  float alpha[N_EDGE];
  double total_rms = 0, total_peak = 0, band_sum[4] = {0,0,0,0};
  size_t i, analysed;
  unsigned ch;
  int heap = num_channels > MAX_STACK_CHANS;
  float *st_heap = NULL;
  double *band_heap = NULL;

  if (!result)
    return "null result pointer";
  memset(result, 0, sizeof(*result));
  if (num_channels == 0)
    return "num_channels must be > 0";
  if (itype != SOXR_FLOAT32_I && itype != SOXR_FLOAT64_I &&
      itype != SOXR_INT32_I && itype != SOXR_INT16_I)
    return "unsupported input datatype";

  e2_c = e2_stack;
  pk_c = pk_stack;
  if (heap) {
    st_heap = (float *)calloc((size_t)num_channels * N_EDGE * N_STAGE,
                              sizeof(float));
    band_heap = (double *)calloc((size_t)num_channels * 4, sizeof(double));
    if (!st_heap || !band_heap) {
      free(st_heap);
      free(band_heap);
      return "out of memory";
    }
    for (ch = 0; ch < num_channels; ++ch) {
      st_c[ch] = st_heap + (size_t)ch * (N_EDGE * N_STAGE);
      band_c[ch] = band_heap + (size_t)ch * 4;
    }
  }
  else {
    memset(st_stack, 0, sizeof(st_stack));
    memset(e2_stack, 0, sizeof(e2_stack));
    memset(pk_stack, 0, sizeof(pk_stack));
    memset(band_stack, 0, sizeof(band_stack));
    for (ch = 0; ch < num_channels; ++ch) {
      st_c[ch] = (float *)st_stack[ch];
      band_c[ch] = band_stack[ch];
    }
  }

  /* One-pole coefficient for unity DC gain; edge frequency is expressed as a
   * fraction of Nyquist, hence fc(relative to fs) = edge/2. */
  for (i = 0; i < N_EDGE; ++i) {
    double w = 2.0 * M_PI * (band_edge[i] / 2.0);
    alpha[i] = (float)(w / (1.0 + w));
  }

  /* Probe a middle window rather than the stream head: it avoids start-up
   * transients and is more likely to contain periodic content (e.g. a steady
   * beat or syllable rate).  The cost remains bounded by PROBE samples. */
  analysed = len < SOXR_AUTO_PROBE_SAMPLES ? len : SOXR_AUTO_PROBE_SAMPLES;
  {
    size_t probe_start = len > analysed ? (len - analysed) / 2 : 0;
    if (data && analysed) {
      for (i = 0; i < analysed; ++i) {
        int settling = (int)i < TRANSIENT_SKIP;
        for (ch = 0; ch < num_channels; ++ch) {
          float x = read_sample(itype, (void const *)data,
              (probe_start + i) * (size_t)num_channels + ch);
        float ax = fabsf(x);
        float lp[N_EDGE];
        float * st = st_c[ch];
        double * be = band_c[ch];
        unsigned k;

        if (!settling) {
          e2_c[ch] += (double)x * x;
          if ((double)ax > pk_c[ch])
            pk_c[ch] = ax;
        }

        for (k = 0; k < N_EDGE; ++k) {
          float a = alpha[k];
          float * s = st + k * N_STAGE;
          s[0] += a * (x - s[0]);
          s[1] += a * (s[0] - s[1]);
          lp[k] = s[1];                 /* two cascaded one-poles */
          if (settling)
            continue;
        }
        if (settling)
          continue;

        /* Nested low-pass energies partition the spectrum into the 4 bands
         * [0,.125],[.125,.25],[.25,.5],[.5,1] (E(x-y)=Ex-Ey for nested LPFs). */
        be[0] += (double)lp[0] * lp[0];
        be[1] += (double)lp[1] * lp[1] - (double)lp[0] * lp[0];
        be[2] += (double)lp[2] * lp[2] - (double)lp[1] * lp[1];
        be[3] += (double)x * x - (double)lp[2] * lp[2];
        }
      }
    }
  }

  /* Average per-channel statistics (all channels carry equal weight). */
  for (ch = 0; ch < num_channels; ++ch) {
    size_t counted = analysed > TRANSIENT_SKIP
        ? analysed - TRANSIENT_SKIP : analysed;
    double mean_e = counted ? e2_c[ch] / (double)counted : 0;
    size_t b;
    total_rms += sqrt(mean_e > 0 ? mean_e : 0);
    total_peak += pk_c[ch];
    for (b = 0; b < 4; ++b) {
      double v = band_c[ch][b];
      if (v < 0)
        v = 0; /* numerical residue from energy differences */
      band_sum[b] += counted ? v / (double)counted : 0;
    }
  }
  total_rms /= (double)num_channels;
  total_peak /= (double)num_channels;

  {
    double total_band = band_sum[0] + band_sum[1] + band_sum[2] + band_sum[3];
    size_t b;
    double rms_db = total_rms > 0 ? 20 * log10(total_rms) : -999.0;
    double peak_db = total_peak > 0 ? 20 * log10(total_peak) : -999.0;

    result->features.n = analysed;
    result->features.num_channels = num_channels;
    result->features.rms_dbfs = rms_db;
    result->features.peak_dbfs = peak_db;
    result->features.crest_factor_db = peak_db - rms_db;
    if (total_band <= 0)
      total_band = 1e-30;
    for (b = 0; b < 4; ++b)
      result->features.band_energy[b] = band_sum[b] / total_band;
    result->features.hf_energy_fraction = band_sum[3] / total_band;

    /* Highest band individually carrying significant energy. */
    result->features.occupied_bw = band_top[0];
    for (b = 0; b < 4; ++b) {
      if (result->features.band_energy[b] >= BAND_ACTIVE_FRACTION)
        result->features.occupied_bw = band_top[b];
    }
  }

  result->flags = 0;
  result->recipe = soxr_auto_quality_select(&result->features,
      result->rule, sizeof(result->rule),
      result->reason, sizeof(result->reason));

  if (heap) {
    free(st_heap);
    free(band_heap);
  }
  return NULL;
  #undef MAX_STACK_CHANS
  #undef N_EDGE
  #undef N_STAGE
}



soxr_t soxr_auto_create(
    double input_rate, double output_rate, unsigned num_channels,
    soxr_error_t * error,
    soxr_auto_quality_t const * auto_quality,
    soxr_io_spec_t const * io_spec,
    soxr_runtime_spec_t const * runtime_spec,
    soxr_quality_spec_t * chosen_spec_out)
{
  soxr_quality_spec_t q_spec;
  soxr_error_t local_error = NULL;

  if (!auto_quality) {
    local_error = "null auto-quality pointer";
    if (error)
      *error = local_error;
    return NULL;
  }
  q_spec = soxr_quality_spec(auto_quality->recipe, auto_quality->flags);
  if (q_spec.e) {
    local_error = (char const *)q_spec.e;
    if (error)
      *error = local_error;
    return NULL;
  }
  if (chosen_spec_out)
    *chosen_spec_out = q_spec;
  return soxr_create(input_rate, output_rate, num_channels,
      error ? error : &local_error, io_spec, &q_spec, runtime_spec);
}
