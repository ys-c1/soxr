/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details.
 *
 * Channel mapping / down-mix extension, implemented entirely on top of the
 * public resampling API in soxr.h (no internal libsoxr headers are used).
 *
 * The module supplies:
 *   - standard channel layouts (mono, stereo, 5.1, 7.1);
 *   - explicit or label-synthesised down-mix weight matrices;
 *   - per-output-channel linear gain and optional energy normalisation;
 *   - phase-consistency (polarity) validation and full config checking;
 *   - a one-shot "map then resample" convenience wrapper.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "soxr-map.h"

#define MAP_VERSION "libsoxr-map-0.1.3"
#define MAP_MAX SOXR_MAP_MAX_CHANNELS

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif



/* ----------------------------- Ready-made layouts ------------------------- */

soxr_channel_layout_t const soxr_layout_mono = {
  1, {SOXR_CH_FC}
};

soxr_channel_layout_t const soxr_layout_stereo = {
  2, {SOXR_CH_FL, SOXR_CH_FR}
};

soxr_channel_layout_t const soxr_layout_5_1 = {
  6, {SOXR_CH_FL, SOXR_CH_FR, SOXR_CH_FC, SOXR_CH_LFE,
      SOXR_CH_BL, SOXR_CH_BR}
};

soxr_channel_layout_t const soxr_layout_7_1 = {
  8, {SOXR_CH_FL, SOXR_CH_FR, SOXR_CH_FC, SOXR_CH_LFE,
      SOXR_CH_BL, SOXR_CH_BR, SOXR_CH_SL, SOXR_CH_SR}
};



/* ------------------------------- Internal state --------------------------- */

struct soxr_channel_map {
  unsigned in_channels, out_channels;
  unsigned long flags;
  /* Sparse representation of the EFFECTIVE matrix (gain + normalisation
   * applied): for each output row, a list of (input index, coefficient)
   * entries whose coefficient is non-zero. */
  unsigned row_count[MAP_MAX];
  unsigned row_index[MAP_MAX][MAP_MAX];
  double   row_coef [MAP_MAX][MAP_MAX];
  double   gain[MAP_MAX];
};



soxr_channel_map_spec_t soxr_channel_map_spec(void)
{
  soxr_channel_map_spec_t spec;
  memset(&spec, 0, sizeof(spec));
  return spec;
}



char const * soxr_map_version(void)
{
  return MAP_VERSION;
}



static int finite_double(double x)
{
  return x == x && x * 0.0 == 0.0; /* !NaN && +/-Inf excluded */
}



/* Coefficient used when routing input label `in' to output label `out';
 * 0.0 means "no contribution".  Exact label matches are handled separately
 * (unity), so only the non-trivial down-mix rules live here. */
static double auto_coefficient(
    soxr_channel_label_t in, soxr_channel_label_t out)
{
  if (in == out)
    return 1.0;

  /* Centre channel folds into the two front legs at constant power. */
  if (in == SOXR_CH_FC &&
      (out == SOXR_CH_FL || out == SOXR_CH_FR))
    return M_SQRT1_2;

  /* Stereo material folds back to a centre (mono) leg. */
  if (out == SOXR_CH_FC &&
      (in == SOXR_CH_FL || in == SOXR_CH_FR))
    return M_SQRT1_2;

  /* Surround material folds to the centre at half power. */
  if (out == SOXR_CH_FC &&
      (in == SOXR_CH_BL || in == SOXR_CH_BR ||
       in == SOXR_CH_SL || in == SOXR_CH_SR || in == SOXR_CH_BC))
    return 0.5;

  /* Centre channel into the two front legs (constant-power). */
  if (in == SOXR_CH_FLC || in == SOXR_CH_FRC)
    return out == SOXR_CH_FC ? M_SQRT1_2 : 0.0;

  /* Surround / back channels fold down to the front leg on the same side. */
  if (out == SOXR_CH_FL) {
    if (in == SOXR_CH_BL || in == SOXR_CH_SL ||
        in == SOXR_CH_TFL || in == SOXR_CH_TSL ||
        in == SOXR_CH_TBL)
      return M_SQRT1_2;
  }
  if (out == SOXR_CH_FR) {
    if (in == SOXR_CH_BR || in == SOXR_CH_SR ||
        in == SOXR_CH_TFR || in == SOXR_CH_TSR ||
        in == SOXR_CH_TBR)
      return M_SQRT1_2;
  }

  /* Back centre splits equally to both back legs / fronts. */
  if (in == SOXR_CH_BC &&
      (out == SOXR_CH_BL || out == SOXR_CH_BR ||
       out == SOXR_CH_FL || out == SOXR_CH_FR))
    return M_SQRT1_2;

  /* Top / overhead centre material folds to the centre leg. */
  if ((in == SOXR_CH_TFC || in == SOXR_CH_TC || in == SOXR_CH_TBC) &&
      (out == SOXR_CH_FC || out == SOXR_CH_BC))
    return M_SQRT1_2;

  /* LFE is not folded into full-range channels by default. */
  return 0.0;
}



/* Build a layout-driven matrix into `matrix' (row-major), given whether the
 * whole input is a mono layout (special-cased). */
static void build_auto_matrix(soxr_channel_layout_t const * in_layout,
    soxr_channel_layout_t const * out_layout, double * matrix)
{
  unsigned o, i;
  int input_is_mono = in_layout->count == 1 &&
      in_layout->labels[0] == SOXR_CH_FC;

  for (o = 0; o < out_layout->count; ++o) {
    soxr_channel_label_t out_label = out_layout->labels[o];
    for (i = 0; i < in_layout->count; ++i) {
      soxr_channel_label_t in_label = in_layout->labels[i];
      double w;
      if (out_label == in_label)
        w = 1.0;
      else if (input_is_mono &&
               (out_label == SOXR_CH_FL || out_label == SOXR_CH_FR))
        w = 1.0;                 /* mono -> stereo: duplicate at unity */
      else
        w = auto_coefficient(in_label, out_label);
      matrix[(size_t)o * in_layout->count + i] = w;
    }
  }
}



soxr_error_t soxr_channel_map_check(
    soxr_channel_map_spec_t const * spec,
    soxr_channel_map_report_t * report)
{
  double local_matrix[MAP_MAX * MAP_MAX];
  double const * matrix;
  double work[MAP_MAX * MAP_MAX];
  double gain[MAP_MAX];
  unsigned in_ch, out_ch, o, i;
  unsigned unused_inputs[MAP_MAX], silent_outputs[MAP_MAX];
  unsigned n_unused = 0, n_silent = 0, n_phase = 0;
  unsigned long flags;
  double scale = 1.0, energy = 0.0;
  soxr_error_t error = 0;

  if (report)
    memset(report, 0, sizeof(*report));
  if (!spec)
    return "null mapping spec";

  in_ch = spec->in_channels;
  out_ch = spec->out_channels;
  flags = spec->flags;

  if (in_ch == 0 || in_ch > MAP_MAX)
    return "invalid # of input channels";
  if (out_ch == 0 || out_ch > MAP_MAX)
    return "invalid # of output channels";
  if (spec->in_layout && spec->in_layout->count != in_ch)
    return "input layout does not match # of input channels";
  if (spec->out_layout && spec->out_layout->count != out_ch)
    return "output layout does not match # of output channels";

  for (i = 0; i < in_ch; ++i)
    if (spec->in_layout && spec->in_layout->labels[i] >= SOXR_CH_MAXLABEL)
      return "invalid input channel label";
  for (o = 0; o < out_ch; ++o)
    if (spec->out_layout && spec->out_layout->labels[o] >= SOXR_CH_MAXLABEL)
      return "invalid output channel label";

  for (o = 0; o < out_ch; ++o) {
    gain[o] = spec->gain ? spec->gain[o] : 1.0;
    if (!finite_double(gain[o]))
      return "non-finite channel gain";
  }

  if (spec->matrix)
    matrix = spec->matrix;
  else {
    if (!spec->in_layout || !spec->out_layout)
      return "no matrix and no input/output layouts given";
    build_auto_matrix(spec->in_layout, spec->out_layout, local_matrix);
    matrix = local_matrix;
  }

  for (o = 0; o < out_ch; ++o)
    for (i = 0; i < in_ch; ++i) {
      double w = matrix[(size_t)o * in_ch + i];
      if (!finite_double(w))
        return "non-finite matrix coefficient";
      if (w < 0.0)
        ++n_phase;
    }

  /* Effective matrix = per-output gain, then optional energy scaling. */
  for (o = 0; o < out_ch; ++o)
    for (i = 0; i < in_ch; ++i) {
      double w = matrix[(size_t)o * in_ch + i] * gain[o];
      work[(size_t)o * in_ch + i] = w;
      energy += w * w;
    }

  if ((flags & SOXR_MX_NORMALISE_ENERGY)) {
    if (energy <= 0.0)
      return "cannot normalise an all-zero matrix";
    /* For independent inputs with equal per-channel energy E, total input
     * energy is in_ch * E and output energy is (sum of squared weights) * E.
     * Scale so that the sum of squared weights equals in_ch, i.e. the total
     * signal energy is preserved. */
    scale = sqrt((double)in_ch / energy);
  }

  for (o = 0; o < out_ch; ++o)
    for (i = 0; i < in_ch; ++i)
      work[(size_t)o * in_ch + i] *= scale;

  /* Structural checks on the effective matrix. */
  for (i = 0; i < in_ch; ++i) {
    int used = 0;
    for (o = 0; o < out_ch; ++o)
      if (work[(size_t)o * in_ch + i] != 0.0) { used = 1; break; }
    if (!used)
      unused_inputs[n_unused++] = i;
  }
  for (o = 0; o < out_ch; ++o) {
    int sounding = 0;
    for (i = 0; i < in_ch; ++i)
      if (work[(size_t)o * in_ch + i] != 0.0) { sounding = 1; break; }
    if (!sounding)
      silent_outputs[n_silent++] = o;
  }

  if (n_phase && !(flags & SOXR_MX_ALLOW_POLARITY_INVERT))
    error = "phase-inverting coefficient not permitted";
  else if ((flags & SOXR_MX_STRICT) && (n_unused || n_silent))
    error = "strict mapping has unused input or silent output channel";

  if (report) {
    report->num_in_channels = in_ch;
    report->num_out_channels = out_ch;
    report->num_unused_inputs = n_unused;
    report->num_silent_outputs = n_silent;
    report->num_phase_conflicts = n_phase;
    for (i = 0; i < n_unused; ++i)
      report->unused_inputs[i] = unused_inputs[i];
    for (o = 0; o < n_silent; ++o)
      report->silent_outputs[o] = silent_outputs[o];
    for (o = 0; o < out_ch; ++o)
      for (i = 0; i < in_ch; ++i)
        report->effective_matrix[(size_t)o * in_ch + i] =
            work[(size_t)o * in_ch + i];
  }
  return error;
}



soxr_channel_map_t soxr_channel_map_create(
    soxr_channel_map_spec_t const * spec,
    soxr_error_t * error0)
{
  soxr_channel_map_t map = 0;
  soxr_channel_map_report_t report;
  soxr_error_t error;
  unsigned in_ch, out_ch, o, i;

  error = soxr_channel_map_check(spec, &report);
  if (!error) {
    map = calloc(1, sizeof(*map));
    if (!map)
      error = "malloc failed";
  }
  if (map) {
    in_ch = spec->in_channels;
    out_ch = spec->out_channels;
    map->in_channels = in_ch;
    map->out_channels = out_ch;
    map->flags = spec->flags;
    for (o = 0; o < out_ch; ++o) {
      map->gain[o] = spec->gain ? spec->gain[o] : 1.0;
      for (i = 0; i < in_ch; ++i) {
        double w = report.effective_matrix[(size_t)o * in_ch + i];
        if (w != 0.0) {
          map->row_index[o][map->row_count[o]] = i;
          map->row_coef [o][map->row_count[o]] = w;
          ++map->row_count[o];
        }
      }
    }
  }
  if (error0)
    *error0 = error;
  return map;
}



void soxr_channel_map_delete(soxr_channel_map_t map)
{
  free(map);
}



soxr_error_t soxr_channel_map_apply(
    soxr_channel_map_t map,
    float const * in, size_t len, int interleaved,
    float * out)
{
  unsigned o, k;
  size_t n;

  if (!map)
    return "null channel map";
  if (!in && len)
    return "null input buffer";
  if (!out && len)
    return "null output buffer";

  if (interleaved) {
    unsigned in_ch = map->in_channels, out_ch = map->out_channels;
    for (n = 0; n < len; ++n) {
      for (o = 0; o < out_ch; ++o) {
        double acc = 0.0;
        for (k = 0; k < map->row_count[o]; ++k)
          acc += map->row_coef[o][k] *
                 (double)in[n * in_ch + map->row_index[o][k]];
        out[n * out_ch + o] = (float)acc;
      }
    }
  }
  else {
    /* Split mode: `in' / `out' are arrays of per-channel buffer pointers. */
    float const * const * in_bufs = (float const * const *)in;
    float * * out_bufs = (float * *)out;
    for (o = 0; o < map->out_channels; ++o) {
      float * dst = out_bufs[o];
      for (n = 0; n < len; ++n) {
        double acc = 0.0;
        for (k = 0; k < map->row_count[o]; ++k)
          acc += map->row_coef[o][k] *
                 (double)in_bufs[map->row_index[o][k]][n];
        dst[n] = (float)acc;
      }
    }
  }
  return 0;
}



soxr_error_t soxr_map_resample_oneshot(
    soxr_channel_map_spec_t const * map_spec,
    double input_rate,
    double output_rate,
    float const * in, size_t ilen,
    float * out, size_t olen, size_t * odone,
    soxr_quality_spec_t const * q_spec,
    soxr_runtime_spec_t const * runtime_spec)
{
  soxr_channel_map_t map = 0;
  soxr_t resampler = 0;
  soxr_io_spec_t io_spec;
  soxr_error_t error;
  float * mapped = 0;
  float * resampled = 0;
  float * in_ptrs[MAP_MAX];
  float * out_ptrs[MAP_MAX];
  unsigned i, o;
  size_t n;
  size_t map_odone = 0;

  if (!in && ilen)
    return "null input buffer";
  if (!out && olen)
    return "null output buffer";

  map = soxr_channel_map_create(map_spec, &error);
  if (error)
    goto done;

  mapped = malloc(sizeof(float) * ilen * map->out_channels);
  resampled = malloc(sizeof(float) * olen * map->out_channels);
  if (!mapped || !resampled) { error = "malloc failed"; goto done; }

  error = soxr_channel_map_apply(map, in, ilen, 1, mapped);
  if (error)
    goto done;

  /* Split (per-channel contiguous) buffers for the resampler. */
  for (i = 0; i < map->out_channels; ++i) {
    in_ptrs[i] = mapped + (size_t)i * ilen;
    out_ptrs[i] = resampled + (size_t)i * olen;
  }

  io_spec = soxr_io_spec(SOXR_FLOAT32_S, SOXR_FLOAT32_S);
  resampler = soxr_create(input_rate, output_rate, map->out_channels,
      &error, &io_spec, q_spec, runtime_spec);
  if (error)
    goto done;

  {
    size_t done = 0;
    error = soxr_process(resampler, in_ptrs, ilen, 0,
        out_ptrs, olen, &done);
    map_odone = done;

    /* Drain the resampler's filter delay: keep calling with NULL input,
     * advancing each split output pointer by what has been produced. */
    while (!error && done && map_odone < olen) {
      size_t remain = olen - map_odone;
      for (i = 0; i < map->out_channels; ++i)
        out_ptrs[i] = resampled + (size_t)i * olen + map_odone;
      done = 0;
      error = soxr_process(resampler, 0, 0, 0, out_ptrs, remain, &done);
      map_odone += done;
    }
  }

  /* Interleave the resampled channels into the caller's output buffer. */
  if (!error) {
    for (n = 0; n < map_odone; ++n)
      for (o = 0; o < map->out_channels; ++o)
        out[n * map->out_channels + o] = resampled[(size_t)o * olen + n];
  }

done:
  if (resampler)
    soxr_delete(resampler);
  free(resampled);
  free(mapped);
  soxr_channel_map_delete(map);
  if (odone)
    *odone = error ? 0 : map_odone;
  return error;
}
