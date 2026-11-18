/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details.
 *
 * Channel mapping / down-mix extension, built on top of the libsoxr
 * resampling API declared in soxr.h.
 *
 * Conventions (same as soxr.h):
 *   Buffer lengths are expressed as the number of samples per channel.
 *   Weight matrices are stored row-major, one row per OUTPUT channel and
 *   one column per INPUT channel:  out[o] = sum_i matrix[o][i] * in[i].
 */

#if !defined soxr_map_included
#define soxr_map_included

#include "soxr.h"

/* soxr.h #undef's its export macro at the end; re-establish an equivalent
 * one for the symbols declared in this header. */
#if defined SOXR_DLL
  #if defined soxr_EXPORTS
    #define SOXR __declspec(dllexport)
  #else
    #define SOXR __declspec(dllimport)
  #endif
#elif defined SOXR_VISIBILITY && defined __GNUC__ && \
      (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1))
  #define SOXR __attribute__ ((visibility("default")))
#else
  #define SOXR
#endif

#if defined __cplusplus
extern "C" {
#endif

/* Maximum number of channels addressable by a mapping configuration. */
#define SOXR_MAP_MAX_CHANNELS 32



/* Speaker / channel labels, in the fixed ITU-R BS.2051 (SMPTE/Windows)
 * channel order.  A "layout" is an ordered sequence of these labels, i.e.
 * it defines both the channel count and the meaning of each slot. */

typedef enum {
  SOXR_CH_UNKNOWN = 0,
  SOXR_CH_FL,                    /* Front left.                          */
  SOXR_CH_FR,                    /* Front right.                         */
  SOXR_CH_FC,                    /* Front centre.                        */
  SOXR_CH_LFE,                   /* Low-frequency effects (sub-woofer).  */
  SOXR_CH_BL,                    /* Back / surround left.                */
  SOXR_CH_BR,                    /* Back / surround right.               */
  SOXR_CH_FLC,                   /* Front left-of-centre.                */
  SOXR_CH_FRC,                   /* Front right-of-centre.               */
  SOXR_CH_BC,                    /* Back centre.                         */
  SOXR_CH_SL,                    /* Side surround left.                  */
  SOXR_CH_SR,                    /* Side surround right.                 */
  SOXR_CH_TFL,                   /* Top front left.                      */
  SOXR_CH_TFR,                   /* Top front right.                     */
  SOXR_CH_TFC,                   /* Top front centre.                    */
  SOXR_CH_TC,                    /* Top centre.                          */
  SOXR_CH_TBL,                   /* Top back left.                       */
  SOXR_CH_TBR,                   /* Top back right.                      */
  SOXR_CH_TBC,                   /* Top back centre.                     */
  SOXR_CH_TSL,                   /* Top side left.                       */
  SOXR_CH_TSR,                   /* Top side right.                      */
  SOXR_CH_LLFE,                  /* Second LFE.                          */
  SOXR_CH_SLFE,                  /* Third LFE.                           */
  SOXR_CH_MAXLABEL
} soxr_channel_label_t;

/* A channel layout: `count' labels in the order the channels are stored. */
typedef struct {
  unsigned count;
  soxr_channel_label_t labels[SOXR_MAP_MAX_CHANNELS];
} soxr_channel_layout_t;

/* Ready-made layouts in standard channel order. */
SOXR extern soxr_channel_layout_t const soxr_layout_mono;    /* 1.0: FC        */
SOXR extern soxr_channel_layout_t const soxr_layout_stereo;  /* 2.0: FL, FR    */
SOXR extern soxr_channel_layout_t const soxr_layout_5_1;     /* 5.1: FL FR FC LFE BL BR */
SOXR extern soxr_channel_layout_t const soxr_layout_7_1;     /* 7.1: FL FR FC LFE BL BR SL SR */



/* Flags for soxr_channel_map_spec_t.flags. */

/* Scale the whole weight matrix so that, for independent equal-energy
 * input channels, the summed output energy equals the total input energy:
 * the sum of all squared effective coefficients becomes in_channels. */
#define SOXR_MX_NORMALISE_ENERGY 1u

/* Reject configurations in which any input channel is unused or any
 * output channel is silent (all-zero row). */
#define SOXR_MX_STRICT           2u

/* By default a negative (phase-inverting) coefficient is a configuration
 * error.  Set this flag to permit polarity inversion; the inversion is
 * still reported by soxr_channel_map_check(). */
#define SOXR_MX_ALLOW_POLARITY_INVERT 4u



/* Configuration of a channel mapping / down-mix.
 *
 * Two ways to define the routing, at least one of which is required:
 *
 *  1. Explicit weight matrix: set `matrix' to a flat array of
 *     out_channels * in_channels doubles, row-major (output-major).
 *     Per-output-channel gains in `gain' are applied on top of it.
 *
 *  2. Layout-driven matrix: leave `matrix' NULL and set `in_layout' /
 *     `out_layout'.  A conventional down-/up-mix matrix is synthesised
 *     from the channel labels (see soxr-map.c for the exact rules).
 *
 * The two may be combined: if `matrix' is given it is always used and the
 * layouts only serve as documentation / label-based validation inputs.
 */
typedef struct {
  unsigned in_channels;
  unsigned out_channels;

  double const * matrix;              /* NULL => derive from the layouts.   */
  soxr_channel_layout_t const * in_layout;   /* Optional. */
  soxr_channel_layout_t const * out_layout;  /* Optional. */

  /* Per-OUTPUT-channel linear gain, applied after the matrix mix.
   * NULL (or 0.0 entries) means unity gain. */
  double const * gain;

  unsigned long flags;
  void * e;                           /* Reserved for internal use. */
} soxr_channel_map_spec_t;

/* Result of validating / inspecting a mapping. */
typedef struct {
  unsigned num_in_channels;
  unsigned num_out_channels;
  unsigned num_unused_inputs;         /* Inputs with an all-zero column.    */
  unsigned num_silent_outputs;        /* Outputs with an all-zero row.      */
  unsigned num_phase_conflicts;       /* Negative effective coefficients.   */
  unsigned unused_inputs[SOXR_MAP_MAX_CHANNELS];
  unsigned silent_outputs[SOXR_MAP_MAX_CHANNELS];
  /* Effective coefficients actually applied (gain and normalisation
   * included); row-major, out_channels rows of in_channels. */
  double effective_matrix[SOXR_MAP_MAX_CHANNELS * SOXR_MAP_MAX_CHANNELS];
} soxr_channel_map_report_t;

typedef struct soxr_channel_map * soxr_channel_map_t;



/* Default-initialise a spec (call first, then override fields). */
SOXR soxr_channel_map_spec_t soxr_channel_map_spec(void);

/* Build a mapping from a spec.  Returns NULL and, if `error' is non-NULL,
 * sets *error to a static string describing the configuration error. */
SOXR soxr_channel_map_t soxr_channel_map_create(
    soxr_channel_map_spec_t const * spec,
    soxr_error_t * error);

SOXR void soxr_channel_map_delete(soxr_channel_map_t);

/* Validate a spec without building a mapping; always fills `report' (when
 * non-NULL).  Returns NULL on success, or a static error string when the
 * configuration is rejected (out-of-range sizes, non-finite weights,
 * polarity inversion without SOXR_MX_ALLOW_POLARITY_INVERT, strict-mode
 * unused/silent channels, ...). */
SOXR soxr_error_t soxr_channel_map_check(
    soxr_channel_map_spec_t const * spec,
    soxr_channel_map_report_t * report);

/* Apply the mapping to a block of 32-bit float samples.
 * `interleaved' != 0: one channel-interleaved buffer, length in samples
 *                     per channel; in  has in_channels  * len floats,
 *                     out has out_channels * len floats (in-place is OK
 *                     only when out_channels <= in_channels).
 * `interleaved' == 0: arrays of per-channel buffers; `in' / `out' point at
 *                     in_channels / out_channels float buffers of `len'. */
SOXR soxr_error_t soxr_channel_map_apply(
    soxr_channel_map_t map,
    float const * in, size_t len, int interleaved,
    float * out);

/* One-shot: channel-map then resample a signal held entirely in memory.
 * Input and output are channel-interleaved 32-bit float.  The mapping and
 * the resampler are created internally and destroyed before returning. */
SOXR soxr_error_t soxr_map_resample_oneshot(
    soxr_channel_map_spec_t const * map_spec,
    double input_rate,
    double output_rate,
    float const * in, size_t ilen,
    float * out, size_t olen, size_t * odone,
    soxr_quality_spec_t const * q_spec,
    soxr_runtime_spec_t const * runtime_spec);

SOXR char const * soxr_map_version(void);

#undef SOXR

#if defined __cplusplus
}
#endif

#endif
