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



/* Auto-quality selector: a thin layer on top of the public soxr API.
 *
 * It analyses input-signal characteristics (occupied bandwidth and dynamic
 * range / crest factor), selects one of the standard soxr quality recipes
 * (SOXR_QQ .. SOXR_VHQ), and records a printable rationale.  It does not
 * modify the resampling kernel and introduces no dependency beyond libc/libm:
 * analysis is a bounded, single-pass O(PROBE) streaming sweep using cascaded
 * first-order IIR band-split filters and running level statistics, so its
 * cost is independent of the total stream length.
 *
 * Typical use:
 *
 *   soxr_auto_quality_t aq;
 *   soxr_auto_quality_analyse(SOXR_FLOAT32_I, data, n, 1, &aq);
 *   soxr_t r = soxr_auto_create(48000, 44100, 1, NULL, &aq, NULL, NULL, NULL);
 *   fputs(soxr_auto_rationale(&aq), stderr);
 *   ... soxr_process(r, ...) as usual ...
 */

#if !defined soxr_auto_included
#define soxr_auto_included

#include "soxr.h"

/* soxr.h #undef's its SOXR export macro at the end; re-derive an equivalent
 * one for the symbols exported by this selector layer. */
#if defined SOXR_DLL
  #if defined soxr_EXPORTS
    #define SOXR_AUTO_API __declspec(dllexport)
  #else
    #define SOXR_AUTO_API __declspec(dllimport)
  #endif
#elif defined SOXR_VISIBILITY && defined __GNUC__ && \
      (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1))
  #define SOXR_AUTO_API __attribute__ ((visibility("default")))
#else
  #define SOXR_AUTO_API
#endif

#if defined __cplusplus
extern "C" {
#endif

#define SOXR_AUTO_THIS_VERSION 1u

/* Names of the standard recipes, indexed by SOXR_QQ .. SOXR_VHQ (0..4). */
SOXR_AUTO_API char const * soxr_auto_recipe_name(unsigned long recipe);

/* Measured characteristics of the analysed input.
 * All level quantities are in dBFS; bandwidth quantities are fractions of
 * the input Nyquist frequency (0..1). */
typedef struct soxr_auto_features {
  size_t n;                 /* Number of samples per channel analysed.     */
  unsigned num_channels;    /* Number of channels analysed.                */
  double rms_dbfs;          /* Mean RMS level, dBFS.                       */
  double peak_dbfs;         /* Mean absolute-peak level, dBFS.             */
  double crest_factor_db;   /* peak_dbfs - rms_dbfs: dynamic-range proxy.  */
  double occupied_bw;       /* Highest band carrying significant energy.   */
  double hf_energy_fraction;/* Energy above 0.5 Nyquist / total energy.    */
  double band_energy[4];    /* Energy per band, fractions summing to 1.    */
} soxr_auto_features_t;

/* Result of a selection. */
typedef struct soxr_auto_quality {
  soxr_auto_features_t features; /* Input features the decision is based on. */
  unsigned long recipe;          /* Selected soxr quality recipe.           */
  unsigned long flags;           /* Flags to pass to soxr_quality_spec().   */
  char reason[320];              /* Human-readable, one-line rationale.     */
  char rule[32];                 /* Machine-readable id of the rule fired.  */
} soxr_auto_quality_t;

/* Analyse interleaved input of the given soxr datatype (FLOAT32_I,
 * FLOAT64_I, INT32_I or INT16_I).  `len' is samples per channel.
 * Returns a soxr error string (NULL == no error).  No heap allocation for
 * 8 channels or fewer; more channels use two temporary calloc buffers. */
SOXR_AUTO_API char const * soxr_auto_quality_analyse(
    soxr_datatype_t itype,
    soxr_cbuf_t data, size_t len, unsigned num_channels,
    soxr_auto_quality_t * result);

/* Select a recipe from already-computed features (separates the analysis
 * cost from the decision; soxr_auto_quality_analyse calls this internally). */
SOXR_AUTO_API unsigned long soxr_auto_quality_select(
    soxr_auto_features_t const * features,
    char * rule, size_t rule_size,
    char * reason, size_t reason_size);

/* Convenience wrapper: build the quality spec from the selection and create
 * the resampler via the public soxr_create().  Any of io/runtime/error may
 * be NULL; the selected spec is passed through soxr_quality_spec() only. */
SOXR_AUTO_API soxr_t soxr_auto_create(
    double input_rate, double output_rate, unsigned num_channels,
    soxr_error_t * error,
    soxr_auto_quality_t const * auto_quality,
    soxr_io_spec_t const * io_spec,
    soxr_runtime_spec_t const * runtime_spec,
    soxr_quality_spec_t * chosen_spec_out); /* may be NULL */

/* Return the rationale line stored in the selection (never NULL). */
SOXR_AUTO_API char const * soxr_auto_rationale(soxr_auto_quality_t const *);

#if defined __cplusplus
}
#endif

#endif
