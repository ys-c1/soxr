/* SoX Resampler Library      Copyright (c) 2007-13 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Streaming statistics layer on top of the public soxr API, for real-time
 * streaming scenarios that need to know processing latency and throughput:
 *
 *   - latency report   : input-to-output delay, in output samples, input
 *                        samples, and milliseconds (per soxr_delay()).
 *   - throughput stats : processed samples per second (input and output).
 *   - block histogram  : wall-clock time spent inside soxr_process() per
 *                        call, accumulated into log2-scaled microsecond bins.
 *
 * All query/reset functions and soxr_stats_process() are thread-safe (a
 * per-object mutex serialises statistic updates and snapshots).  Note that
 * thread-safety of the *statistics* does not lift the underlying soxr
 * restriction: concurrent soxr_stats_process() calls on objects attached to
 * the *same* soxr_t instance are not permitted; attach one soxr_stats_t per
 * soxr_t stream.  Concurrent queries from any thread are always safe.
 *
 * This file is part of the tests/ additions; it is not installed. */

#ifndef SOXR_STATS_H
#define SOXR_STATS_H

#include <soxr.h>

#if defined __cplusplus
extern "C" {
#endif

#define SOXR_STATS_HIST_BINS 24   /* Bin i (i>0): [2^(i-1), 2^i) us; 0: <1us. */

typedef struct soxr_stats soxr_stats_t;

/* Attach a statistics collector to an existing stream resampler.  The
 * resampler is NOT owned by the returned object (caller still soxr_delete()s
 * it).  irate/orate must match those given to soxr_create(); they are used
 * for the latency/throughput unit conversions.  Returns 0 on malloc failure. */
soxr_stats_t * soxr_stats_create(soxr_t resampler, double irate, double orate);
void           soxr_stats_delete(soxr_stats_t *);

/* Timing wrapper around soxr_process(): identical arguments and semantics,
 * plus per-call wall-clock accounting.  idone/odone may be 0 as per
 * soxr_process(). */
soxr_error_t soxr_stats_process(soxr_stats_t *,
    void const * in, size_t ilen, size_t * idone,
    void       * out, size_t olen, size_t * odone);

/* Reset all accumulated counters/timers (latency is unaffected: it is a
 * property of the resampler, not of the history). */
void soxr_stats_reset(soxr_stats_t *);

typedef struct {
  size_t blocks;                /* Number of timed soxr_process() calls.     */
  double in_samples;            /* Total input  samples consumed (per ch.).  */
  double out_samples;           /* Total output samples produced (per ch.).  */
  double process_seconds;       /* Cumulative time inside soxr_process().    */
  double latency_out_samples;   /* Steady-state input->output delay, in      */
                                /* output samples (soxr_delay() sampled while */
                                /* input is being consumed; falls back to the */
                                /* current soxr_delay() before any input).    */
  double latency_in_samples;    /* Same, expressed in input samples.         */
  double latency_ms;            /* Same, in milliseconds at the output rate. */
  double throughput_in_sps;     /* in_samples  / process_seconds.            */
  double throughput_out_sps;    /* out_samples / process_seconds.            */
  unsigned long hist[SOXR_STATS_HIST_BINS];  /* Per-block time histogram.    */
} soxr_stats_report_t;

/* Take a consistent snapshot of all statistics (single lock acquisition). */
void soxr_stats_report(soxr_stats_t *, soxr_stats_report_t * report);

/* Microsecond range covered by histogram bin `bin' (for display). */
void soxr_stats_hist_bin_us(unsigned bin, double * lo_us, double * hi_us);

#if defined __cplusplus
} /* extern "C" */
#endif
#endif
