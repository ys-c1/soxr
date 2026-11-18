# Channel mapping / down-mix extension and acceptance tests

This directory holds the acceptance test for the channel-mapping extension
that lives alongside the public resampling API:

- `../src/soxr-map.h` — public API of the extension (layouts, mapping spec,
  validation report, apply / one-shot functions).
- `../src/soxr-map.c` — implementation, built only on top of the public
  `soxr.h` API (no internal libsoxr headers, no new third-party code).
- `channel-map-test.c` — acceptance program (all test signals are generated
  in code; nothing is downloaded).
- `run-tests` — configures libsoxr into `../artifacts/soxr-build`, builds the
  extension and this test, and runs everything in one shot.

## Run

    ./run-tests

Environment:

- `SOXR_MAP_TIME_LIMIT_MS` overrides the throughput gate (default 50 ms);
  useful when running under a sanitizer, e.g.
  `SOXR_MAP_TIME_LIMIT_MS=500 ./run-tests`.

## What is verified

1. Channel mapping for mono, stereo, 5.1 and 7.1 against hand-computed
   values (explicit, configurable weight matrices).
2. Layout-driven (label-based) matrix synthesis and the validation report
   (unused inputs, silent outputs, phase conflicts, effective matrix).
3. Down-mix total-energy deviation; the energy-normalised configuration must
   stay within +/-0.5 dB of the total input energy.
4. Throughput: 10 s @ 48 kHz, 5.1 -> stereo, well under 50 ms.
5. Rejection of invalid configurations (bad channel counts, missing matrix,
   layout mismatch, non-finite weights/gains, polarity inversion, strict
   unused/silent channels, normalisation of a zero matrix).
6. Interleaved vs split-buffer equivalence and an end-to-end
   map + resample one-shot call.
