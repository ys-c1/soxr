/* Smoke tests for the auto-quality selector public API:
 *  - all four supported input datatypes analyse without error and select a
 *    valid standard recipe;
 *  - edge cases (NULL data, zero length, zero channels, bad datatype, NULL
 *    result) are handled with documented error strings, no crashes;
 *  - the rationale is printable;
 *  - soxr_auto_create produces a working resampler whose output matches the
 *    equivalent manual soxr_create with the selected recipe.
 */

#include <soxr.h>
#include "../src/soxr-auto.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(cond) do { \
  if (!(cond)) { fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #cond); \
                 return 1; } } while (0)

static size_t mono_resample(soxr_t r, double const * in, size_t ilen,
                            double * out, size_t ocap)
{
  size_t idone = 0, odone = 0, od2 = 0;
  soxr_process(r, (soxr_in_t)in, ilen, &idone, out, ocap, &odone);
  soxr_process(r, NULL, 0, NULL, out + odone, ocap - odone, &od2);
  return odone + od2;
}

int main(void)
{
  enum { N = 4096 };
  size_t i;
  double d[N];
  float f[N];
  int32_t i32[N];
  int16_t i16[N];
  soxr_auto_quality_t a;
  char const * err;

  /* Build a non-trivial broad-band signal. */
  for (i = 0; i < N; ++i) {
    double t = (double)i / 48000.0;
    d[i] = 0.4 * sin(2 * M_PI * 440 * t) + 0.2 * sin(2 * M_PI * 9000 * t);
    f[i] = (float)d[i];
    i32[i] = (int32_t)(d[i] * 2147483647.0);
    i16[i] = (int16_t)(d[i] * 32767.0);
  }
  /* All four datatypes analyse successfully and pick a standard recipe. */
  err = soxr_auto_quality_analyse(SOXR_FLOAT32_I, f, N, 1, &a);
  CHECK(err == NULL && a.recipe <= SOXR_VHQ && a.reason[0]);
  err = soxr_auto_quality_analyse(SOXR_FLOAT64_I, d, N, 1, &a);
  CHECK(err == NULL && a.recipe <= SOXR_VHQ);
  err = soxr_auto_quality_analyse(SOXR_INT32_I, i32, N, 1, &a);
  CHECK(err == NULL && a.recipe <= SOXR_VHQ);
  err = soxr_auto_quality_analyse(SOXR_INT16_I, i16, N, 1, &a);
  CHECK(err == NULL && a.recipe <= SOXR_VHQ);
  printf("datatypes OK; f64 rationale: %s\n", soxr_auto_rationale(&a));

  /* Multichannel (interleaved stereo). */
  {
    float st[N * 2];
    for (i = 0; i < N; ++i) { st[2*i] = f[i]; st[2*i+1] = f[i]; }
    err = soxr_auto_quality_analyse(SOXR_FLOAT32_I, st, N, 2, &a);
    CHECK(err == NULL && a.features.num_channels == 2);
  }

  /* Edge cases. */
  CHECK(soxr_auto_quality_analyse(SOXR_FLOAT32_I, f, N, 1, NULL) != NULL);
  CHECK(soxr_auto_quality_analyse(SOXR_FLOAT32_I, f, N, 0, &a) != NULL);
  CHECK(soxr_auto_quality_analyse((soxr_datatype_t)99, f, N, 1, &a) != NULL);
  err = soxr_auto_quality_analyse(SOXR_FLOAT32_I, NULL, N, 1, &a);
  CHECK(err == NULL);                                   /* NULL data allowed */
  err = soxr_auto_quality_analyse(SOXR_FLOAT32_I, f, 0, 1, &a);
  CHECK(err == NULL && a.features.n == 0);
  CHECK(soxr_auto_create(48000, 44100, 1, NULL, NULL, NULL, NULL, NULL)
        == NULL);                                       /* NULL selection */

  /* Silence selects QQ. */
  memset(d, 0, sizeof(d));
  err = soxr_auto_quality_analyse(SOXR_FLOAT64_I, d, N, 1, &a);
  CHECK(err == NULL && a.recipe == SOXR_QQ);
  CHECK(strcmp(a.rule, "silence") == 0);

  /* End-to-end: auto-create output equals a manual create with the recipe. */
  for (i = 0; i < N; ++i)
    d[i] = 0.4 * sin(2 * M_PI * 440 * (double)i / 48000.0)
         + 0.2 * sin(2 * M_PI * 9000 * (double)i / 48000.0);
  {
    soxr_auto_quality_t sel;
    soxr_quality_spec_t qs_auto, qs_man;
    soxr_io_spec_t io = soxr_io_spec(SOXR_FLOAT64_I, SOXR_FLOAT64_I);
    soxr_error_t e1 = NULL, e2 = NULL;
    soxr_t r1, r2;
    double o1[N], o2[N];
    size_t n1, n2, k, mism = 0;

    CHECK(soxr_auto_quality_analyse(SOXR_FLOAT64_I, d, N, 1, &sel) == NULL);
    r1 = soxr_auto_create(48000, 16000, 1, &e1, &sel, &io, NULL, &qs_auto);
    qs_man = soxr_quality_spec(sel.recipe, sel.flags);
    r2 = soxr_create(48000, 16000, 1, &e2, &io, &qs_man, NULL);
    CHECK(r1 != NULL && r2 != NULL && e1 == NULL && e2 == NULL);
    CHECK(qs_auto.precision == qs_man.precision &&
          qs_auto.passband_end == qs_man.passband_end);
    n1 = mono_resample(r1, d, N, o1, N);
    n2 = mono_resample(r2, d, N, o2, N);
    CHECK(n1 == n2);
    for (k = 0; k < n1; ++k)
      if (o1[k] != o2[k])
        ++mism;
    CHECK(mism == 0);
    soxr_delete(r1);
    soxr_delete(r2);
  }

  printf("autoq-smoke: all checks passed\n");
  return 0;
}
