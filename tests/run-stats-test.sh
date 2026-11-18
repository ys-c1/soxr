#!/bin/sh
# Build soxr (if needed) and run the soxr-stats acceptance test, in one go.
# Usage: tests/run-stats-test.sh [build-dir]
set -e

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BUILD=${1:-"$ROOT/artifacts/build"}

if [ ! -f "$BUILD/src/libsoxr.so" ]; then
  mkdir -p "$BUILD"
  cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF -DWITH_OPENMP=OFF \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5
  cmake --build "$BUILD" -j"$(nproc)"
fi

cc -O2 -Wall -Wextra -I"$ROOT/src" -I"$ROOT/tests" \
  "$ROOT/tests/soxr-stats.c" "$ROOT/tests/soxr-stats-test.c" \
  -o "$BUILD/soxr-stats-test" \
  -L"$BUILD/src" -lsoxr -lm -lpthread -Wl,-rpath,"$BUILD/src"

exec "$BUILD/soxr-stats-test"
