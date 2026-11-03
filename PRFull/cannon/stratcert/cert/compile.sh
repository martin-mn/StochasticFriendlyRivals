#!/bin/bash
# Build the certification engine.  macOS: universal dylib (the Anaconda python is x86_64 and runs
# under Rosetta; the arm64 slice serves a native python).  Linux (Cannon): libcert.so.
# No -ffast-math and no FP contraction: the interval engine's outward rounding assumes correctly
# rounded single operations (see the proof at dn()/up() in certeng.c).
set -e
cd "$(dirname "$0")"
if [ "$(uname)" = "Darwin" ]; then
  cc -O2 -ffp-contract=off -fno-fast-math -arch arm64 -arch x86_64 -fPIC -shared \
     -o libcert.dylib certeng.c -lm -lpthread
  echo "built $(pwd)/libcert.dylib"
else
  ${CC:-gcc} -O2 -ffp-contract=off -fno-fast-math -fPIC -shared -o libcert.so certeng.c -lm -lpthread
  echo "built $(pwd)/libcert.so"
fi
