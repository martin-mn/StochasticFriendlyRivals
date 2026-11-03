#!/bin/bash
# Build the W-separation driver.  Same flags as ../cert/compile.sh: no fast-math, no FP contraction
# (the interval engine's outward rounding assumes correctly rounded single operations).
set -e
cd "$(dirname "$0")"
${CC:-cc} -O2 -ffp-contract=off -fno-fast-math -o wsep wsep.c -lm -lpthread
echo "built $(pwd)/wsep"
