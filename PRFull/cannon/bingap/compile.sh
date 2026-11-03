#!/bin/bash
# Cannon build (gcc).  No -ffast-math / FP contraction: the interval engine needs correctly rounded ops.
set -e
cd "$(dirname "$0")"
${CC:-gcc} -O2 -ffp-contract=off -fno-fast-math -o bingap bingap.c certeng.c -lm -lpthread
echo "built $(pwd)/bingap"
