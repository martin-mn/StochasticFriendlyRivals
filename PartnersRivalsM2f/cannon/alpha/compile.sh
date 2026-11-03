#!/bin/bash
set -e
cc -O3 -march=native -o census census.c m2f.c -lm
echo "built census (psweep)"
