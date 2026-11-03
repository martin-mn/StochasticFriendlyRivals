#!/bin/bash
# Build the leading-order GTH engine.  Plain C99, libm only, no modules needed.
# No -march=native: the login node and the sapphire compute nodes are not the same chip.
set -e
cc -O3 -o m1f_pairs m1f_pairs.c -lm
echo "built m1f_pairs"
