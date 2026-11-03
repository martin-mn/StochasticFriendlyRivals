#!/bin/bash
# Cannon build. gcc is fine; -march=native is safe because every partition used here
# (sapphire / shared / test) compiles and runs on the same node.
set -e
cc -O3 -march=native -o census census.c m2f.c -lm
echo "built census"
