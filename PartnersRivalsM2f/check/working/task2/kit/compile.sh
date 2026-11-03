#!/bin/bash
set -e
cc -O3 -march=native -o t2 t2.c m2f.c -lm
echo "built t2"
