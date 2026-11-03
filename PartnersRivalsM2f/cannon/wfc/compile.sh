#!/bin/bash
# wfc  : the large-population run F1 (N = 1000, beta = 100, mu = 1e-2), the defaults of wfc.c
# wfc2 : the small-population run F2 (N = 100,  beta = 3,   mu = 1e-4), the same source
set -e
cc -O3 -march=native -o wfc  wfc.c m2f.c -lm
cc -O3 -march=native -DNPOP=100 -DBETAV=3.0 -DMUV=1e-4 -o wfc2 wfc.c m2f.c -lm
echo "built wfc (F1) and wfc2 (F2)"
