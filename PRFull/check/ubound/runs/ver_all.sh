#!/bin/bash
cd /Users/martin/Documents/CL1/PRFull/check/ubound
for h in A B; do nice -n 10 ./ub verify runs/tri_all_$h.bin runs/ver_all_$h 3 > runs/ver_all_$h.out 2>&1; done
