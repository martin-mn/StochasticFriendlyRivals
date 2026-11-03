#!/bin/bash
# Reproduce the upper-bound certificate from scratch (`./run_all.sh`; `./run_all.sh all` adds the by-product) (Martin's Mac: 3 threads, nice 10; about 20 min
# in all, no single step over 7 min).  Each step's output goes to results/.  PYTHON=... and CC=... select the
# interpreter and the C compiler (default python3 and cc).
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python3}; CC=${CC:-cc}
mkdir -p runs results
$CC -O2 -ffp-contract=off -fno-fast-math -o ub ub.c -lm -lpthread      # certeng.c = check/cert/certeng.c
cmp certeng.c ../cert/certeng.c
N="nice -n 10"
# 1. exact efficiency class of all 3^16 open strata (15 s)
$N ./ub classes cls.bin 3                                       | tee results/classes.out
$N $PY counts.py                                                | tee results/counts.out
# 2. the certificate: the 12 co-players, a ternary DP over all 3^16 strata per co-player (1 s)
rm -f runs/x.idx runs/x.qp runs/x.qm runs/x.log
$N ./ub cover cls.bin qstar.txt runs/x 3                        | tee results/cover12.out
$N $PY status.py runs/x                                         | tee results/status_cover12.out
# 3. independent re-check of every (stratum, sign, q): direct 2^d corner enumeration + interval centre (6 min)
$N $PY mktriples.py runs/x runs/tri_star.bin                    | tee results/triples.out
$N ./ub verify runs/tri_star.bin runs/ver_star 3                | tee results/verify_cover12.out
# 4. the 12 rows in exact rational arithmetic, compared byte by byte with the interval rows (80 s)
$N $PY exact_rows.py qstar.txt                                  | tee results/exact_rows.out
# 5. exact spot checks (interior margins, classes, binary pairs) and the census consistency (40 s)
ln -sf tri_star.bin runs/x_tri.bin
$N $PY exact_check.py runs/x 20000 2000 4000                    | tee results/exact_check.out
$N $PY census_check.py runs/x                                   | tee results/census_check_cover12.out
$N $PY export.py                                                | tee results/export.out
# ---------------------------------------------------------------- optional: the all-dimension by-product
# 6. the 1681-co-player pool (results/pool_c2_q.txt, found adaptively by greedy.py): what is NOT excluded
#    in W, N, E, the T>S ray and the face p6=p8=0, in every dimension (2 min)
if [ "$1" = "all" ]; then
  export UB_BOUNDS='{"W":0,"W68":0,"N":0,"E":0,"ray+":0}'
  rm -f runs/c2.idx runs/c2.qp runs/c2.qm runs/c2.log
  $N ./ub cover cls.bin results/pool_c2_q.txt runs/c2 3         | tee results/cover_c2.out
  $N $PY status.py runs/c2                                      | tee results/status_c2_alldims.out
  $N $PY census_check.py runs/c2                                | tee results/census_check_c2.out
  $N $PY analysis.py runs/c2                                    | tee results/analysis_c2.out
  $N $PY extra.py runs/c2                                       | tee results/extra_c2.out
  # 7. independent re-check of all 14.2 M exclusions, in two halves of about 6.5 min each
  $N $PY mktriples_all.py runs/c2 runs/tri_all
  for h in A B; do $N ./ub verify runs/tri_all_$h.bin runs/ver_all_$h 3; done | tee results/verify_alldims_c2.out
fi
