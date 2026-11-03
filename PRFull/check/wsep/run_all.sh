#!/bin/bash
# Reproduce the W-separation certificate end to end (3 cores, nice 10; about 5 min, plus about 8 min for the cover).
# PYTHON=... and CC=... select the interpreter and the C compiler (default python3 and cc; compile.sh reads CC too).
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python3}; CC=${CC:-cc}
./compile.sh
$CC -O2 -ffp-contract=off -fno-fast-math -o cover cover.c -lm -lpthread
mkdir -p out/extra out/cover
N="nice -n 10"
# negative control: the 364 known friendly-rival strata of W must all survive (none excluded)
$N ./wsep file:ctrl_W364.txt 3 ctrl_W364.out > ctrl_W364.summary
for c in c6 c8 c68 c9a c9b; do $N ./wsep $c 3 out/$c.tsv > out/$c.summary; done
$N ./wsep c0 3 out/c0.tsv 1 0 1 > out/c0.summary                    # p0 < 1: efficiency only
cd out/extra
for spec in "x4 1...[1*]..........." "x13 1............[1*].." "x5 1....1.........." "x15 1..............[1*]" \
            "f1p11 1.....0.0..[0*]...." "f1p14 1.....0.0.....[0*]." "f2p11 1.....1.1..1...." "f2p14 1.....1.1.....1."; do
  set -- $spec; $N ../../wsep "t:$2" 3 $1.tsv > $1.summary
done
cd ../..
$N $PY verify.py 3 c9a c9b c68 c8 c6
$N $PY verify.py 3 extra/x4 extra/x13 extra/x5 extra/x15 extra/f1p11 extra/f1p14 extra/f2p11 extra/f2p14
$N $PY lemma_check.py
$N $PY exact_corners.py 600 9
$N $PY corner_pairs.py
$N $PY exact_table.py
$PY stats.py
if [ "${COVER:-1}" = 1 ]; then
  $PY cover_prep.py
  $N ./cover out/cover/sep.q out/cover/sep.pat 3 out/cover/sep.mask > out/cover/sep.log
  $PY cover_greedy.py
fi
