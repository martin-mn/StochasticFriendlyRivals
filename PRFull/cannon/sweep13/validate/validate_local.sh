#!/bin/bash
# validate_local.sh -- the local checks of README "Validation", on at most 2 cores at nice 10 (7 minutes on an M2 Pro).
#   bash validate/validate_local.sh [ARCHIVE_OUT]
# ARCHIVE_OUT: the archived 45831015/45830542 outputs, i.e. PartnersRivalsM2f_alpha/out extracted from
# PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz (extracted to a temporary directory if not given).
# Writes only under validate/build/.
set -euo pipefail
cd "$(dirname "$0")"
B=build; rm -rf "$B"; mkdir -p "$B"
CC=${CC:-cc}
A=${1:-}
if [ -z "$A" ]; then
  mkdir -p "$B/arch"; tar xzf /Users/martin/Documents/CL1/PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz \
    -C "$B/arch" 'PartnersRivalsM2f_alpha/out/*'; A=$B/arch/PartnersRivalsM2f_alpha/out
fi
$CC -O3 -o $B/census ../census.c ../m2f.c -lm
for t in replay kthr; do $CC -O3 -w -o $B/$t $t.c ../m2f.c -lm; done
$CC -O2 -o $B/invert invert.c
N="nice -n 10"
echo "== 1. probe 45830542 task 4 (seed 90004, 1e5 draws), permille and rational forms"
( $N $B/census psweep 100000 90004 333 $B/pt_p4.txt > $B/a_p4.tsv 2> $B/a_p4.err &
  $N $B/census psweep 100000 90004 333/1000 $B/pt_p4r.txt > $B/a_p4r.tsv 2> $B/a_p4r.err & wait )
cmp $B/a_p4.tsv ../ref/probe45830542_task4.tsv && cmp $B/a_p4.err ../ref/probe45830542_task4.err && echo "   333: ALPHA line and summary IDENTICAL"
diff <(cut -f3- $B/a_p4r.tsv) <(cut -f3- ../ref/probe45830542_task4.tsv) && echo "   333/1000: counts IDENTICAL"
echo "== 1b. the probe seed, recovered from the one W point of probe task 6 (alpha 0.4)"
$B/invert 400 3000000 $(cut -d' ' -f2- "$A/pt_6.txt")
echo "== 2. replay of the 512 task streams of 45831015 against the archived W points"
seq 0 511 | xargs -P2 -I{} sh -c "$N $B/replay 3190000 \$((433000+{})) 333 $A/pt_333_{}.txt" > $B/replay_all.txt || true
awk '{np+=$2; nf+=$3; if($2!=$3) bad++} END{print "   W points archived", np, "found", nf, "tasks with a miss", bad+0}' $B/replay_all.txt
( $N $B/census psweep 9000 433055 333 $B/pt_t55.txt > /dev/null 2>&1 & $N $B/census psweep 3000 433232 333 $B/pt_t232.txt > /dev/null 2>&1 & wait )
cmp $B/pt_t55.txt <(head -2 $A/pt_333_55.txt) && cmp $B/pt_t232.txt <(head -1 $A/pt_333_232.txt) && echo "   W-file prefixes of tasks 55 and 232 IDENTICAL"
echo "== 3. the law at exactly 1/3"
$N $B/census pdraw 1000000 3000000 1/3 | head -1
$N $B/census pdraw 100000000 3000000 1/3 | head -1
$B/kthr 1/3
echo "== 5. collect.py --reproduce"
( cd .. && $N python3 collect.py --reproduce | tail -1; rm -rf __pycache__ )
