#!/bin/bash
# rcheck.sh -- the replay of the efficiency test (README "Replay of the efficiency test"): build rcheck.c, replay the
# listed tasks of the run at alpha = 1/3 (seed 3000000 + task, 3 190 000 draws each), and compare with the deposit.
#   bash validate/rcheck.sh [TASK ...]        default: all 512 tasks, 0..511 (47 core-minutes on an M2 Pro; about
#                                             17 min at JOBS=3 for the x86_64 binary that cc builds when it is called
#                                             from an x86_64 process under Rosetta, 5.8 s of CPU per task)
#   JOBS=N replays N tasks at a time (default 2); CC=... and PYTHON=... select the compiler and the interpreter.
# Writes only under validate/build/: the program rcheck and replay_<task>.out per task (its RC line and one FLAG line per
# flagged draw), which rcheck_compare.py then compares with rc.tsv, flags.tsv and ../out_13/a_1o3_<task>.tsv.
set -euo pipefail
cd "$(dirname "$0")"
B=build; mkdir -p "$B"
CC=${CC:-cc}; PY=${PYTHON:-python3}; JOBS=${JOBS:-2}
T=("$@"); [ $# -gt 0 ] || T=($(seq 0 511))
$CC -O2 -ffp-contract=off -fno-fast-math -o $B/rcheck rcheck.c -lm
printf '%s\n' "${T[@]}" | xargs -P "$JOBS" -I{} sh -c "nice -n 10 $B/rcheck {} 3190000 classify > $B/replay_{}.out"
$PY rcheck_compare.py --replay $B "${T[@]}"
