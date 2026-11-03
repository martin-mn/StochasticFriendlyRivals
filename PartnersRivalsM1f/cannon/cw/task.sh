#!/bin/bash
# Run ONE cw cell (one task = one game at eps = 1e-4 = 10 replicates), idempotently.
#
#   bash task.sh <isl> [rundir] [exe]        defaults: rundir run_local, exe cw.x
#
#   locally      bash task.sh 3                      (cw.x: a full production task, ~50 min here)
#                bash task.sh 3 run_local cwt.x      (a short test build, see compile.sh)
#   on Cannon    the Slurm scripts call it:  bash task.sh $SLURM_ARRAY_TASK_ID run1 cw.x
#
# isl is sf1's task index: cell ipt (0..511) at eps = 1e-4 is isl = 3*ipt+3.  The run directory must
# hold games.dat (submit.sh stages it; locally it is copied here if absent).
#
# COMPLETE = the summary <isl> has 10 lines AND h<isl> has its '# total' line and exactly 105 count
# lines (8 atoms + 81 strata + 16 corners).  cw.f, like sf1.f, flushes <isl> after every replicate but
# writes the h-file tail in one burst at the end, so 10 summary lines alone do NOT mean done: a task
# killed in that window would otherwise be skipped for ever (dw1's lesson, kept from prod1.slurm).
# A complete cell is skipped; an incomplete one is rerun from scratch (the program truncates its own
# three files).  After the run the same test is applied and a failure exits non-zero, so sacct shows
# the task FAILED rather than COMPLETED.
set -u
isl=${1:?usage: task.sh <isl> [rundir] [exe]}
rd=${2:-run_local}
exe=${3:-cw.x}
kit=$(cd "$(dirname "$0")" && pwd)
[ -x "$kit/$exe" ] || { echo "$kit/$exe not built -- run: bash compile.sh"; exit 2; }
mkdir -p "$rd"
if [ ! -f "$rd/games.dat" ]; then
  if [ -n "${SLURM_JOB_ID:-}" ]; then echo "$rd/games.dat missing -- use submit.sh"; exit 2; fi
  cp "$kit/games.dat" "$rd/.games.dat.new" && mv -f "$rd/.games.dat.new" "$rd/games.dat"
fi
cd "$rd"

complete () {
  local lines=0 nc=0
  [ -f "$isl" ] && lines=$(wc -l < "$isl")
  [ -f "h$isl" ] && grep -q '^# total' "h$isl" || return 1
  nc=$(grep -vc '^#' "h$isl") || nc=0
  [ "$lines" -eq 10 ] && [ "$nc" -eq 105 ]
}

if complete; then
  echo "task $isl already complete in $rd (10 replicate lines + 105-line h-file), skipping"
  exit 0
fi
t0=$(date +%s)
"$kit/$exe" "$isl"
rc=$?
t1=$(date +%s)
if [ $rc -ne 0 ]; then echo "task $isl: $exe exited $rc"; exit $rc; fi
if ! complete; then echo "task $isl: $exe returned 0 but the products are INCOMPLETE"; exit 3; fi
echo "task $isl complete in $rd, $((t1-t0)) s wall ($exe)"
