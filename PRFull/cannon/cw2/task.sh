#!/bin/bash
# PRFull cannon/cw2 -- run ONE task (one line of a task file), idempotently.
#
#   bash task.sh <taskfile> [k] [--check]     k defaults to $SLURM_ARRAY_TASK_ID
#
# Task-file columns (whitespace separated, '#' lines and blank lines ignored; k counts data lines from 1):
#   arm  exe  isl  ipt  irep0  irep1  seedbase  eps  itend  outdir
#     arm       E3 | E5 | TS | VAL | VALTS | P* (probe rows) | KA (kernel masses, see below)
#     exe       the binary in the kit (cw.x, cw_e5.x, cw_ts.x, cwp.x, ...)
#     isl       cw.f's task index: 3*ipt + ie (ie = eps slot: 1 -> EPS1, 2 -> 1e-3, 3 -> 1e-4)
#     irep0..irep1   the replicates run (1..10 = the whole cell)
#     seedbase  '-' = the binary's own (40000000, production's base); otherwise passed as the 4th argument.
#               Replicate irep's seed is seedbase + (isl-1)*10 + irep.
#     eps itend what the h-file must say (checked: a binary of the wrong build cannot pass for the arm)
#     outdir    relative to the kit; the program's files <isl>, w<isl>, h<isl> (and t<isl>) land there
#   If irep0..irep1 = 1..10 and seedbase = '-', the program is called exactly as in production: `exe <isl>`.
#   KA rows: arm KA, exe kappa.x|kappa_e5.x, isl = the eps slot ie, ipt = a games file (games71.txt),
#            irep0 = M draws per game, outdir; one kappa line k<isl> per game (isl = 3*ipt+ie).
#
# COMPLETE (cw.f tasks) = <isl> has irep1-irep0+1 lines AND h<isl> has '# total' and exactly 105 count lines
# AND the h-file names the task's eps and itend AND (TSERIES binaries) t<isl> has nblk 'TS' lines and one
# '# end irun' line per replicate.  A complete task is skipped; an incomplete one is rerun from scratch (the
# program truncates its own files); after the run the same test is applied and a failure exits non-zero, so
# sacct shows FAILED rather than COMPLETED.  --check only reports (exit 0 complete, 1 not).
set -u
tf=${1:?usage: task.sh <taskfile> [k] [--check]}
k=${2:-${SLURM_ARRAY_TASK_ID:-}}
[ -n "$k" ] || { echo "no task index (argument 2 or SLURM_ARRAY_TASK_ID)"; exit 2; }
chk=${3:-}
kit=$(cd "$(dirname "$0")" && pwd)
line=$(grep -v '^[[:space:]]*#' "$kit/$tf" | grep -v '^[[:space:]]*$' | sed -n "${k}p")
[ -n "$line" ] || { echo "no task $k in $tf"; exit 2; }
line=${line%%#*}                                  # a trailing '# ...' is a comment
read -r arm exe isl ipt irep0 irep1 seedb eps itend outd extra <<< "$line"
[ -z "${extra:-}" ] || { echo "task $k: too many columns"; exit 2; }
[ "$chk" = --check ] || [ -x "$kit/$exe" ] || { echo "$kit/$exe not built -- run: bash compile.sh (or compile.sh probe)"; exit 2; }

if [ "$chk" = --check ] && [ ! -d "$kit/$outd" ]; then echo "task $k ($arm isl $isl) NOT complete: no $outd"; exit 1; fi
mkdir -p "$kit/$outd"
cd "$kit/$outd"
if [ ! -f games.dat ]; then
  # [verifier fix] --check must not write; and the 71 E3/E5 tasks (and the probe's shared dirs) start together in
  # ONE directory: a shared temp name let one task rename another's half-written copy into place.  A temp name
  # unique per host/job/process, then an atomic rename, means games.dat is always a complete copy.
  [ "$chk" = --check ] && { echo "task $k ($arm isl $isl) NOT complete: no games.dat in $outd"; exit 1; }
  tmp=".games.dat.new.$(hostname -s 2>/dev/null || echo h).${SLURM_JOB_ID:-0}.$$"
  cp "$kit/games.dat" "$tmp" && mv -f "$tmp" games.dat
fi
cmp -s games.dat "$kit/games.dat" || { echo "$outd/games.dat differs from the kit's"; exit 2; }

# ------------------------------------------------------------------ KA: kernel masses at the arm's eps
if [ "$arm" = KA ]; then
  ie=$isl; gf=$ipt; M=$irep0
  ipts=$(grep -v '^#' "$kit/$gf" | awk '{print $1}')
  n=0; nd=0
  for p in $ipts; do
    s=$((3*p+ie))
    if [ -s "k$s" ] && [ "$(wc -w < "k$s")" -eq 98 ]; then nd=$((nd+1)); continue; fi
    [ "$chk" = --check ] && { echo "KA k$s missing"; exit 1; }
    "$kit/$exe" "$s" "$M" > ".k$s.new" && mv -f ".k$s.new" "k$s" || { echo "kappa $s failed"; exit 3; }
    n=$((n+1))
  done
  [ "$chk" = --check ] && { echo "KA complete ($nd cells)"; exit 0; }
  echo "KA task $k: $n computed, $nd already present, in $outd"
  exit 0
fi

# ------------------------------------------------------------------ cw.f tasks
nrep=$((irep1-irep0+1))
case $eps in
  1e-2) es=0.100E-01 ;; 1e-3) es=0.100E-02 ;; 1e-4) es=0.100E-03 ;; 1e-5) es=0.100E-04 ;;
  *) echo "task $k: eps $eps not in the table"; exit 2 ;;
esac
args="$isl"
if [ "$irep0" != 1 ] || [ "$irep1" != 10 ] || [ "$seedb" != - ]; then
  args="$isl $irep0 $irep1"
  [ "$seedb" != - ] && args="$args $seedb"
fi
ts=0; case $exe in *ts*.x) ts=1 ;; esac      # [verifier fix] *_ts* missed the local test build cw_tts.x

complete () {
  local lines=0 nc=0 e it nb
  [ -f "$isl" ] && lines=$(wc -l < "$isl")
  [ "$lines" -eq "$nrep" ] || return 1
  [ -f "h$isl" ] && grep -q '^# total' "h$isl" || return 1
  nc=$(grep -vc '^#' "h$isl") || nc=0
  [ "$nc" -eq 105 ] || return 1
  e=$(awk '/^# cell:/{for(i=1;i<=NF;i++) if($i=="eps") print $(i+1)}' "h$isl")
  it=$(awk '/^# run:/{for(i=1;i<=NF;i++) if($i=="itend") print $(i+1)}' "h$isl")
  [ "$e" = "$es" ] && [ "$it" = "$itend" ] || { echo "task $k: h$isl says eps $e itend $it, the task wants $es $itend"; return 1; }
  if [ $ts -eq 1 ]; then
    [ -f "t$isl" ] || return 1
    nb=$(awk '/^# ts:/{for(i=1;i<=NF;i++) if($i=="nblk") print $(i+1)}' "t$isl")
    [ -n "$nb" ] || return 1
    [ "$(grep -c '^TS' "t$isl")" -eq $((nb*nrep)) ] || return 1
    [ "$(grep -c '^# end irun' "t$isl")" -eq "$nrep" ] || return 1
  fi
  return 0
}

if complete; then
  [ "$chk" = --check ] && { echo "task $k ($arm isl $isl reps $irep0-$irep1) complete in $outd"; exit 0; }
  echo "task $k ($arm isl $isl reps $irep0-$irep1) already complete in $outd, skipping"
  exit 0
fi
[ "$chk" = --check ] && { echo "task $k ($arm isl $isl reps $irep0-$irep1) NOT complete in $outd"; exit 1; }
t0=$(date +%s)
"$kit/$exe" $args
rc=$?
t1=$(date +%s)
if [ $rc -ne 0 ]; then echo "task $k: $exe $args exited $rc"; exit $rc; fi
if ! complete; then echo "task $k: $exe $args returned 0 but the products are INCOMPLETE"; exit 3; fi
echo "task $k ($arm $exe $args) complete in $outd, $((t1-t0)) s wall, $nrep replicate(s)"
