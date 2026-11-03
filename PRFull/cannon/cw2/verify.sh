#!/bin/bash
# PRFull cannon/cw2 -- is an arm done, and done right?   (Cannon login node, in the kit directory)
#
#   bash verify.sh probe|E3|E5|TS|TSX
#
# Completion is decided PER TASK from sacct on the ID in JOBID_<arm>, never from squeue -u:
#   1. sacct: every task that is not COMPLETED is listed (an empty list means clean);
#   2. every non-empty stderr log of the array is listed;
#   3. task.sh --check on every task line: the products are complete and name the right eps and itend;
#   4. core-hours and the elapsed range from sacct.
# For the probe it also checks bit-identity against the archived production run (VAL, VALTS), the
# consistency of the kernel masses with production's (KA), prints the measured per-replicate rates and the
# -t each arm's sbatch should carry (2x the worst measured task), and writes PROBE_OK when every check passes.
set -u
cd "$(dirname "$0")"
arm=${1:?usage: verify.sh probe|E3|E5|TS|TSX}
case $arm in probe) tf=tasks_probe.txt ;; E3) tf=tasks_e3.txt ;; E5) tf=tasks_e5.txt ;; TS) tf=tasks_ts.txt ;; TSX) tf=tasks_tsx.txt ;;
  *) echo "unknown arm $arm"; exit 2 ;; esac
[ -f "JOBID_$arm" ] || { echo "no JOBID_$arm: this kit does not know of any $arm array (nothing inferred from squeue)"; exit 1; }
jid=$(cat "JOBID_$arm")
nt=$(grep -v '^[[:space:]]*#' "$tf" | grep -vc '^[[:space:]]*$')
ok=1
echo "== $arm: array $jid, $nt tasks ($tf)"
if command -v sacct >/dev/null 2>&1; then
  nq=$(squeue -j "$jid" -h 2>/dev/null | wc -l)
  echo "   queued or running now: $nq"
  sacct -j "$jid" --format=JobID,State,ExitCode,Elapsed -X --noheader > ".sacct_$arm" || true
  bad=$(awk '$2 != "COMPLETED"' ".sacct_$arm")
  ndone=$(awk '$2 == "COMPLETED"' ".sacct_$arm" | wc -l)
  echo "   sacct: $ndone of $nt tasks COMPLETED"
  [ -n "$bad" ] && { echo "   NOT COMPLETED:"; echo "$bad" | sed 's/^/      /'; ok=0; }
  [ "$ndone" -eq "$nt" ] || ok=0
  sacct -j "$jid" --format=ElapsedRaw,NCPUS -X --noheader |
    awk '{s+=$1*$2; if($1>mx)mx=$1; if(mn==""||$1<mn)mn=$1} END{printf "   core-hours %.2f, elapsed per task %d .. %d s\n", s/3600, mn, mx}'
else
  echo "   (no sacct here: steps 1 and 4 skipped -- run this on Cannon)"; ok=0
fi
ne=0
for f in out/log/${arm}_${jid}_*.err; do [ -s "$f" ] && { echo "   stderr not empty: $f"; ne=$((ne+1)); }; done
[ $ne -eq 0 ] && echo "   stderr: all empty" || ok=0
ninc=0
for k in $(seq 1 "$nt"); do bash task.sh "$tf" "$k" --check >/dev/null 2>&1 || { ninc=$((ninc+1)); [ $ninc -le 10 ] && bash task.sh "$tf" "$k" --check | sed 's/^/   /'; }; done
echo "   products: $((nt-ninc)) of $nt tasks complete"
[ $ninc -eq 0 ] || ok=0

if [ "$arm" = probe ]; then
  echo "== probe checks"
  pc=1
  for d in VAL VALTS; do
    if [ -f "out/probe/$d/3" ]; then
      if [ "$(sed -n 1p out/probe/$d/3)" = "$(sed -n 1p ref/prod/3)" ] &&
         [ "$(sed -n 3p out/probe/$d/w3)" = "$(sed -n 3p ref/prod/w3)" ]; then
        echo "   $d: REP line and w-line of isl 3 replicate 1 IDENTICAL to production (ref/prod)"
      else echo "   $d: DIFFERS from production:"; diff <(sed -n 1p out/probe/$d/3) <(sed -n 1p ref/prod/3); pc=0; fi
    else echo "   $d: no output"; pc=0; fi
  done
  if [ -f out/probe/VALTS/t3 ]; then
    echo "   VALTS: $(grep -c '^TS' out/probe/VALTS/t3) TS lines (the in-program check blocks == atom counts passed)"
  fi
  echo "   kernel masses (KA): max |kappa(eps) - kappa_e4| over atoms and games, in units of its binomial s.e.:"
  python3 -B collect.py --kappa-check out/probe/KA3 out/probe/KA5 || pc=0
  echo "   rates (wall per replicate, from the task logs) and the --time each arm needs:"
  python3 -B collect.py --probe-times out/log "probe_${jid}_" || true
  if [ $ok -eq 1 ] && [ $pc -eq 1 ]; then date > PROBE_OK; echo "== PROBE_OK written"; else rm -f PROBE_OK; echo "== probe NOT ok"; fi
fi
case $arm in   # the probe ran replicate 1 of one task of each arm with the arm's own seed: it must reappear
  E3) pf=out/probe/E3/1400; af=out/E3/1400 ;;
  E5) pf=out/probe/E5/1399; af=out/E5/1399 ;;
  TS) pf=out/probe/TS/1401; af=out/TS/p466_r01/1401 ;;
  *)  pf= ;;
esac
if [ -n "$pf" ] && [ -f "$pf" ] && [ -f "$af" ]; then
  if [ "$(sed -n 1p "$pf")" = "$(sed -n 1p "$af")" ]; then
    echo "   replicate 1 at ipt 466: IDENTICAL to the probe's run of it ($pf)"
  else echo "   replicate 1 at ipt 466 DIFFERS from the probe's run of it ($pf)"; ok=0; fi
  if [ "$arm" = TS ]; then cmp -s out/probe/TS/t1401 out/TS/p466_r01/t1401 && echo "   and its t-file too" ||
    { echo "   its t-file DIFFERS from the probe's"; ok=0; }; fi
fi
[ $ok -eq 1 ] && echo "== $arm clean" || echo "== $arm NOT clean"
