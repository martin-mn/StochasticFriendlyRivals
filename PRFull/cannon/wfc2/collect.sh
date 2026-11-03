#!/bin/bash
# collect.sh -- after an arm has left the queue: job states per task (sacct), copy the outputs back, check, analyse.
#     bash collect.sh <arm> [<arm> ...]         arm: probe A1 A2 K0 B0 E3 E5 E3F2 E5F2 TS
# For each arm: (1) needs JOBID_<arm> (else it says it does not know of a job); (2) squeue -j: still queued/running?
# (3) sacct -j per task: every line whose State is not COMPLETED is printed; (4) rsync of out/<arm>/ (out/probe/ for the
# probe) and of the Slurm logs; (5) the number of REP lines against the task file, the neutral files, empty .err files;
# (6) analyze.py (analyze.py --probe for the probe).  Results are reported only if (3) and (5) are clean.
set -uo pipefail
cd "$(dirname "$0")"
REMOTE=${REMOTE:-CL1/PRFull/cannon/wfc2}
ssh -o BatchMode=yes cannon true || { echo "the Cannon session is closed: ask Martin to open ssh cannon" >&2; exit 1; }
for arm in "$@"; do
  echo "================ $arm"
  [ -f "JOBID_$arm" ] || { echo "$arm: no JOBID_$arm here, so no job of this kit is known for it"; continue; }
  jid=$(cat "JOBID_$arm")
  q=$(ssh -o BatchMode=yes cannon "squeue -j $jid -h 2>/dev/null | wc -l")
  [ "$q" -eq 0 ] || { echo "$arm ($jid): $q tasks still pending or running; not collecting"; continue; }
  bad=$(ssh -o BatchMode=yes cannon "sacct -j $jid --format=JobID,State,ExitCode,Elapsed -X --noheader" | awk '$2 != "COMPLETED"')
  tf=tasks_${arm}.txt; [ "$arm" = probe ] && tf=tasks_probe.txt
  sub=$arm; [ "$arm" = probe ] && sub=probe
  [ "$arm" = B0 ] && sub="B0 B0K0"
  mkdir -p out
  for s in $sub; do rsync -a "cannon:$REMOTE/out/$s" out/; done
  rsync -a --include "${arm}_${jid}_*" --exclude '*' "cannon:$REMOTE/out/" out/
  exp=$(grep -v '^#' "$tf" | awk '{n+=$6} END{print n}')
  if [ "$arm" = probe ]; then got=$(cat out/probe/*/r_*.txt 2>/dev/null | grep -c '^REP '); else
    got=$(for s in $sub; do cat out/$s/r_*.txt; done 2>/dev/null | grep -c '^REP '); fi
  nerr=$(find out -maxdepth 1 -name "${arm}_${jid}_*.err" -size +0 | wc -l | tr -d ' ')
  echo "$arm ($jid): REP lines $got of $exp expected; non-empty .err files: $nerr"
  clean=1
  [ -z "$bad" ] || { echo "tasks not COMPLETED (investigate against their outputs before any rerun):"; echo "$bad"; clean=0; }
  [ "$got" = "$exp" ] || clean=0
  [ "$nerr" = 0 ] || { find out -maxdepth 1 -name "${arm}_${jid}_*.err" -size +0 | head -5 | xargs -I{} sh -c 'echo "--- {}"; head -3 {}'; clean=0; }
  if [ "$arm" = A1 ] || [ "$arm" = A2 ] || [ "$arm" = A3 ]; then
    nn=$(ls out/$arm/n_*.txt 2>/dev/null | wc -l | tr -d ' '); ng=$(grep -vc '^#' "$tf")
    echo "$arm: neutral files $nn of $ng"; [ "$nn" = "$ng" ] || clean=0
  fi
  [ $clean = 1 ] || { echo "$arm: NOT clean -- no results reported"; continue; }
  if [ "$arm" = probe ]; then python3 analyze.py --probe; else python3 analyze.py $sub; fi
done
