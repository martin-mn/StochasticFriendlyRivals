#!/bin/bash
# PRFull cannon/cw2 -- submit ONE arm, following the Cannon job protocol (~/.claude/CLAUDE.md):
#
#   bash submit.sh probe|E3|E5|TS|TSX        (run on a Cannon login node, in the kit directory; TSX optional)
#
#  * sbatch --parsable, the array job ID written to JOBID_<arm> BEFORE anything else, and one line
#      date  PRFull  cw2  jobid  ntasks  [arm] note
#    appended to ~/jobs.tsv in the same step;
#  * refuses if JOBID_<arm> exists (a campaign is identified by its JOBID file; RESUBMIT=1 keeps the old
#    one as JOBID_<arm>.prev -- only after `bash verify.sh <arm>` has shown why);
#  * refuses if the arm's binaries are missing or were built from other sources (the h-file of every task
#    is checked by task.sh anyway);
#  * for arrays of more than 100 tasks, prints the user's queued+running count and refuses unless QUEUE_OK=1
#    (the fair-share courtesy check: report the number to Martin and ask before adding to a big queue);
#  * refuses an arm before its probe has been verified (PROBE_OK file, written by `verify.sh probe`
#    when every probe check passes), unless FORCE_NOPROBE=1.
set -eu
cd "$(dirname "$0")"
arm=${1:?usage: submit.sh probe|E3|E5|TS|TSX}
case $arm in
  probe) sb=probe.sbatch; tf=tasks_probe.txt; bins="cw.x cw_ts1.x cwp.x cw_e5.x cwp_e5.x cw_ts.x cwp_ts.x kappa.x kappa_e5.x"
         note="probe: VAL/VALTS bit-identity vs production isl 3 rep 1; per-replicate rates E3/E5/TS; kappa at 1e-3,1e-5 (71 games)" ;;
  E3)    sb=e3.sbatch; tf=tasks_e3.txt; bins="cw.x"
         note="eps=1e-3 F1 N=1000 beta=100 u=1e-2 1e7 gen 10 reps x 71 games; cw.x unchanged, isl=3ipt+2, seeds 40000000+" ;;
  E5)    sb=e5.sbatch; tf=tasks_e5.txt; bins="cw_e5.x"
         note="eps=1e-5 F1 10 reps x 71 games; cw_e5.x (-DEPS1=1.0d-5), isl=3ipt+1, seeds 40000000+" ;;
  TS)    sb=ts.sbatch; tf=tasks_ts.txt; bins="cw_ts.x"
         note="time series F1 eps=1e-4 4e7 gen, 2000 blocks, 24 games x 10 reps, 1 rep/task, seeds 43000000+" ;;
  TSX)   sb=tsx.sbatch; tf=tasks_tsx.txt; bins="cw_ts.x"
         note="OPTIONAL time series at 4 memory-one split games (246 75 311 185) x 10 reps, as TS, seeds 43000000+" ;;
  *) echo "unknown arm $arm"; exit 2 ;;
esac

python3 -B mktasks.py --check || { echo "REFUSED: task files differ from mktasks.py's output"; exit 1; }
for b in $bins; do [ -x "$b" ] || { echo "REFUSED: $b not built (bash compile.sh / compile.sh probe)"; exit 1; }; done
for b in $bins; do case $b in kappa*) ;; *)
  strings "$b" | grep -q 'usage: cwx.x' || { echo "REFUSED: $b is not a cwx.f build"; exit 1; } ;; esac; done
if [ "$arm" != probe ] && [ ! -f PROBE_OK ] && [ -z "${FORCE_NOPROBE:-}" ]; then
  echo "REFUSED: no PROBE_OK -- submit and verify the probe first (bash submit.sh probe; bash verify.sh probe)"; exit 1
fi
if [ -f "JOBID_$arm" ]; then
  if [ -z "${RESUBMIT:-}" ]; then
    echo "REFUSED: JOBID_$arm exists ($(cat JOBID_$arm)); check it with: bash verify.sh $arm"; exit 1
  fi
  mv -f "JOBID_$arm" "JOBID_$arm.prev"
fi

arr=$(grep '^#SBATCH --array=' "$sb" | sed 's/.*--array=//; s/%.*//')
lo=${arr%-*}; hi=${arr#*-}; ntasks=$((hi-lo+1))
nt=$(grep -v '^[[:space:]]*#' "$tf" | grep -vc '^[[:space:]]*$')
[ "$nt" -eq "$ntasks" ] || { echo "REFUSED: $sb --array=$arr ($ntasks) but $tf has $nt tasks"; exit 1; }

if [ "$ntasks" -gt 100 ]; then
  q=$(squeue -u "$USER" -h -t PENDING,RUNNING -o '%A' | wc -l)
  echo "fair share: $q of $USER's jobs are queued or running now (all projects)."
  if [ -z "${QUEUE_OK:-}" ]; then
    echo "REFUSED: $ntasks tasks; report the number above and resubmit with QUEUE_OK=1 once agreed"; exit 1
  fi
fi

mkdir -p out/log
jid=$(sbatch --parsable "$sb")
jid=${jid%%;*}
echo "$jid" > "JOBID_$arm"
printf '%s\tPRFull\tcw2\t%s\t%s\t[%s] %s\n' "$(date +%F)" "$jid" "$ntasks" "$arm" "$note" >> ~/jobs.tsv
echo "submitted $arm: array job $jid ($ntasks tasks) -> JOBID_$arm, ledger ~/jobs.tsv"
