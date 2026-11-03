#!/bin/bash
# submit.sh -- submit one arm of the wfc2 kit to Cannon under the JOBID protocol (~/.claude/CLAUDE.md), from the Mac:
#     bash submit.sh <arm> [--yes]        arm: probe A1 A2 K0 B0 E3 E5 E3F2 E5F2 TS  (run_<arm>.sbatch / probe.sbatch)
# Needs the open ControlMaster session to Cannon (ssh cannon) and the kit staged at $REMOTE (stage.sh).
#   1. refuses if JOBID_<arm> already exists (never resubmit an array because the queue looks empty)
#   2. for an array of more than 100 tasks, or once this kit's submitted arms pass 100 tasks, prints how many of
#      mnowak's jobs are pending/running and stops unless --yes
#   3. sbatch --parsable in the remote kit; writes the array job ID to JOBID_<arm> here and there
#   4. appends the ledger line  date  PRFull  wfc2.<arm>  jobid  ntasks  note  to ~/jobs.tsv on Cannon
set -euo pipefail
cd "$(dirname "$0")"
REMOTE=${REMOTE:-CL1/PRFull/cannon/wfc2}
arm=${1:?usage: submit.sh <arm> [--yes]}; yes=${2:-}
script=run_${arm}.sbatch; [ "$arm" = probe ] && script=probe.sbatch
[ -f "$script" ] || { echo "no $script" >&2; exit 1; }
[ ! -e "JOBID_$arm" ] || { echo "JOBID_$arm exists ($(cat JOBID_$arm)): this arm was submitted; check it with sacct" >&2; exit 1; }
ntask=$(sed -n 's/^#SBATCH --array=1-\([0-9]*\).*/\1/p' "$script")
ssh -o BatchMode=yes cannon true || { echo "the Cannon session is closed: ask Martin to open ssh cannon" >&2; exit 1; }
# fair share: the rule is per array (> ~100 tasks), but this kit's arms are one campaign (611 sapphire tasks in arrays
# of 16-240), so the check also applies once the kit's submitted arms (JOBID_* here, the probe aside) plus this one
# pass 100 tasks.  mnowak's count is printed for a human to read, never used to decide anything about a job.
prev=0
for f in JOBID_*; do
  [ -e "$f" ] || continue; a=${f#JOBID_}; [ "$a" = probe ] && continue
  if [ -f "run_$a.sbatch" ]; then prev=$((prev + $(sed -n 's/^#SBATCH --array=1-\([0-9]*\).*/\1/p' "run_$a.sbatch"))); fi
done
if { [ "$ntask" -gt 100 ] || [ $((prev + ntask)) -gt 100 ]; } && [ "$arm" != probe ] && [ "$yes" != --yes ]; then
  n=$(ssh -o BatchMode=yes cannon "squeue -u mnowak -h -t PENDING,RUNNING -o '%A' | wc -l")
  echo "$arm has $ntask tasks (this kit has already submitted $prev in other arms); mnowak has $n pending/running jobs" \
       "across all projects.  Report that number, and re-run with --yes once adding to the queue is agreed." >&2
  exit 3
fi
# PARTITION=sapphire bash submit.sh probe  overrides the script's partition (e.g. if test refuses a 55-task array)
jid=$(ssh -o BatchMode=yes cannon "cd ~/$REMOTE && sbatch --parsable ${PARTITION:+--partition=$PARTITION} $script")
[[ "$jid" =~ ^[0-9]+$ ]] || { echo "sbatch did not return a job ID: $jid" >&2; exit 1; }
echo "$jid" > "JOBID_$arm"
case $arm in
  probe) note="probe on test: bit-identical checks + per-arm cost";;  A1) note="alpha=1/10 F1 71 games";;  A3) note="alpha=2/5 F1 71 games";;
  A2) note="alpha=1/5 F1 71 games";;  K0) note="local kernel K0 F1 71 games";;  B0) note="beta=0 neutral check, global+K0, 8 games";;
  E3) note="eps=1e-3 F1 71 games";;  E5) note="eps=1e-5 F1 71 games";;  E3F2) note="eps=1e-3 F2 24 games";;
  E5F2) note="eps=1e-5 F2 24 games";;  TS) note="time series 4x length 24 games";;  *) note="$arm";;
esac
ssh -o BatchMode=yes cannon "cd ~/$REMOTE && echo $jid > JOBID_$arm && printf '%s\tPRFull\twfc2.%s\t%s\t%s\t%s\n' \
  \$(date +%F) $arm $jid $ntask 'Ref1 6-8: $note' >> ~/jobs.tsv"
echo "submitted $arm: array job $jid ($ntask tasks); JOBID_$arm written here and in ~/$REMOTE; ledger line appended"
