#!/bin/bash
# submit.sh -- submit the probe or the full array of sweep13 under the JOBID protocol (~/.claude/CLAUDE.md):
#     bash submit.sh probe          probe.sbatch   -> JOBID_probe   (partition test, 4 tasks)
#     bash submit.sh 13 [--yes]     run_13.sbatch  -> JOBID_13      (sapphire, 512 tasks)
#   1. refuses if the JOBID file already exists (never resubmit because the queue looks empty: check sacct)
#   2. for more than 100 tasks, prints how many of mnowak's jobs are pending/running (for a human to read) and
#      stops unless --yes
#   3. sbatch --parsable in the remote kit; the array job ID goes to the JOBID file here and there
#   4. appends  date  PRFull  sweep13[.probe]  jobid  ntasks  note  to ~/jobs.tsv on Cannon
set -euo pipefail
cd "$(dirname "$0")"
REMOTE=${REMOTE:-CL1/PRFull/cannon/sweep13}
arm=${1:?usage: submit.sh probe|13 [--yes]}; yes=${2:-}
case $arm in
  probe) script=probe.sbatch; kit=sweep13.probe; note="probe on test: repro of 45830542 task 4 + cost at alpha=1/3";;
  13)    script=run_13.sbatch; kit=sweep13; note="rarity sweep at exactly alpha=1/3, 512x3.19e6 draws, seeds 3000000+i";;
  *) echo "usage: submit.sh probe|13 [--yes]" >&2; exit 1;;
esac
[ ! -e "JOBID_$arm" ] || { echo "JOBID_$arm exists ($(cat JOBID_$arm)): already submitted; check it with sacct" >&2; exit 1; }
range=$(sed -n 's/^#SBATCH --array=\([0-9]*\)-\([0-9]*\).*/\1 \2/p' "$script"); read -r a0 a1 <<< "$range"
ntask=$((a1 - a0 + 1))
ssh -o BatchMode=yes cannon true || { echo "the Cannon session is closed: ask Martin to open ssh cannon" >&2; exit 1; }
if [ "$ntask" -gt 100 ] && [ "$yes" != --yes ]; then
  n=$(ssh -o BatchMode=yes cannon "squeue -u mnowak -h -t PENDING,RUNNING -o '%A' | wc -l")
  echo "$script has $ntask tasks; mnowak has $n pending/running jobs across all projects." \
       "Report that number, and re-run with --yes once adding to the queue is agreed." >&2
  exit 3
fi
jid=$(ssh -o BatchMode=yes cannon "cd ~/$REMOTE && sbatch --parsable $script")
[[ "$jid" =~ ^[0-9]+$ ]] || { echo "sbatch did not return a job ID: $jid" >&2; exit 1; }
echo "$jid" > "JOBID_$arm"
ssh -o BatchMode=yes cannon "cd ~/$REMOTE && echo $jid > JOBID_$arm && printf '%s\tPRFull\t%s\t%s\t%s\t%s\n' \
  \$(date +%F) $kit $jid $ntask '$note' >> ~/jobs.tsv"
echo "submitted $arm: array job $jid ($ntask tasks); JOBID_$arm written here and in ~/$REMOTE; ledger line appended"
