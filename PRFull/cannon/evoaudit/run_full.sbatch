#!/bin/bash
# PRFull evoaudit FULL RUN: chunks 4-127 of the 128 chunks of sample/sample.npy, one core each (124 tasks, ~213
# strategies per chunk).  Chunks 0-3 are the probe's (JOBID_probe = 48975385, all COMPLETED, in out_probe/); the
# analysis merges out_probe and out:  python3 src/analyze.py out_probe,out results
# NOT SUBMITTED.  Submit only by the protocol in README.md (JOBID_full + one ~/jobs.tsv line).
# Cost from the probe (4 x 213 strategies on sapphire, one core): 1.39 s per strategy, 277-318 s per chunk
# (sacct Elapsed 4:40-5:22).  Projected: ~300 s per chunk, slowest of 124 ~6-7 min (bootstrap 99.9%: 422 s),
# ~10.3 core-hours in all.  --time = 2 x the worst probe task (2 x 5:22 = 10:44), rounded up.
# Idempotent: a timed-out or killed task is completed by resubmitting its index (records already written are skipped).
#SBATCH --job-name=PRFull.evoaudit.full
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=4-127
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:11:00
#SBATCH --output=logs/full_%A_%a.out
#SBATCH --error=logs/full_%A_%a.err
cd "$SLURM_SUBMIT_DIR"
bash task.sh "$SLURM_ARRAY_TASK_ID" 128 out
