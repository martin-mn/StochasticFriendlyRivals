#!/bin/bash
#SBATCH --job-name=PRFull.stratcert.full
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=4-63
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/full_%A_%a.out
#SBATCH --error=logs/full_%A_%a.err
# FULL RUN: tasks 4-63 of the 64-task partition uid % 64 (12748 units = 7726 sign + and 5022 sign -).
# Tasks 0-3 were run by the probe 48975432 into the same out/ (all 800 units certified, t_0..3.done written);
# a rerun of them would find every unit done and exit (idempotent), so they are left out of the array.
# Sizing from the probe (test partition, Sapphire Rapids, 8 cpus): worst task 3:36 wall, 1651 core-s,
# MaxRSS 0.79 GB per task.  Projected per task 195-203 s.  --time: 2 x worst = 7:12, raised to 15 min to cover
# a unit whose census points all fail (up to 24 certify_point calls of <= 16 s each, about 6.5 min) -- the
# probe never used that fallback.  --mem: 5 x the probe's MaxRSS.  A task killed at the limit resumes where
# it stopped when resubmitted with the same array index.
cd "$SLURM_SUBMIT_DIR"
mkdir -p out
./task.sh "$SLURM_ARRAY_TASK_ID" 64 out
