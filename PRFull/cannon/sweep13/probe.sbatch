#!/bin/bash
# Probe of sweep13 (partition test: at most 5 jobs per user, so 4 tasks).  tasks_probe.txt:
#   0  the archived probe 45830542 task 4 (alpha 333 permille, seed 90004, 1e5 draws) with THIS build:
#      must reproduce ref/probe45830542_task4.tsv bit for bit
#   1  the same draws through the exact-rational path, alpha = 333/1000: the same counts
#   2,3  cost at exactly alpha = 1/3, 2e4 draws each, seeds 3900000, 3900001 (never used by the full run)
#SBATCH --job-name=PRFull.sweep13-probe
#SBATCH --comment=PRFull
#SBATCH --partition=test
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-00:30:00
#SBATCH --output=out_probe/probe_%A_%a.out
#SBATCH --error=out_probe/probe_%A_%a.err
export TASKFILE=tasks_probe.txt OUTDIR=out_probe
srun ./task.sh
