#!/bin/bash
#SBATCH --job-name=PRFull.stratcert.probe
#SBATCH --comment=PRFull
#SBATCH --partition=test
#SBATCH --array=0-3
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/probe_%A_%a.out
#SBATCH --error=logs/probe_%A_%a.err
# PROBE: tasks 0-3 of the 64-task partition (uid % 64), i.e. 4/64 of the full run, all four loci, dims 0-8.
cd "$SLURM_SUBMIT_DIR"
mkdir -p out
./task.sh "$SLURM_ARRAY_TASK_ID" 64 out
