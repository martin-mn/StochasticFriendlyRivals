#!/bin/bash
# PRFull evoaudit PROBE: chunks 0-3 of the 128 chunks of sample/sample.npy (4 x ~213 of the 27 261 strategies), one core each.
# Measures the per-strategy cost on Cannon; the outputs are complete chunks of the full run (copy out_probe/ to out/ and
# the full run skips them).  Expected ~15 min per task; --time 1:30 is a 6x margin.
#SBATCH --job-name=PRFull.evoaudit.probe
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:30:00
#SBATCH --output=logs/probe_%A_%a.out
#SBATCH --error=logs/probe_%A_%a.err
cd "$SLURM_SUBMIT_DIR"
bash task.sh "$SLURM_ARRAY_TASK_ID" 128 out_probe
