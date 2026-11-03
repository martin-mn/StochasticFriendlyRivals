#!/bin/bash
#SBATCH --job-name=PRFull.wfc2.B0
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=1-16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-08:00:00
#SBATCH --output=out/B0_%A_%a.out
#SBATCH --error=out/B0_%A_%a.err
export TASKFILE=tasks_B0.txt
srun ./task.sh
