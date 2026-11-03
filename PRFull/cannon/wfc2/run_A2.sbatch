#!/bin/bash
#SBATCH --job-name=PRFull.wfc2.A2
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=1-71
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-04:00:00
#SBATCH --output=out/A2_%A_%a.out
#SBATCH --error=out/A2_%A_%a.err
export TASKFILE=tasks_A2.txt
srun ./task.sh
