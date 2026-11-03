#!/bin/bash
#SBATCH --job-name=PRFull.wfc2.K0
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=1-71
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-08:00:00
#SBATCH --output=out/K0_%A_%a.out
#SBATCH --error=out/K0_%A_%a.err
export TASKFILE=tasks_K0.txt
srun ./task.sh
