#!/bin/bash
#SBATCH --job-name=PRFull.wfc2.E5F2
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-04:00:00
#SBATCH --output=out/E5F2_%A_%a.out
#SBATCH --error=out/E5F2_%A_%a.err
export TASKFILE=tasks_E5F2.txt
srun ./task.sh
