#!/bin/bash
#SBATCH --job-name=PRFull.wfc2.TS
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=1-240
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-02:00:00
#SBATCH --output=out/TS_%A_%a.out
#SBATCH --error=out/TS_%A_%a.err
export TASKFILE=tasks_TS.txt TS_NT=2000 TS_NI=10 TS_NB=400
srun ./task.sh
