#!/bin/bash
#SBATCH --job-name=PRFull.wfc2.probe
#SBATCH --comment=PRFull
#SBATCH --partition=test
#SBATCH --array=1-60
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-03:00:00
#SBATCH --output=out/probe_%A_%a.out
#SBATCH --error=out/probe_%A_%a.err
export TASKFILE=tasks_probe.txt TS_NT=2000 TS_NI=10 TS_NB=400
srun ./task.sh
