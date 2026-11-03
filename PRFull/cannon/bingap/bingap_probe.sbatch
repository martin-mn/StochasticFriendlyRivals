#!/bin/bash
#SBATCH --job-name=PRFull.bingap
#SBATCH --comment=PRFull
#SBATCH --partition=test
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M
#SBATCH --time=00:10:00
#SBATCH --output=logs/probe_%a.out
#SBATCH --error=logs/probe_%a.err
# probe: tasks 0-3 of 64 (the test partition allows at most 5 jobs per user), restricted to p < 8192 (128 focal genomes each, ~7.9M pairs per task)
export NTASKS=64 EXTRA="--pmax 8192"
srun ./task.sh
