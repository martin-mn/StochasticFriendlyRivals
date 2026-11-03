#!/bin/bash
#SBATCH --job-name=PRFull.bingap
#SBATCH --comment=PRFull
#SBATCH --partition=shared
#SBATCH --array=0-63
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M
#SBATCH --time=00:30:00
#SBATCH --output=logs/full_%a.out
#SBATCH --error=logs/full_%a.err
# full census: task t handles every focal genome p with p % 64 == t, co-players q >= p
export NTASKS=64 EXTRA=""
srun ./task.sh
