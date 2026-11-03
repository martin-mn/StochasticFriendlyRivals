#!/bin/bash
#SBATCH --job-name=PartnersRivalsM1f.fr1probe
#SBATCH --comment=PartnersRivalsM1f
#SBATCH --partition=test
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-00:25:00
#SBATCH --output=out/frprobe_%A_%a.out
#SBATCH --error=out/frprobe_%A_%a.err
export SHARDS=1 N5=100 N6=20000000
id=$SLURM_ARRAY_TASK_ID
export SLURM_ARRAY_TASK_ID=$(( id + 4 ))
srun ./task.sh
