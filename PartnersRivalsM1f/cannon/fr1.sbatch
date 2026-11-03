#!/bin/bash
#SBATCH --job-name=PartnersRivalsM1f.fr1
#SBATCH --comment=PartnersRivalsM1f
#SBATCH --partition=sapphire
#SBATCH --array=0-255
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-04:00:00
#SBATCH --output=out/fr1_%A_%a.out
#SBATCH --error=out/fr1_%A_%a.err
# segments 4 and 5 only; SHARDS=128 so task id 0..127 -> fr, 128..255 -> frand
export SHARDS=128 N5=500 N6=100000000000
id=$SLURM_ARRAY_TASK_ID
export SLURM_ARRAY_TASK_ID=$(( id + 4*128 ))
srun ./task.sh
