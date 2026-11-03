#!/bin/bash
#SBATCH --job-name=PartnersRivalsM1f.riv1
#SBATCH --comment=PartnersRivalsM1f
#SBATCH --partition=sapphire
#SBATCH --array=0-511
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-04:00:00
#SBATCH --output=out/riv1_%A_%a.out
#SBATCH --error=out/riv1_%A_%a.err
export SHARDS=128 N1=500 N2=10000000000 N3=30 N4=1000000
srun ./task.sh
