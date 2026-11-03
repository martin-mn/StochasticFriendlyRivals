#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.ca1
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-511
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-12:00:00
#SBATCH --output=out/st_%A_%a.out
#SBATCH --error=out/st_%A_%a.err
export NSTRIPES=512 KEFF=64 KOTH=8
srun ./task.sh strata
