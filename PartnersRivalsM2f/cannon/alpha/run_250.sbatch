#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.alpha250
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-255
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-12:00:00
#SBATCH --output=out/j_250_%A_%a.out
#SBATCH --error=out/j_250_%A_%a.err
export AM=250 NSWEEP=5860000
srun ./task.sh
