#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.alpha450
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-255
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-04:00:00
#SBATCH --output=out/j_450_%A_%a.out
#SBATCH --error=out/j_450_%A_%a.err
export AM=450 NSWEEP=390000
srun ./task.sh
