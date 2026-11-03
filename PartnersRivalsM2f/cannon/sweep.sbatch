#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.ca1sw
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-255
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-12:00:00
#SBATCH --output=out/sw_%A_%a.out
#SBATCH --error=out/sw_%A_%a.err
export NSWEEP=4000000
srun ./task.sh sweep
