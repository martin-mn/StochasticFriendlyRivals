#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.alpha333
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-511
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-08:00:00
#SBATCH --output=out/j_333_%A_%a.out
#SBATCH --error=out/j_333_%A_%a.err
export AM=333 NSWEEP=3190000
srun ./task.sh
