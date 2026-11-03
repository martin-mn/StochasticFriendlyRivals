#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.alpha-probe
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --output=out/probe_%A_%a.out
#SBATCH --error=out/probe_%A_%a.err
export NSWEEP=100000
srun ./task.sh
