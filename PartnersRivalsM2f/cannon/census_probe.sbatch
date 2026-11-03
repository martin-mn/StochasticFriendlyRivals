#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.ca1
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=test
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-02:00:00
#SBATCH --output=out/probe_%A_%a.out
#SBATCH --error=out/probe_%A_%a.err
export NSTRIPES=20000 KEFF=64 KOTH=8
srun ./task.sh strata
