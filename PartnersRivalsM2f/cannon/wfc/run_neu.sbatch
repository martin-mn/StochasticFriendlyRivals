#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.wfc-neu
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=1-84
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-04:00:00
#SBATCH --output=out/neu_%A_%a.out
#SBATCH --error=out/neu_%A_%a.err
export TASKFILE=tasks_neu.txt
srun ./task_neu.sh
