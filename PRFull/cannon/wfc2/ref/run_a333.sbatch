#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.wfc-a333
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=1-512
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-06:00:00
#SBATCH --output=out/a333_%A_%a.out
#SBATCH --error=out/a333_%A_%a.err
export TASKFILE=tasks_a333.txt ITEND=10000000
srun ./task.sh
