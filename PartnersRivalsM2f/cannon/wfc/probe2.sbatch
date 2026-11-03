#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.wfc-probe2
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=1-28
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-02:00:00
#SBATCH --output=out/probe2_%A_%a.out
#SBATCH --error=out/probe2_%A_%a.err
export TASKFILE=tasks_probe2.txt ITEND=10000000
srun ./task.sh
