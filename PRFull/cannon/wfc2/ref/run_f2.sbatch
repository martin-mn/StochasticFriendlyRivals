#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.wfc2-f2
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=1-512
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-04:00:00
#SBATCH --output=out_f2/f2_%A_%a.out
#SBATCH --error=out_f2/f2_%A_%a.err
export TASKFILE=tasks_f2.txt ITEND=100000000 WFC=./wfc2 OUTDIR=out_f2 NT=4000 NI=1
srun ./task.sh
