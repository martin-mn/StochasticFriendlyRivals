#!/bin/bash
#SBATCH --job-name=M2f.t2list
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-63
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-06:00:00
#SBATCH --output=/dev/null
export NSHARD=64 SEED=20260909
srun ./task.sh list
