#!/bin/bash
#SBATCH --job-name=M2f.t2ncensus
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-NARR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-12:00:00
#SBATCH --output=/dev/null
export CHUNK=CHUNKV SEED=20260909
srun ./task.sh ncensusP
