#!/bin/bash
#SBATCH --job-name=M2f.t2seed
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-95
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-04:00:00
#SBATCH --output=/dev/null
set -e
cd "$SLURM_SUBMIT_DIR"; mkdir -p out err
id=$SLURM_ARRAY_TASK_ID
# 96 x 500 = 48000 emulated per-stratum seeds, on the very stratum the census missed.
# Sign +1, stratum 10561381 (wedge N, dim 7, true hit fraction ~0.107).
srun ./t2 seedtest 10561381 1 $(( 1000000 + id*500 )) 500 64 > out/seed_$id.txt 2> err/seed_$id.err
