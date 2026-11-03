#!/bin/bash
#SBATCH --job-name=M2f.t2nprobe
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=test
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-00:40:00
#SBATCH --output=/dev/null
set -e
cd "$SLURM_SUBMIT_DIR"; mkdir -p out err
id=$SLURM_ARRAY_TASK_ID
lo=$(( id * 3587226 )); hi=$(( lo + 20000 ))
srun ./t2 ncensus 1 $lo $hi 64 4000 20260909 > out/nprobe_$id.txt 2> err/nprobe_$id.err
