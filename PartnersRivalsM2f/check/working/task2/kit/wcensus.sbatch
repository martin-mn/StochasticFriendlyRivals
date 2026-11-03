#!/bin/bash
#SBATCH --job-name=M2f.t2wcen
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-127
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-08:00:00
#SBATCH --output=/dev/null
set -e
cd "$SLURM_SUBMIT_DIR"; mkdir -p out err
id=$SLURM_ARRAY_TASK_ID; C=336303; TOT=43046721
lo=$(( id*C )); hi=$(( lo + C )); [ $hi -gt $TOT ] && hi=$TOT
[ $lo -ge $TOT ] && exit 0
srun ./t2b ccensus 1 $lo $hi 64 4000 20260909 > out/wc_$id.txt 2> err/wc_$id.err
