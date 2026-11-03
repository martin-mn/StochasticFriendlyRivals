#!/bin/bash
#SBATCH --job-name=PartnersRivalsM2f.ca1dim
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-383
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-12:00:00
#SBATCH --output=out/dim_%A_%a.out
#SBATCH --error=out/dim_%A_%a.err
# 6 (class,sign) combinations x 64 stripes = 384 tasks.
# combination index c = SLURM_ARRAY_TASK_ID / 64 ; stripe = SLURM_ARRAY_TASK_ID % 64
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; mkdir -p out
i=${SLURM_ARRAY_TASK_ID}; c=$((i/64)); s=$((i%64))
case $c in
  0) CLS=C; SGN=1  ;;   # wedge W
  1) CLS=C; SGN=-1 ;;   # wedge S
  2) CLS=A; SGN=1  ;;   # wedge N
  3) CLS=A; SGN=-1 ;;   # wedge E
  4) CLS=Z; SGN=1  ;;   # switch-line ray, T>S  (efficiency W_DD=0 only)
  5) CLS=Z; SGN=-1 ;;   # switch-line ray, T<S
esac
srun ./dimprobe 5 512 $SGN $CLS $s 64 > out/dim_${CLS}_${SGN}_${s}.txt
