#!/bin/bash
# The rarity sweep at EXACTLY alpha = 1/3: 512 x 3 190 000 = 1 633 280 000 draws, the size of the alpha = 0.333
# run 45831015, with seeds 3000000 + i (that run: 433000 + i; its probe 45830542: 90000 + i).
#SBATCH --job-name=PRFull.sweep13
#SBATCH --comment=PRFull
#SBATCH --partition=sapphire
#SBATCH --array=0-511
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-08:00:00
#SBATCH --output=out/j_13_%A_%a.out
#SBATCH --error=out/j_13_%A_%a.err
export ALPHA=1/3 NSWEEP=3190000 SEEDBASE=3000000
srun ./task.sh
