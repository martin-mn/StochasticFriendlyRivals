#!/bin/bash
# Probe: two shards of each segment at reduced size, to measure the true per-unit rate on
# Cannon before the full array is sized.  Never rescale a local timing to guess this.
#SBATCH --job-name=PartnersRivalsM1f.riv1probe
#SBATCH --comment=PartnersRivalsM1f
#SBATCH --partition=test
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-00:30:00
#SBATCH --output=out/probe_%A_%a.out
#SBATCH --error=out/probe_%A_%a.err
export SHARDS=1 N1=100 N2=10000000 N3=8 N4=20000
srun ./task.sh
