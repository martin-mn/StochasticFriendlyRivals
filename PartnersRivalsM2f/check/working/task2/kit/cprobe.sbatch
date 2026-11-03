#!/bin/bash
#SBATCH --job-name=M2f.t2cprobe
#SBATCH --comment=PartnersRivalsM2f
#SBATCH --partition=sapphire
#SBATCH --array=0-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-00:40:00
#SBATCH --output=/dev/null
set -e
cd "$SLURM_SUBMIT_DIR"; mkdir -p out err
id=$SLURM_ARRAY_TASK_ID
# ids 0,1,2 = wedge W (sign +1, g < 3^14 after the d0=2,d15=0 filter is applied inside);
# ids 3,4,5 = wedge S (sign -1).  Three widely separated 40000-wide windows each.
if [ $id -lt 3 ]; then s=1; lo=$(( id * 4782969 + 1000000 )); else s=-1; lo=$(( (id-3) * 14348907 + 1000000 )); fi
hi=$(( lo + 40000 ))
srun ./t2b ccensus $s $lo $hi 64 4000 20260909 > out/cprobe_$id.txt 2> err/cprobe_$id.err
