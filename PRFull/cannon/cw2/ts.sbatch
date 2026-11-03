#!/bin/bash
# PRFull cw2 arm TS: the time series over the memory-one cube, F1 regime, eps = 1e-4, 4x production length
# (itend 4e7, sampled half (2e7,4e7]), 24 games (games24.txt) x 10 replicates, ONE replicate per task (240 tasks),
# tasks_ts.txt.  cw_ts.x = cwx.f -DTSERIES -DITEND=40000000: 2000 blocks of 2e4 generations, exact per-block atom
# counts (lineages classified at birth), block efficiency, live lineages, mutants -> t<isl>.  No extra random
# number: the trajectory is the default build's (validated bit for bit against production at itend 1e7).
# Seeds 43000000 + (isl-1)*10 + irep: disjoint from production.
# --time: production F1 replicates at these 24 games took 258-309 s per 1e7 generations on sapphire (o-logs of
# 48366968); -DTSERIES costs +9.1% (measured locally, 334.6 vs 306.7 s per 1e7 at ipt 0); so 4e7 generations
# take <= 4 x 309 x 1.091 = 1350 s; 2x = 45 min.  Check against the probe's PTS task (ipt 466) before submitting.
# 240 tasks > 100: submit.sh prints the queue count and wants QUEUE_OK=1 (fair-share courtesy check).
#SBATCH --partition=sapphire
#SBATCH --job-name=PRFull.cw2.TS
#SBATCH --comment=PRFull
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:50:00
#SBATCH --mem=400M
#SBATCH --output=out/log/TS_%A_%a.out
#SBATCH --error=out/log/TS_%A_%a.err
#SBATCH --array=1-240

module load intel/25.2.1-fasrc01
cd "$SLURM_SUBMIT_DIR"
bash task.sh tasks_ts.txt "$SLURM_ARRAY_TASK_ID"
