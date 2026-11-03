#!/bin/bash
# PRFull cw2 arm E5: eps = 1e-5 over the memory-one cube, F1 regime (N = 1000, beta = 100, u = 1e-2, itend 1e7,
# second half sampled), 10 replicates per task = one game per task, 71 games (games71.txt), tasks_e5.txt.
# cw_e5.x is cwx.f built with -DEPS1=1.0d-5: eps slot ie = 1 holds 1e-5 (isl = 3*ipt+1); nothing else changes.
# Seeds 40000000 + (isl-1)*10 + irep: disjoint from production (isl = 3*ipt+3) and from E3 (3*ipt+2).
# --time: production F1 tasks (same program, 10 x 1e7 generations) took 2311-3136 s on sapphire (o-logs of
# 48366968); 2 x 3136 s = 105 min.  If the probe's E5 full replicate (PE5, ipt 466) exceeds 314 s by more than
# 10%, raise --time to 2 x 10 x that replicate before submitting.  MaxRSS measured 34 MB; --mem 400M as production.
#SBATCH --partition=sapphire
#SBATCH --job-name=PRFull.cw2.E5
#SBATCH --comment=PRFull
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:45:00
#SBATCH --mem=400M
#SBATCH --output=out/log/E5_%A_%a.out
#SBATCH --error=out/log/E5_%A_%a.err
#SBATCH --array=1-71

module load intel/25.2.1-fasrc01
cd "$SLURM_SUBMIT_DIR"
bash task.sh tasks_e5.txt "$SLURM_ARRAY_TASK_ID"
