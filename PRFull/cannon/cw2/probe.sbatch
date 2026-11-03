#!/bin/bash
# PRFull cw2 PROBE on partition sapphire (13 tasks, tasks_probe.txt; test refuses arrays of more than 5 jobs per user):
#   1 VAL    cw.x      isl 3 rep 1, production seed, 1e7: REP line and w-line must equal ref/prod (production run1)
#   2 VALTS  cw_ts1.x  the same with -DTSERIES (itend 1e7): the same two lines, and the in-program block check
#   3 PE3    cw.x      E3 full replicate (1e7) at ipt 466, the slowest production game      -> E3 rate
#   4-5 PE3s cwp.x     E3 1e6 generations at ipt 0 and 480                                    -> E3 rate spread
#   6 PE5    cw_e5.x   E5 full replicate at ipt 466                                           -> E5 rate; eps 1e-5 build
#   7-8 PE5s cwp_e5.x  E5 1e6 at ipt 0 and 480
#   9 PTS    cw_ts.x   one full TS task (4e7) at ipt 466 = the TS arm's task (466, rep 1)     -> TS rate
#   10-11 PTSs cwp_ts.x TS 4e6 at ipt 0 and 184
#   12-13 KA kappa.x / kappa_e5.x: kernel masses at eps 1e-3 / 1e-5 at the 71 games (M = 1e7)
# Expected: ~1 core-hour; the longest task (PTS) ~23 min, so --time 1 h is a 2.6x margin.
# Then: bash verify.sh probe  (writes PROBE_OK when every check passes; submit.sh wants it for the arms).
#SBATCH --partition=sapphire
#SBATCH --job-name=PRFull.cw2.probe
#SBATCH --comment=PRFull
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=400M
#SBATCH --output=out/log/probe_%A_%a.out
#SBATCH --error=out/log/probe_%A_%a.err
#SBATCH --array=1-13

module load intel/25.2.1-fasrc01
cd "$SLURM_SUBMIT_DIR"
bash task.sh tasks_probe.txt "$SLURM_ARRAY_TASK_ID"
