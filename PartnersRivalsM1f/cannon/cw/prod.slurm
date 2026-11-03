#! /bin/bash
# cw / prod: the WF process of DiskM1WF's sf1.f (run F1) over the CONTINUOUS memory-one strategy space,
# at the 512 Voronoi cells of BinM2Ev ca6's sunflower, eps = 1e-4.  Kit cannon/cw; README.md there.
#
#     N = 1000, m = 1, beta = 100, u = 1e-2, itend = 1e7, second half sampled, 10 replicates per cell
#     seeds 40000000 + (isl-1)*10 + irep;  cw.x = cw.f + m1cls.c (the real classifier)
#
# Task decode is sf1's: ie = mod(isl-1,3)+1 (eps 1e-2, 1e-3, 1e-4), ig = (isl-1)/3+1 = games.dat row,
# so cell ipt is task 3*ipt+3 and eps = 1e-4 is the stride-3 array 3-1536:3.
#
# Idempotent (task.sh): a cell counts as done only with 10 replicate lines AND an h-file carrying its
# '# total' line and exactly 105 count lines -- cw.f flushes the summary after every replicate but writes
# the h-file tail in one burst at the end (prod1.slurm's lesson, kept).  A finished cell is skipped, an
# unfinished one rerun from scratch; task.sh exits non-zero if the products are incomplete afterwards.
#
# -t and --mem ARE PROVISIONAL.  LOCAL estimate only (gfortran -O2 on the Mac, single thread): 0.29 s per
# 1e4 generations -> ~2900 s = 48 min per task (1e8 generations), m1cls 2.5-4 % of it.  sf1.x (the binary
# run) measured 0.228 s locally and 42-69 min per task on Cannon (shared).  -t 360 is 7.5x the local
# estimate.  REPLACE IT with ~3x the probe's worst task x 10 once `bash verify.sh probe` has run.
# Peak RSS measured locally 66 MB (the pair cache is 3 x 2002^2 x 8 B = 96 MB at most); --mem=400.
# PARTITION sapphire, not mw's shared: the probe runs on test, whose nodes are Sapphire Rapids, so the
# probe's rate is a measurement of THIS hardware only if production runs on Sapphire Rapids too (shared
# is Cascade Lake).  It also keeps -xHost safe: code tuned on a login node runs on the newer SPR cores.
#SBATCH -p sapphire
#SBATCH -J PartnersRivalsM1f.cwprod
#SBATCH --comment=PartnersRivalsM1f
#SBATCH -t 150
#SBATCH -c 1
#SBATCH --mem=400
#SBATCH -o log/o_prod_%a.txt
#SBATCH -e log/e_prod_%a.txt
#SBATCH --array=3-1536:3

module load intel/25.2.1-fasrc01
bash task.sh "$SLURM_ARRAY_TASK_ID" run1 cw.x
