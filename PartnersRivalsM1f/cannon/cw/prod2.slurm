#! /bin/bash
# cw / prod2: run F2 -- the WF process of DiskM1WF's sf2.f over the CONTINUOUS memory-one strategy space,
# at the 512 Voronoi cells of BinM2Ev ca6's sunflower, eps = 1e-4.  The data for the N = 100 versions of
# Figures 4 and 5 (SI Figs 4 and 5).  Kit cannon/cw; README.md there, section F2.
#
#     N = 100, m = 1, beta = 3, u = 1e-4, itend = 1e8, second half sampled, 10 replicates per cell
#     seeds 42000000 + (isl-1)*10 + irep;  cw2.x = cw2.f (= mkcw2.py(cw.f)) + m1cls.c (the real classifier)
#
# Task decode is sf1's = sf2's: cell ipt is task 3*ipt+3, the stride-3 array 3-1536:3.  Idempotent through
# the UNCHANGED task.sh (10 replicate lines AND an h-file with '# total' and 105 count lines), in run2/.
#
# -t IS PROVISIONAL.  LOCAL estimate only (gfortran -O2 on the Mac, single thread): 0.27-0.29 s per 1e5
# generations -> ~2700-2900 s = 45-48 min per task (1e9 generations).  The binary F2 run (mw prod2) took
# 56.6-71.7 min per task on Cascade Lake.  -t 240 is ~5x the local estimate.  REPLACE IT with ~3x the probe2's
# worst task x 10 once `bash verify.sh probe2` has run.  --mem as prod.slurm (400; F2's tables are
# 100x smaller: maxc = 202).
# PARTITION sapphire, as prod.slurm: the probe runs on test (Sapphire Rapids), and the login node builds
# -xHost Sapphire Rapids code.  NEVER run it on shared (Cascade Lake).
#SBATCH -p sapphire
#SBATCH -J PartnersRivalsM1f.cwprod2
#SBATCH --comment=PartnersRivalsM1f
#SBATCH -t 150
#SBATCH -c 1
#SBATCH --mem=400
#SBATCH -o log/o_f2prod_%a.txt
#SBATCH -e log/e_f2prod_%a.txt
#SBATCH --array=3-1536:3

module load intel/25.2.1-fasrc01
bash task.sh "$SLURM_ARRAY_TASK_ID" run2 cw2.x
