#! /bin/bash
# cw / probe: the RATE PROBE of the continuous memory-one WF run (cw.f, kit cannon/cw).
#
#     N = 1000, m = 1, beta = 100, u = 1e-2, eps = 1e-4, 10 replicates -- but itend = 1e6, 1/10 of production
#
# WHY A SEPARATE BUILD (cwp.x = cw.f with -DITEND=1000000), NOT A COMMAND-LINE itend:
#  * cw.x stays exactly the production program: no argument can shorten a production task by accident,
#    and the program's interface is sf1's (one argument, the task index).
#  * the probe writes into its OWN directory probe/, never run1/.  Its h-files say itend 1000000, so
#    neither the idempotence guard of prod.slurm (which only counts lines) nor verify.sh could ever mistake
#    a probe cell for a finished production cell, and pack.py refuses them outright (itend != 1e7).
#  * the cost per generation does not depend on itend: init = itend/2 in both, so the unsampled and the
#    sampled (classifying) halves are the same fractions of the run, and the lineage population
#    equilibrates within ~1e3 generations (meanlive ~ 28 by itend = 2e4).  Production wall = 10 x probe
#    wall, to the extent the per-generation rate is stationary -- which the probe's own per-replicate
#    timings (the o-log) let you check.
#
# THE CELLS: the same five as DiskM1WF's mw probe1, so the continuum rate can be read against the binary
# run's measured rate at the identical games (mw README: F1 42.0-69.4 min per full task):
#   isl    3  ipt   0  r =  0.12  PD   the centre
#   isl 1527  ipt 508  r = 30.72  PD
#   isl 1530  ipt 509  r = 32.71  HA
#   isl 1533  ipt 510  r = 35.13  SD
#   isl 1536  ipt 511  r = 38.16  SH   the farthest of all
# Partition test: 5 tasks is exactly its job cap, 12 h its limit.
#
# LOCAL estimate (gfortran -O2, Martin's Mac, single thread, measured 2026-09-24 at isl 3, 45, 9):
# 0.284-0.290 s per 1e4 generations, so ~290 s per probe task here.  NOT the Cannon rate -- that is
# what this probe measures.  -t 60 is a 12x margin on the local number.
#SBATCH -p test
#SBATCH -J PartnersRivalsM1f.cwprobe
#SBATCH --comment=PartnersRivalsM1f
#SBATCH -t 60
#SBATCH -c 1
#SBATCH --mem=400
#SBATCH -o log/o_probe_%a.txt
#SBATCH -e log/e_probe_%a.txt
#SBATCH --array=3,1527,1530,1533,1536

module load intel/25.2.1-fasrc01
bash task.sh "$SLURM_ARRAY_TASK_ID" probe cwp.x
