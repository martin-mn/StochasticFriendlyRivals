#! /bin/bash
# cw / probe2: the RATE PROBE of run F2 of the continuous memory-one WF run (cw2.f, kit cannon/cw).
#
#     N = 100, m = 1, beta = 3, u = 1e-4, eps = 1e-4, 10 replicates -- but itend = 1e7, 1/10 of F2's 1e8
#     seeds 42000000 + (isl-1)*10 + irep;  cwp2.x = cw2.f -DITEND=10000000 + m1cls.c
#
# cw2.f is GENERATED from cw.f by mkcw2.py (diff cw.f cw2.f = the five sf1 -> sf2 parameter lines, the w-file maxc string, the two
# seed bases, and comments).  The probe is its own build (cwp2.x) in its own directory probe2/, for the reasons
# probe.slurm gives: cw2.x stays exactly the F2 production program, and pack.py refuses a probe cell
# (itend != 1e8).  The cost per generation does not depend on itend (init = itend/2 in both; the lineage
# population equilibrates at ~1 live lineage within ~1e3 generations at N*u = 0.01), so F2 production wall
# = 10 x probe2 wall.
#
# THE CELLS: the same five as probe.slurm and DiskM1WF's mw probe2 (the binary F2 probe, 46973350):
#   isl    3  PD the centre,   1527 PD,   1530 HA,   1533 SD,   1536 SH the farthest of all
# The binary F2 run (mw prod2, 46984101) took 526 core-h for 512 cells, 56.6-71.7 min per full task on
# Cascade Lake (shared), so ~6-7 min for 1/10 of a task there.
#
# LOCAL estimate (gfortran -O2, Martin's Mac, single thread, cwt2.x at isl 3, 45, 9, itend 2e5 x 10):
# 0.27-0.29 s per 1e5 generations, so ~280 s per probe2 task (1e8 generations) here, ~47 min per full F2 task.
# NOT the Cannon rate -- that is what this probe measures.  -t 60 is a ~13x margin on the local number.
#
# Log names o_f2probe / e_f2probe, so that verify.sh probe's glob e_probe* never picks them up.
#SBATCH -p test
#SBATCH -J PartnersRivalsM1f.cwprobe2
#SBATCH --comment=PartnersRivalsM1f
#SBATCH -t 60
#SBATCH -c 1
#SBATCH --mem=400
#SBATCH -o log/o_f2probe_%a.txt
#SBATCH -e log/e_f2probe_%a.txt
#SBATCH --array=3,1527,1530,1533,1536

module load intel/25.2.1-fasrc01
bash task.sh "$SLURM_ARRAY_TASK_ID" probe2 cwp2.x
