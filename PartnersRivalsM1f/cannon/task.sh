#!/bin/bash
# One array element.  SLURM_ARRAY_TASK_ID in 0..(4*SHARDS-1) selects a segment and a shard.
#   segment 0  m1     exhaustive grid of spacing 1/N1 vs the 16 memory-one vertices
#   segment 1  rand   random points drawn from all 81 strata, off-lattice, vs the 16
#   segment 2  m2     grid of spacing 1/N3 vs ALL 65536 memory-two co-players
#   segment 3  rand2  random points from all 81 strata vs all 65536 memory-two co-players
#   segment 4  fr     ITEM 2: predicate-free friendly-rival search on a grid, per wedge
#   segment 5  frand  ITEM 2: the same on random points drawn from all 81 strata
set -e
SHARDS=${SHARDS:-128}
N1=${N1:-400}          # segment 0 grid resolution
N2=${N2:-1000000000}   # segment 1 sample size, total across shards
N3=${N3:-20}           # segment 2 grid resolution
N4=${N4:-1000000}      # segment 3 sample size, total across shards
id=$SLURM_ARRAY_TASK_ID
seg=$(( id / SHARDS ))
shard=$(( id % SHARDS ))
case $seg in
  0) mode=m1;    n=$N1 ;;
  1) mode=rand;  n=$N2 ;;
  2) mode=m2;    n=$N3 ;;
  3) mode=rand2; n=$N4 ;;
  4) mode=fr;    n=$N5 ;;
  5) mode=frand; n=$N6 ;;
  *) echo "bad segment $seg" >&2; exit 2 ;;
esac
echo "# task $id  segment $seg  mode $mode  n $n  shard $shard/$SHARDS  host $(hostname)"
/usr/bin/time -v ./m1f_pairs "$mode" "$n" "$shard" "$SHARDS"
