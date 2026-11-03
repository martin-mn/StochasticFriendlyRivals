#!/bin/bash
# $1 = mode ; remaining args passed through.  Writes out/<mode>_<id>.txt and err/<mode>_<id>.err
set -e
cd "$SLURM_SUBMIT_DIR"
mkdir -p out err
id=${SLURM_ARRAY_TASK_ID:-0}
mode=$1; shift
case "$mode" in
  list)    ./t2 list list.txt "$id" "$NSHARD" "$SEED"  > out/list_$id.txt   2> err/list_$id.err ;;
  ncensusP) lo=$(( id * CHUNK )); hi=$(( lo + CHUNK ));
           ./t2 ncensus  1 "$lo" "$hi" 64 4000 "$SEED" > out/nP_$id.txt 2> err/nP_$id.err ;;
  ncensusM) lo=$(( id * CHUNK )); hi=$(( lo + CHUNK ));
           ./t2 ncensus -1 "$lo" "$hi" 64 4000 "$SEED" > out/nM_$id.txt 2> err/nM_$id.err ;;
  seedtest) ./t2 seedtest "$SG" "$SSIGN" $(( id * SNS )) "$SNS" "$SK" > out/seed_$id.txt 2> err/seed_$id.err ;;
esac
