#!/bin/bash
# task.sh -- one array element of the memory-two continuum census.
#   $1 = mode (strata | sweep)
# SLURM_ARRAY_TASK_ID selects the stripe; NSTRIPES, KEFF, KOTH, NSWEEP come from the environment.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p out
i=${SLURM_ARRAY_TASK_ID}
case "$1" in
  strata)
    ./census strata "$i" "${NSTRIPES}" "${KEFF}" out/st_${i}.txt out/pt_${i}.txt "${KOTH}" ;;
  sweep)
    ./census sweep "${NSWEEP}" $((i+1000)) out/sw_${i}.txt ;;
  *) echo "unknown mode $1" >&2; exit 2 ;;
esac
