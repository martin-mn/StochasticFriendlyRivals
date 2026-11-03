#!/bin/bash
# task.sh -- one array element of an alpha-sweep.  alpha (permille) and n come from the
# environment, so each alpha is its own array with its own JOBID and its own sizing.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p out
i=${SLURM_ARRAY_TASK_ID}
./census psweep "${NSWEEP}" $(( 100000 + 1000*AM + i )) "${AM}" out/pt_${AM}_${i}.txt \
        > out/a_${AM}_${i}.tsv 2> out/a_${AM}_${i}.err
