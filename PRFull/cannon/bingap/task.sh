#!/bin/bash
# one array element: task index = SLURM_ARRAY_TASK_ID, NTASKS and extra flags from the batch script
set -e
cd "$(dirname "$0")"
T=${NTASKS:-64}; t=${SLURM_ARRAY_TASK_ID:?}
mkdir -p out
./bingap "$t" "$T" --checksym $EXTRA > "out/task_${t}.txt" 2> "out/task_${t}.err"
