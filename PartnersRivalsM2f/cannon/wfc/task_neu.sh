#!/bin/bash
# task_neu.sh -- one line of $TASKFILE per array element:   alpha  game  M  seed  chunk
#   a further M independent draws of the alpha-law at game line <game>, each classified exactly: refines the
#   alpha-measure of the atoms where the 20 000 draws of the main run left the most-enriched atom undetermined.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p out
k=${SLURM_ARRAY_TASK_ID}
read am g m sd ch < <(grep -v '^#' "$TASKFILE" | sed -n "${k}p")
at=${am//\//o}
./wfc neutral "$g" "$am" "$m" "$sd" > "out/nr_${at}_${g}_${ch}.txt"
echo "done $k ${at}_${g}_${ch}" >&2
