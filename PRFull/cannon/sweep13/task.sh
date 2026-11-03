#!/bin/bash
# task.sh -- one array element of the sweep13 rarity sweep (census psweep).
#   Full run (run_13.sbatch): ALPHA (exact "1/3", or a permille "333"), NSWEEP draws, seed SEEDBASE + index;
#                             tag <alpha with / as o>_<index>, e.g. 1o3_17.
#   Probe (probe.sbatch):     TASKFILE; its line index+1 (comment lines skipped) is   alpha  n  seed  tag
# Writes, in OUTDIR (default out):
#   a_<tag>.tsv   the ALPHA line (the collector's input)
#   a_<tag>.log   "task <i> alpha <a> n <n> seed <s> host <h> start <epoch>", census's summary line(s), "end <epoch>"
#   pt_<tag>.txt  the friendly rivals of W found, one per line, in draw order
# Nothing goes to the Slurm stderr unless something fails, so j_*.err must be empty.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
OUTDIR=${OUTDIR:-out}
mkdir -p "$OUTDIR"
i=${SLURM_ARRAY_TASK_ID}
if [ -n "${TASKFILE:-}" ]; then
  line=$(grep -v '^#' "$TASKFILE" | sed -n "$((i + 1))p")
  [ -n "$line" ] || { echo "task $i: no line in $TASKFILE" >&2; exit 2; }
  read -r al n seed tag <<< "$line"
else
  al=$ALPHA; n=$NSWEEP; seed=$((SEEDBASE + i)); tag="${al//\//o}_${i}"
fi
log="$OUTDIR/a_${tag}.log"
echo "task $i alpha $al n $n seed $seed host $(hostname) start $(date +%s)" > "$log"
./census psweep "$n" "$seed" "$al" "$OUTDIR/pt_${tag}.txt" > "$OUTDIR/a_${tag}.tsv" 2>> "$log"
echo "end $(date +%s)" >> "$log"
