#!/bin/bash
# task.sh -- one line of $TASKFILE per array element:   alpha  game  rep0  nrep  neutralM  seedbase
#   alpha is exact: "1/3", or permille ("500").  Runs replicates rep0..rep0+nrep-1 of the continuum Wright-Fisher
#   process at game line <game> of games.dat, then, if neutralM > 0, the alpha-measure of the seven atoms there.
#   Every output file carries alpha in its name (1/3 -> 1o3), so arrays at different alpha can share a directory.
#   Environment: ITEND (generations), and optionally WFC (the binary, default ./wfc = F1; ./wfc2 = F2),
#   OUTDIR (default out), NT and NI (sample times and individuals per time, default 1000 and 4).
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
WFC=${WFC:-./wfc}; OUTDIR=${OUTDIR:-out}; NT=${NT:-1000}; NI=${NI:-4}
mkdir -p "$OUTDIR"
k=${SLURM_ARRAY_TASK_ID}
read am g r0 nr nm sb < <(grep -v '^#' "$TASKFILE" | sed -n "${k}p")
at=${am//\//o}
tag="${at}_${g}_${r0}"
"$WFC" run "$g" "$am" "${ITEND}" "$NT" "$NI" "$nr" "$sb" "$r0" "$OUTDIR/S_${tag}.bin" > "$OUTDIR/r_${tag}.txt"
if [ "$nm" -gt 0 ]; then "$WFC" neutral "$g" "$am" "$nm" > "$OUTDIR/n_${at}_${g}.txt"; fi
echo "done $k $tag" >&2
