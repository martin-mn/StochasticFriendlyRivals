#!/bin/bash
# task.sh -- one line of $TASKFILE (lines starting with # are skipped) per Slurm array element:
#     outdir  bin  alpha  game  rep0  nrep  itend  NT  NI  seedbase  neutralM  neutralseed
#   Runs replicates rep0..rep0+nrep-1 of <bin> (a binary of compile.sh) at line <game> of games.dat, replicate r seeded
#   seedbase + 10 (game-1) + r (as production), itend generations, NT sample times x NI individuals over the second half,
#   then, if neutralM > 0, neutralM kernel draws at <alpha> seeded neutralseed (the alpha-measure of the atoms).
#   alpha is exact ("1/4", "2/5", "1/3").  Every output goes to out/<outdir>/ and carries alpha in its name (1/4 -> 1o4):
#     r_<a>_<game>_<rep0>.txt   stdout: '#' header lines, REP lines (and MUT, or TB/TS/TSR, lines for those builds)
#     S_<a>_<game>_<rep0>.bin   every distinct sampled strategy (int32 irep, count, atom; double p[16])
#     T_<a>_<game>_<rep0>.bin[.k]  wfc_ts only: the time-series strategies and their per-individual indices
#     n_<a>_<game>.txt          the NEU line of the kernel draws (neutralM > 0 only)
#   A line of "done" goes to stdout, so a successful task leaves its Slurm .err file empty.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
k=${SLURM_ARRAY_TASK_ID:?}
line=$(grep -v '^#' "${TASKFILE:?}" | sed -n "${k}p")
read -r od bin am g r0 nr it nt ni sb nm ns <<< "$line"
[ -n "${ns:-}" ] || { echo "task $k: line '$line' of $TASKFILE has fewer than 12 fields" >&2; exit 2; }
at=${am//\//o}; tag="${at}_${g}_${r0}"
mkdir -p "out/$od"
if [ "$bin" = wfc_ts ]; then export TS_DUMP="out/$od/T_${tag}.bin"; fi
"./$bin" run "$g" "$am" "$it" "$nt" "$ni" "$nr" "$sb" "$r0" "out/$od/S_${tag}.bin" > "out/$od/r_${tag}.txt"
if [ "$nm" -gt 0 ]; then "./$bin" neutral "$g" "$am" "$nm" "$ns" > "out/$od/n_${at}_${g}.txt"; fi
echo "done $k $od $tag"
