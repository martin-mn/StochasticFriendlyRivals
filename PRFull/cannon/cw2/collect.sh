#!/bin/bash
# PRFull cannon/cw2 -- LOCAL: completeness of every task, then the tables (collect.py).
#
#   bash collect.sh                       analyse what is under out/ (E3, E5, TS as present)
#   bash collect.sh fetch <remote-kit>    first copy out/ and the JOBID files back from Cannon (read-only there;
#                                         only when Martin has said the arrays may be collected), then analyse
# Cluster facts (sacct COMPLETED per task, empty stderr) are established on Cannon by `bash verify.sh <arm>`;
# collect.py reports on the products only.  Tables: results/{E3,E5}_games.tsv, results/TS_games.tsv,
# results/TS_series.npz, results/summary_*.txt.
set -eu
cd "$(dirname "$0")"
if [ "${1:-}" = fetch ]; then
  rk=${2:?usage: collect.sh fetch <remote-kit path on cannon>}
  rsync -a cannon:"$rk"/out/ out/
  rsync -a cannon:"$rk"/'JOBID_*' . 2>/dev/null || true
  shift 2
fi
for a in E3:tasks_e3.txt E5:tasks_e5.txt TS:tasks_ts.txt; do
  arm=${a%%:*}; tf=${a#*:}
  [ -d "out/$arm" ] || { echo "$arm: no output yet"; continue; }
  nt=$(grep -v '^[[:space:]]*#' "$tf" | grep -vc '^[[:space:]]*$'); nc=0
  for k in $(seq 1 "$nt"); do bash task.sh "$tf" "$k" --check >/dev/null 2>&1 && nc=$((nc+1)); done
  echo "$arm: $nc of $nt tasks complete"
done
PYTHONDONTWRITEBYTECODE=1 python3 collect.py "$@"
