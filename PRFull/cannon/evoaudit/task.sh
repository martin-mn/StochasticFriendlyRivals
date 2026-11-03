#!/bin/bash
# task.sh <chunk> <nchunks> <outdir> -- audit one chunk of sample/sample.npy (the rows i with i % nchunks == chunk) into
# <outdir>/chunk_<chunk>.jsonl.  Idempotent: records already written are skipped, so a rerun only completes a chunk.
# One core: CERT_MAX_THREADS=1.  A successful task writes nothing to stderr.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
export CERT_MAX_THREADS=${CERT_MAX_THREADS:-1} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$3"
PY=./venv/bin/python; [ -x "$PY" ] || PY=python3
"$PY" src/audit.py run sample/sample.npy "$1" "$2" "$3"
echo "done chunk $1 of $2 -> $3"
