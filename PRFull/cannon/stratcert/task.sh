#!/bin/bash
# One array task of the stratcert kit:  ./task.sh TASK NTASKS OUTDIR
# Runs stratcert.py over the units with uid % NTASKS == TASK on $SLURM_CPUS_PER_TASK worker processes
# (each single-threaded: CERT_MAX_THREADS=1).  Idempotent: units already in OUTDIR/t_TASK.tsv are skipped.
set -e
cd "$(dirname "$(readlink -f "$0")")"
TASK=$1; NTASKS=$2; OUT=$3
NCPU=${SLURM_CPUS_PER_TASK:-1}
export CERT_MAX_THREADS=1 OMP_NUM_THREADS=1
if [ "$(uname)" != "Darwin" ]; then
  PY=${PY:-/n/sw/Mambaforge-22.11.1-4/bin/python3}; export PYTHONPATH=$PWD/pylib
else
  PY=${PY:-python3}
fi
mkdir -p "$OUT"
echo "host $(hostname) task $TASK/$NTASKS ncpu $NCPU start $(date +%FT%T)"
if [ -x /usr/bin/time ] && [ "$(uname)" != "Darwin" ]; then
  /usr/bin/time -v -o "$OUT/time_$TASK.txt" $PY stratcert.py "$TASK" "$NTASKS" "$OUT" --nproc "$NCPU"
else
  $PY stratcert.py "$TASK" "$NTASKS" "$OUT" --nproc "$NCPU"
fi
echo "end $(date +%FT%T)"
