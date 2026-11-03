#!/bin/sh
# Prover 2: symbolic proof, exact numerical cross-check, root search.  Total ~1 min on 12 CPUs.
set -e
cd "$(dirname "$0")"
python3 p2_proof.py
python3 p2_numcheck.py
python3 p2_rootsearch.py
python3 p2_tables.py
