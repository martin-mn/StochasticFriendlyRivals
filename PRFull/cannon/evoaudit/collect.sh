#!/bin/bash
# collect.sh <outdir> -- copy <outdir> (out_probe or out) and the logs back from Cannon.
set -euo pipefail
cd "$(dirname "$0")"
REMOTE=CL1/PRFull/cannon/evoaudit
rsync -a "cannon:$REMOTE/${1:?outdir}/" "$1/"
rsync -a "cannon:$REMOTE/logs/" logs/
