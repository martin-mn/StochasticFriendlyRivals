#!/bin/bash
# stage.sh -- copy the kit to Cannon (~/CL1/PRFull/cannon/evoaudit) without data, outputs, Mac binaries or the venv,
# and build it there (NOBUILD=1 bash stage.sh: copy only, keep the binaries already built on Cannon).  The venv (numpy,
# gmpy2) is made once on Cannon (README).  Never touches JOBID files, outputs, logs or the venv.
set -euo pipefail
cd "$(dirname "$0")"
REMOTE=CL1/PRFull/cannon/evoaudit
ssh -o BatchMode=yes cannon true || { echo "the Cannon session is closed: ask Martin to open ssh cannon" >&2; exit 1; }
rsync -a --exclude /data --exclude '/out*' --exclude /logs --exclude '/results*' --exclude /src/review --exclude /venv --exclude '/JOBID_*' \
      --exclude '*.dylib' --exclude '*.so' --exclude __pycache__ ./ "cannon:$REMOTE/"
if [ "${NOBUILD:-0}" = 1 ]; then ssh -o BatchMode=yes cannon "cd ~/$REMOTE && mkdir -p logs out out_probe"
else ssh -o BatchMode=yes cannon "cd ~/$REMOTE && mkdir -p logs out out_probe && bash compile.sh"; fi
