#!/bin/bash
# stage.sh -- copy the kit to Cannon (~/$REMOTE, default ~/PRFull_wfc2) and build it there (compile.sh checks the
# token identity with the production source first).  Needs the open ControlMaster session (ssh cannon).
# Never overwrites the remote JOBID files or outputs.
set -euo pipefail
cd "$(dirname "$0")"
REMOTE=${REMOTE:-CL1/PRFull/cannon/wfc2}
ssh -o BatchMode=yes cannon true || { echo "the Cannon session is closed: ask Martin to open ssh cannon" >&2; exit 1; }
rsync -a --exclude /out --exclude /analysis --exclude '/JOBID_*' --exclude '/wfc' --exclude '/wfc2' \
      --exclude '/wfc_*' --exclude '/wfc2_*' --exclude __pycache__ ./ "cannon:$REMOTE/"
ssh -o BatchMode=yes cannon "cd ~/$REMOTE && mkdir -p out && bash compile.sh"
