#!/bin/bash
# stage.sh -- copy the kit to Cannon (~/$REMOTE) and build it there.  Needs the open ControlMaster session
# (ssh cannon).  Never copies or overwrites outputs, JOBID files or binaries.
set -euo pipefail
cd "$(dirname "$0")"
REMOTE=${REMOTE:-CL1/PRFull/cannon/sweep13}
ssh -o BatchMode=yes cannon true || { echo "the Cannon session is closed: ask Martin to open ssh cannon" >&2; exit 1; }
ssh -o BatchMode=yes cannon "mkdir -p ~/$REMOTE"
rsync -a --exclude /out --exclude /out_probe --exclude '/JOBID_*' --exclude /census --exclude /BUILD.txt \
      --exclude __pycache__ --exclude /validate/build ./ "cannon:$REMOTE/"
ssh -o BatchMode=yes cannon "cd ~/$REMOTE && bash compile.sh"
