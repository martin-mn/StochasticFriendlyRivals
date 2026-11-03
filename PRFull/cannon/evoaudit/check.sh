#!/bin/bash
# check.sh <arm> -- completion check of one array of this kit (arm = probe or full), from its JOBID_<arm> file only:
# every task COMPLETED in sacct, every .err log empty, every chunk_<k>.done present with all its records and no
# errors.  Prints PASS or the offending items.  Never looks at the queue as a whole.
set -euo pipefail
cd "$(dirname "$0")"
arm=${1:?probe or full}; [ -f "JOBID_$arm" ] || { echo "no JOBID_$arm: this kit has no known $arm array"; exit 1; }
jid=$(cat "JOBID_$arm"); out=$([ "$arm" = probe ] && echo out_probe || echo out)
ssh -o BatchMode=yes cannon "cd ~/CL1/PRFull/cannon/evoaudit && \
  echo '--- sacct tasks not COMPLETED:'; sacct -j $jid --format=JobID,State,ExitCode,Elapsed -X --noheader | awk '\$2 != \"COMPLETED\"'; \
  echo '--- tasks:'; sacct -j $jid --format=JobID -X --noheader | wc -l; \
  echo '--- non-empty stderr:'; find logs -name '${arm}_${jid}_*.err' -size +0 ; \
  echo '--- chunks incomplete or with errors:'; for f in $out/chunk_*.done; do awk -v f=\$f '{ n = split(\$0, a, \"errors \"); e = (n > 1) ? a[2] + 0 : 0; if (\$1 != \$3 || e != 0) print f, \$0 }' \$f; done; \
  echo '--- done files:'; ls $out/chunk_*.done | wc -l; \
  echo '--- elapsed (max):'; sacct -j $jid --format=Elapsed -X --noheader | sort | tail -1"
