#!/bin/bash
# Pull the completed Cannon census and run the analysis.  Only the efficient-class lines are
# transferred (class C, A or Z); the '-' class is ~95% of the volume and is not needed for the
# friendly-rival tables.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
echo "== checking every task completed (sacct, per task, not squeue) =="
ssh -o BatchMode=yes cannon "sacct -j 45751662 --format=JobID,State,ExitCode -X --noheader | awk '\$2 != \"COMPLETED\" {print}'" > /tmp/census_bad.txt
if [ -s /tmp/census_bad.txt ]; then echo "NOT CLEAN:"; cat /tmp/census_bad.txt; exit 1; fi
echo "  all 512 tasks COMPLETED, none otherwise"
echo "== transferring the efficient-class lines =="
ssh -o BatchMode=yes cannon "cd ~/CL1/PartnersRivalsM2f/ca1/out && cat st_*.txt | awk '\$2!=\"-\"'" > "$HERE/st_eff.txt"
wc -l < "$HERE/st_eff.txt"
echo "== census analysis =="
python3 "$HERE/census.py" "$HERE/st_eff.txt" | tee "$HERE/census.out"
echo "== predicate-free sweep totals =="
ssh -o BatchMode=yes cannon "cd ~/CL1/PartnersRivalsM2f/ca1/out && grep -h 'sweep n=' sw_*.err | sed 's/.*: //' | awk '
{ for(i=1;i<NF;i++){ k=\$i; v=\$(i+1); if(v ~ /^[0-9]+\$/) t[k]+=v } n++ }
END{ printf \"tasks %d  (%.4g strategies)\n\", n, n*4e6
     for(k in t) printf \"  %-8s %d\n\", k, t[k] }'" | tee "$HERE/sweep.out"
