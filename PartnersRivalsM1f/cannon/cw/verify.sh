#!/bin/bash
# Completion of a cw array, established PER TASK from sacct -- never from squeue -u, whose output spans
# every project under this user; an empty squeue means "not running", nothing more.
#   bash verify.sh {probe|prod|kappa|probe2|prod2}          (probe2/prod2: run F2, N = 100, itend 1e8)
# Reports: every task not COMPLETED (empty = clean), non-empty stderr logs, core-hours and wall time,
# and the products of THIS array's own index set (IDX_<tag>), complete = 10 replicate lines + an h-file
# with '# total' and 105 count lines (kappa: 512 files of one 98-field line).  For the probe it also
# prints the measured rate and the production -t it implies (probe itend = production / 10).
set -e
cd "$(dirname "$0")"
TAG=${1:?usage: verify.sh probe|prod|kappa|probe2|prod2}
J="JOBID_$TAG"; [ -f "$J" ] || { echo "$J not found -- this kit has no record of that array"; exit 1; }
jid=$(cat "$J")
case "$TAG" in probe) RD=probe; LG=probe ;; prod) RD=run1; LG=prod ;; kappa) RD=kappa; LG=kappa ;;
               probe2) RD=probe2; LG=f2probe ;; prod2) RD=run2; LG=f2prod ;;
               *) echo "unknown tag $TAG"; exit 1 ;; esac
echo "array $jid ($TAG, $RD/)"
echo "--- tasks that are not COMPLETED (empty is clean) ---"
sacct -j "$jid" --format=JobID,State,ExitCode,Elapsed -X --noheader | awk '$2 != "COMPLETED"'
echo "--- non-empty stderr logs ---"
for f in log/e_${LG}*.txt; do [ -s "$f" ] && echo "$f: $(wc -l < "$f") lines, first: $(head -1 "$f")"; done 2>/dev/null || true
echo "--- core-hours, and per-task wall time ---"
sacct -j "$jid" --format=Elapsed,AllocCPUS,State -X --noheader | awk \
 '$3=="COMPLETED"{split($1,a,"[:-]"); if (length(a)==4) s=((a[1]*24+a[2])*60+a[3])*60+a[4]; else s=(a[1]*60+a[2])*60+a[3];
   t+=s*$2; n++; if (s>mx) mx=s; if (mn=="" || s<mn) mn=s}
  END{if(n) printf "%.1f core-h over %d COMPLETED tasks; wall %.1f-%.1f min, mean %.1f\n", t/3600, n, mn/60, mx/60, t/60/n;
      else print "no COMPLETED accounting rows yet"}'
if [ "$TAG" = probe ] || [ "$TAG" = probe2 ]; then
  # F1 probe: 10 x 1e6 = 1e7 generations per task (s/1000 per 1e4); F2 probe2: 10 x 1e7 = 1e8 (s/10000 per 1e4).
  # Either way the probe is 1/10 of its production task, so production = 10 x probe.
  if [ "$TAG" = probe ]; then GDIV=1000; PS=prod.slurm
    echo "--- the rate (probe: 10 replicates x itend 1e6 = 1e7 generations per task) ---"
  else GDIV=10000; PS=prod2.slurm
    echo "--- the rate (probe2, run F2: 10 replicates x itend 1e7 = 1e8 generations per task; F2 task = 1e9) ---"
  fi
  sacct -j "$jid" --format=JobID,Elapsed,State,MaxRSS -X --noheader | awk -v gdiv=$GDIV -v ps=$PS \
   '$3=="COMPLETED"{split($2,a,"[:-]"); if (length(a)==4) s=((a[1]*24+a[2])*60+a[3])*60+a[4]; else s=(a[1]*60+a[2])*60+a[3];
     printf "  %-18s %6d s   %.3f s per 1e4 generations   -> production task %.0f min\n", $1, s, s/gdiv, s*10/60;
     if (s>mx) mx=s}
    END{if (mx) printf "  worst task x 10 = %.0f min; set %s -t to about %d (3x margin)\n", mx*10/60, ps, int(mx*10*3/60+1)}'
  sacct -j "$jid" --format=JobID,MaxRSS,State --noheader | awk '$2!=""{print "  maxRSS", $1, $2}' | head -10
fi
echo "--- products of THIS array ---"
if [ "$TAG" = kappa ]; then
  nk=0; for isl in $(seq 3 3 1536); do
    [ -f "kappa/k$isl" ] && [ "$(awk '{print NF}' "kappa/k$isl")" -eq 98 ] && nk=$((nk+1)); done
  echo "$nk of 512 kappa cells complete (one 98-field line)"; exit 0
fi
I="IDX_$TAG"; [ -f "$I" ] || { echo "$I not found -- cannot scope the count to this array"; exit 0; }
ndone=0; nshort=0; nnoh=0; nmiss=0; miss=""
while read -r isl; do
  [ -n "$isl" ] || continue
  if [ ! -f "$RD/$isl" ]; then nmiss=$((nmiss+1)); miss="$miss $isl"; continue; fi
  L=$(wc -l < "$RD/$isl"); nc=0; tot=0
  if [ -f "$RD/h$isl" ]; then
    grep -q '^# total' "$RD/h$isl" && tot=1
    nc=$(grep -vc '^#' "$RD/h$isl") || nc=0
  fi
  if [ "$L" -eq 10 ] && [ "$tot" -eq 1 ] && [ "$nc" -eq 105 ]; then ndone=$((ndone+1))
  elif [ "$L" -lt 10 ]; then nshort=$((nshort+1)); miss="$miss $isl(${L}rep)"
  else nnoh=$((nnoh+1)); miss="$miss $isl(h $nc/105)"; fi
done < "$I"
echo "$ndone complete of $(wc -l < "$I" | tr -d ' ') in this array"
[ "$nmiss" -gt 0 ] && echo "  $nmiss not started"
[ "$nshort" -gt 0 ] && echo "  $nshort short on replicates"
[ "$nnoh" -gt 0 ] && echo "  $nnoh have 10 replicates but no complete h-file -- RESUBMIT THESE"
[ -n "$miss" ] && echo " incomplete:$miss"
exit 0
