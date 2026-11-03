#!/bin/bash
# collect.sh <JOBIDfile> -- aggregate one alpha array.  REFUSES to report unless every task is
# COMPLETED per sacct (an empty queue is not evidence; see the project job protocol).
set -euo pipefail
jid=$(cat "$1")
bad=$(ssh -o BatchMode=yes cannon "sacct -j $jid --format=JobID,State,ExitCode -X --noheader | awk '\$2 != \"COMPLETED\"'")
if [ -n "$bad" ]; then echo "NOT CLEAN -- refusing to aggregate:"; echo "$bad"; exit 1; fi
ssh -o BatchMode=yes cannon "cd ~/PartnersRivalsM2f_alpha && cat out/a_*.tsv" | awk -F'\t' '
  $1=="ALPHA" { a=$2; n[a]+=$3; effC[a]+=$4; effA[a]+=$5; effZ[a]+=$6; rP[a]+=$7; rM[a]+=$8;
                W[a]+=$9; S[a]+=$10; N[a]+=$11; E[a]+=$12; yP[a]+=$13; yM[a]+=$14; t[a]++ }
  END { printf "%-7s %6s %14s %10s %10s %10s %10s %12s %12s %12s %12s\n",
               "alpha","tasks","n","effCC","effALT","rivP","rivM","FR_W","FR_S","FR_N","FR_E"
        for (a in n) printf "%-7.3f %6d %14d %10d %10d %10d %10d %12d %12d %12d %12d\n",
               a/1000, t[a], n[a], effC[a], effA[a], rP[a], rM[a], W[a], S[a], N[a], E[a] }' | sort -k1,1n
