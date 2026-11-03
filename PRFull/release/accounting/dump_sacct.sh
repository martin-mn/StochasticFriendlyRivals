#!/bin/bash
# dump_sacct.sh -- the Slurm accounting of every array of CS Table 2 of the Computational Supplement (SI Table 8
# before it, and SI Table 7 before the second referee round; the file names keep that old number), from Harvard
# Cannon (read-only queries).
#
#   bash dump_sacct.sh            (needs "ssh cannon" to work non-interactively; writes only in this directory)
#
# Writes, for the job ids of table7_jobs.tsv:
#   sacct_table7.psv    sacct --parsable2, every array task and every job step (batch, extern, srun steps), with
#                       JobID, JobIDRaw, JobName, Partition, State, ExitCode, Submit, Start, End, ElapsedRaw,
#                       CPUTimeRAW, TotalCPU, MaxRSS, ReqMem, NCPUS, NodeList
#   jobs_tsv_lines.tsv  the lines of the cross-project submission ledger ~/jobs.tsv on Cannon that name these ids
#   dump_sacct.log      when, and how many lines
# check_table7.py then recomputes CS Table 2 from sacct_table7.psv.  Nothing is written on Cannon.
set -euo pipefail
cd "$(dirname "$0")"
ids=$(grep -v '^#' table7_jobs.tsv | cut -f1 | paste -sd, -)
n=$(grep -vc '^#' table7_jobs.tsv)
F=JobID,JobIDRaw,JobName%60,Partition,State,ExitCode,Submit,Start,End,ElapsedRaw,CPUTimeRAW,TotalCPU,MaxRSS,ReqMem,NCPUS,NodeList
ssh -o BatchMode=yes cannon "sacct -j $ids --format=$F --parsable2" > sacct_table7.psv.new
pat=$(echo "$ids" | tr ',' '|')
ssh -o BatchMode=yes cannon "grep -P '\t($pat)\t' ~/jobs.tsv" > jobs_tsv_lines.tsv.new || true
mv sacct_table7.psv.new sacct_table7.psv
mv jobs_tsv_lines.tsv.new jobs_tsv_lines.tsv
{ echo "dumped $(date -u +%FT%TZ) from cannon: $n arrays"
  echo "sacct_table7.psv: $(($(wc -l < sacct_table7.psv) - 1)) lines (tasks and steps)"
  echo "jobs_tsv_lines.tsv: $(wc -l < jobs_tsv_lines.tsv | tr -d ' ') ledger lines"; } | tee dump_sacct.log
