#!/usr/bin/env bash
# reproduce_figures.sh -- redraw the figures of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity", recompute
# every cell of SI Table 1, write SI Table 4 and SI Data 1-3, and recompute the numbers behind SI Tables 3 and 5 and CS Table 2, from the
# deposited data, on a working copy; compare every output with the deposited file.
#
#   bash PRFull/release/reproduce_figures.sh            everything below (about 6 min, one core; 3 for the receipts)
#   bash PRFull/release/reproduce_figures.sh --quick    the steps marked *: every figure, SI Tables 1, 3 and 4, CS Table 2, SI Data (about 2 min)
#   options:  --out DIR       working directory (default PRFull/release/work); logs in DIR/logs/figures
#             --only A,B,...  run only these steps (names as in the list below), in the existing working copy
#             --record        after the run, write PRFull/release/figure_inputs.json: every file under the tree that
#                             each figure script read (traced, tools/trace_reads.py); assemble_release.py uses it
#             --list          print the steps, run nothing
#   PYTHON=/path/to/python3 selects the interpreter (numpy, scipy, matplotlib; gmpy2 for fig3's exact checks).
#
# Nothing in the deposit is changed.  PRFull/Figures (the scripts only, no outputs) and PRFull/Text (the sources,
# no build products) are copied to DIR/figures/PRFull; PRFull/cannon, PRFull/check and PRFull/notes are linked there
# read-only, and the data of PartnersRivalsM1f, PartnersRivalsM2f, PartnersRivals1, PartnersRivals and DiskM1WF are
# read in place.  Every script runs in the copy under nice -n 10, traced, and every file it writes is compared with
# the deposited file at the same place (tools/compare.py: PDFs apart from their creation date, PNGs by pixels, text
# apart from run times, dates and paths).  A step fails if it exits non-zero, prints a failure, or does not write
# the files it must write.  The kit analyses (SI Tables 3 and 5) run in copies of their kits with the raw output
# linked, and compare their tables with the deposited ones; CS Table 2 (the cost, in the Computational Supplement
# cs.tex; it was SI Table 8, and SI Table 7 before the second referee round, whence the file names check_table7.py,
# sacct_table7.psv) is recomputed from the deposited Slurm accounting.
#
# Steps (quick mode runs those marked *):
#   * fig1         Figures/fig1.py        Figure 1 (and fig1data.py, its numbers with their provenance)
#   * fig2         Figures/fig2.py        Figure 2
#   * fig3         Figures/fig3.py        Figure 3
#   * fig45        Figures/fig45.py       Figures 4 and 5, SI Figures 1 and 2, fig45_stats.txt
#   * sitable1     Figures/sitable1.py    SI Table 1 (sitable1.out) from the census and the dimension probe
#   * sitable4     Figures/sitable4.py    SI Table 4 (Text/si_table4.tex) from fig45_stats.txt (--only: the deposited one)
#   * sifig3       Figures/sifig3.py      SI Figure 3, sifig3_stats.txt
#   * sidata1      Figures/mk_sidata1.py  SI Data 1 (Text/SIData1_games512.csv), then sidata1_check.py
#     receipts     fig2_check.py, fig2_check_f.py, fig3_witness_check.py (checks of Figures 2 and 3 against raw tests)
#   * table3       cannon/sweep13/collect.py  SI Table 3 (sweep13.json) from the 512 task outputs at alpha = 1/3
#     table5       cannon/cw2/collect.py and cannon/wfc2/analyze.py  SI Table 5 (the per-game tables)
#   * sidata23     release/tools/mk_sidata23.py  SI Data 2 and 3 from the kits' per-game tables (after table5: the
#                  tables just rewritten; otherwise the deposited ones)
#   * table7       release/accounting/check_table7.py  CS Table 2 (the cost) from the Slurm accounting
#
# Exit status: 0 every step passed and every rewritten file is the deposited one; 2 every step passed but some
# file differs; 1 a step failed.  The summary is DIR/logs/figures/summary.txt.

set -u
REL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRF="$(dirname "$REL")"
ROOT="$(dirname "$PRF")"
PY="${PYTHON:-python3}"
MODE=full; OUT="$REL/work"; ONLY=""; RECORD=0; LIST=0
while [ $# -gt 0 ]; do
  case "$1" in
    --quick) MODE=quick ;;
    --out) OUT="$2"; shift ;;
    --only) ONLY=",$2,"; shift ;;
    --record) RECORD=1 ;;
    --list) LIST=1 ;;
    -h|--help) sed -n '2,43p' "$0"; exit 0 ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
  shift
done
if [ "$LIST" = 1 ]; then grep -E '^#   [ *] [a-z][A-Za-z0-9_-]* ' "$0" | sed 's/^#   //'; exit 0; fi
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
W="$OUT/figures"; WP="$W/PRFull"; LOG="$OUT/logs/figures"; SUM="$LOG/summary.jsonl"
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg PYTHONPYCACHEPREFIX="$W/.pycache"   # read every module from its source (traced), write no bytecode

command -v "$PY" >/dev/null 2>&1 || { echo "no Python interpreter '$PY'"; exit 1; }
"$PY" -c "import numpy, scipy, matplotlib" 2>/dev/null || { echo "$PY cannot import numpy, scipy and matplotlib (see PRFull/release/RELEASE.md)"; exit 1; }
T0=$(date +%s)
echo "reproduce_figures.sh ($MODE${ONLY:+, only ${ONLY//,/ }}) in $ROOT"
echo "  $("$PY" -c 'import sys, numpy, scipy, matplotlib; print("Python " + sys.version.split()[0] + ", numpy " + numpy.__version__ + ", scipy " + scipy.__version__ + ", matplotlib " + matplotlib.__version__)')"
echo "  working copy $WP, logs $LOG"

"$PY" "$REL/relocate.py" | tail -1
if [ -z "$ONLY" ] || [ ! -d "$WP/Figures" ]; then       # a fresh working copy (--only reuses the last one)
  rm -rf "$W" "$LOG"; mkdir -p "$WP/Figures" "$WP/Text" "$LOG"
  rsync -a --include '*.py' --exclude '*' "$PRF/Figures/" "$WP/Figures/"
  cp -p "$PRF/Figures/fig1_reference_values.json" "$WP/Figures/"   # an input of fig1data.py (Table 1, deposited)
  [ -f "$PRF/Figures/sitable4_legend.tex" ] && cp -p "$PRF/Figures/sitable4_legend.tex" "$WP/Figures/"   # the legend of SI Table 4 (part of the paper: release 1.1.0)
  rsync -a --exclude 'ms.aux' --exclude 'ms.log' --exclude 'ms.out' --exclude 'ms.pdf' --exclude 'ms.fls' \
        --exclude 'ms.fdb_latexmk' --exclude 'ms.synctex.gz' --exclude 'cs.aux' --exclude 'cs.log' --exclude 'cs.out' \
        --exclude 'cs.pdf' --exclude 'cs.fls' --exclude 'cs.fdb_latexmk' --exclude 'cs.synctex.gz' --exclude '*.bak' \
        "$PRF/Text/" "$WP/Text/"
  for d in cannon check notes; do ln -s "$PRF/$d" "$WP/$d"; done
  "$PY" "$REL/relocate.py" --map "$PRF/Figures=$WP/Figures" --map "$PRF/Text=$WP/Text" "$WP/Figures" | sed 's/^/  /'
fi
mkdir -p "$LOG"

NFAIL=0
want() {   # is step $1 selected?
  local s=$1 q=$2
  if [ -n "$ONLY" ]; then case "$ONLY" in *",$s,"*) return 0 ;; *) return 1 ;; esac; fi
  [ "$MODE" = full ] || [ "$q" = "*" ]
}
# the deposited outputs were written under the author's tree; paths under it in them are compared as paths, so that the
# release verifies wherever it is unpacked
AUTHOR_ROOT=/Users/martin/Documents/CL1
step() {   # step NAME LABEL DIR [runstep options] -- COMMAND ...
  local name=$1 label=$2 dir=$3; shift 3
  "$PY" "$REL/tools/runstep.py" --name "$name" --label "$label" --cwd "$dir" --logdir "$LOG" --summary "$SUM" \
      --strip "$WP" --strip "$PRF" --strip "$ROOT" --strip "$AUTHOR_ROOT/PRFull" --strip "$AUTHOR_ROOT" "$@"
  [ $? -eq 1 ] && NFAIL=$((NFAIL + 1))
  return 0
}
F="$WP/Figures"
SC=(--scan "$WP=$PRF")
traced() { # traced NAME SCRIPT ARGS...: the command that runs a figure script traced
  local n=$1; shift
  echo "$PY" "$REL/tools/trace_reads.py" --root "$ROOT" --log "$LOG/trace_$n.json" --alias "$WP=$PRF" --cwd "$F" "$@"
}

if want fig1 "*"; then
  step fig1 "Figure 1 (fig1.py; fig1data.py prints its numbers with their provenance)" "$F" "${SC[@]}" \
       --must-write Fig1.pdf --must-write Fig1.png -- bash -c "$(traced fig1data fig1data.py) > fig1data.out && $(traced fig1 fig1.py)"
fi
if want fig2 "*"; then
  step fig2 "Figure 2 (fig2.py)" "$F" "${SC[@]}" --must-write Fig2.pdf --must-write Fig2.png -- $(traced fig2 fig2.py)
fi
if want fig3 "*"; then
  step fig3 "Figure 3 (fig3.py)" "$F" "${SC[@]}" --must-write Fig3.pdf --must-write Fig3.png -- $(traced fig3 fig3.py)
fi
# the inputs of fig45 and sitable4 that the GitHub repository does not carry: they are in the Zenodo record's archive
# PRFull-1.0.0_binary-paper-dumps.tar.gz (fetch_zenodo.sh); without them the two steps are skipped and say so
DUMPS="$ROOT/DiskM2WF/Data/dw1_e4.bin $ROOT/DiskM2WF/Data/dw2_e4.bin"
DUMPS1="$ROOT/PartnersRivals1/Data/f1_e4.bin $ROOT/PartnersRivals1/Data/f2_e4.bin $ROOT/PartnersRivals1/Cannon/atoms512.bin"
NEEDD=(--need-what "they are in the archive PRFull-1.0.0_binary-paper-dumps.tar.gz of the Zenodo record; bash PRFull/release/fetch_zenodo.sh fetches it")
for f in $DUMPS; do NEED45+=(--need-file "$f"); done
NEEDT4=("${NEED45[@]}"); for f in $DUMPS1; do NEEDT4+=(--need-file "$f"); done
if want fig45 "*"; then
  step fig45 "Figures 4 and 5, SI Figures 1 and 2 (fig45.py)" "$F" "${SC[@]}" "${NEED45[@]}" "${NEEDD[@]}" --must-write Fig4.pdf --must-write Fig5.pdf \
       --must-write SIFig1.pdf --must-write SIFig2.pdf --must-write fig45_stats.txt -- $(traced fig45 fig45.py)
fi
# the raw output of the dimension probe, which sitable1 reads, is in the GitHub repository and the Zenodo record's code
# archive (hosting/hosting_plan.py; in the author's tree it is where assemble_release.py harvested it); without it the
# step is skipped and says so
CA1="$ROOT/PartnersRivalsM2f/cannon/ca1_out"
[ -d "$CA1" ] || [ ! -d "$PRF/release/harvest/ca1_out" ] || CA1="$PRF/release/harvest/ca1_out"
NEEDT1=(--need-what "they are in the GitHub repository and in the archive PRFull-1.0.0_code.tar.gz of the Zenodo record")
i=0
for cs in C_1 C_-1 A_1 A_-1 Z_1 Z_-1; do                       # dim.sbatch: task i = 64 * block + stripe
  for s in $(seq 0 63); do NEEDT1+=(--need-file "$CA1/dim_${cs}_$s.txt" --need-file "$CA1/dim_45752239_$i.err"); i=$((i + 1)); done
done
if want sitable1 "*"; then
  step sitable1 "SI Table 1 (sitable1.py: the census and the dimension probe)" "$F" "${SC[@]}" "${NEEDT1[@]}" \
       --must-write sitable1.out --expect '^checks failed: 0 of [0-9]+$' -- $(traced sitable1 sitable1.py)
fi
if want sitable4 "*"; then
  [ -f "$F/fig45_stats.txt" ] || cp -p "$PRF/Figures/fig45_stats.txt" "$F/fig45_stats.txt"   # --only: the deposited statistics
  step sitable4 "SI Table 4 (sitable4.py -> Text/si_table4.tex)" "$F" "${SC[@]}" "${NEEDT4[@]}" "${NEEDD[@]}" --must-write ../Text/si_table4.tex \
       --must-write sitable4.out -- $(traced sitable4 sitable4.py)
fi
if want sifig3 "*"; then
  step sifig3 "SI Figure 3 (sifig3.py)" "$F" "${SC[@]}" --must-write SIFig3.pdf --must-write sifig3_stats.txt -- $(traced sifig3 sifig3.py)
fi
if want sidata1 "*"; then
  step sidata1 "SI Data 1 (mk_sidata1.py -> Text/SIData1_games512.csv; sidata1_check.py)" "$F" "${SC[@]}" \
       --must-write ../Text/SIData1_games512.csv -- \
       bash -c "$(traced sidata1 mk_sidata1.py ../Text/SIData1_games512.csv) && $PY sidata1_check.py > sidata1_check.out"
fi
if want receipts ""; then
  # fig3_witness_check.py runs 6 processes; the copy runs 3 (the local cap).  Its arg-min co-player is one of 27 840
  # exact ties and depends on the order the processes finish, so it is not compared.
  sed -i.orig 's/^NPROC = 6$/NPROC = 3/' "$F/fig3_witness_check.py" && rm -f "$F/fig3_witness_check.py.orig"
  step receipts "the receipts of Figures 2 and 3 (fig2_check.py, fig2_check_f.py, fig3_witness_check.py)" "$F" "${SC[@]}" \
       --normalize 'at q = +[0-9]+ +(?=\[[0-9]+ co-players at 0\])' \
       --must-write cache/fig2_check.txt --must-write cache/fig2_check_f.txt --must-write fig3_witness_check.out -- \
       bash -c "$(traced fig2_check fig2_check.py) && $(traced fig2_check_f fig2_check_f.py) && $PY fig3_witness_check.py > fig3_witness_check.out"
fi

kitcopy() {   # kitcopy KIT BIGDIR...: a copy of PRFull/cannon/KIT in $W/kits/KIT, its raw output directories linked
  local k=$1; shift
  local ex=() d
  for d in "$@"; do ex+=(--exclude "/$d"); done
  rm -rf "$W/kits/$k"; mkdir -p "$W/kits/$k"                   # rm -r removes the links, never what they point to
  rsync -a --exclude __pycache__ --exclude .DS_Store "${ex[@]}" "$PRF/cannon/$k/" "$W/kits/$k/"
  for d in "$@"; do [ -e "$PRF/cannon/$k/$d" ] && ln -sn "$PRF/cannon/$k/$d" "$W/kits/$k/$d"; done
  "$PY" "$REL/relocate.py" --map "$PRF/cannon/$k=$W/kits/$k" "$W/kits/$k" > /dev/null
}
if want table3 "*"; then
  kitcopy sweep13 out_13 out_probe orig ref validate
  K="$W/kits/sweep13"
  step table3 "SI Table 3 (cannon/sweep13/collect.py on the 512 task outputs at alpha = 1/3)" "$K" \
       --scan "$K=$PRF/cannon/sweep13" --must-write sweep13.json -- \
       "$PY" collect.py out_13 --alpha 1/3 --seedbase 3000000 \
       --compare "$ROOT/PartnersRivalsM2f/check/sweep_a333.tsv" --json sweep13.json
fi
if want table5 ""; then
  # the raw output of the sensitivity runs is in the Zenodo record's archive PRFull-1.0.0_evolution-runs.tar.gz; the
  # GitHub repository carries only the 80 run-length series of SI Figure 3.  One file of each arm the step reads:
  NEED5=(--need-what "they are in the archive PRFull-1.0.0_evolution-runs.tar.gz of the Zenodo record; bash PRFull/release/fetch_zenodo.sh evolution-runs fetches it")
  for f in cw2/out/E3/1013 cw2/out/E5/1 cw2/out/TS/p017_r01/t54 wfc2/out/TS/r_1o3_110_1.txt \
           wfc2/out/{A1/r_1o10,A2/r_1o5,A3/r_2o5,K0/r_1o3,E3/r_1o3,E5/r_1o3,E3F2/r_1o3,E5F2/r_1o3}_1_1.txt \
           wfc2/out/{B0,B0b,B0K0,B0K0b}/r_1o3_30_1.txt; do NEED5+=(--need-file "$PRF/cannon/$f"); done
  kitcopy cw2 out ref validation build results_local
  kitcopy wfc2 out ref
  step table5 "SI Table 5 (cannon/cw2/collect.py, memory one; cannon/wfc2/analyze.py and b0pool.py, memory two)" "$W/kits" \
       "${NEED5[@]}" --scan "$W/kits/cw2=$PRF/cannon/cw2" --scan "$W/kits/wfc2=$PRF/cannon/wfc2" --ignore '^analysis/summary\.txt$' \
       --must-write cw2/results/E3_games.tsv --must-write cw2/results/TS_games.tsv --must-write wfc2/analysis/K0_games.tsv \
       --must-write wfc2/analysis/b0pool.txt -- bash -c "
       set -e; cd cw2; $PY collect.py > collect.out; cd ../wfc2
       $PY analyze.py A1 A2 A3 K0 E3 E5 E3F2 E5F2 TS > analyze.out; $PY b0pool.py > b0pool.out"
fi
if want sidata23 "*"; then
  K5="$W/kits"; SD1="$WP/Text/SIData1_games512.csv"; [ -f "$SD1" ] || SD1="$PRF/Text/SIData1_games512.csv"
  if [ -f "$K5/wfc2/analyze.out" ] && [ -f "$K5/cw2/results/TS_games.tsv" ]; then
    SRC=(--wfc2 "$K5/wfc2" --cw2 "$K5/cw2" --wfc2-summary "$K5/wfc2/analyze.out"); WHAT="the tables table5 just rewrote"
  else
    SRC=(--wfc2 "$PRF/cannon/wfc2" --cw2 "$PRF/cannon/cw2"); WHAT="the deposited per-game tables"
  fi
  step sidata23 "SI Data 2 and 3 (release/tools/mk_sidata23.py, from $WHAT)" "$WP/Text" --scan "$WP/Text=$PRF/Text" \
       --must-write SIData2_sensitivity.csv --must-write SIData3_runlength.csv \
       --expect 'SIData2_sensitivity.csv: 616 rows, byte-identical' --expect 'SIData3_runlength.csv: 48 rows, byte-identical' -- \
       "$PY" "$REL/tools/mk_sidata23.py" "${SRC[@]}" --sidata1 "$SD1" --out .
fi
if want table7 "*"; then
  step table7 "CS Table 2 (release/accounting/check_table7.py on the deposited sacct dump)" "$REL/accounting" \
       --no-fail-words --nocounts --expect 'rows or arrays that differ: 0' -- "$PY" check_table7.py
fi

if [ "$RECORD" = 1 ]; then
  "$PY" - "$LOG" "$ROOT" "$REL/figure_inputs.json" <<'EOF'
import glob, json, os, sys
log, root, dst = sys.argv[1:]
old = json.load(open(dst)) if os.path.exists(dst) else {"steps": {}}
for f in sorted(glob.glob(os.path.join(log, "trace_*.json"))):
    t = json.load(open(f))
    n = os.path.basename(f)[6:-5]
    old["steps"][n] = {"script": os.path.relpath(t["script"], root) if t["script"].startswith(root) else t["script"],
                       "exit": t["exit"], "seconds": t["seconds"],
                       "reads": [os.path.relpath(p, root) for p in t["reads"]],
                       "via_links": {os.path.relpath(k, root): os.path.relpath(v, root) for k, v in t.get("realpath", {}).items()}}
old["note"] = ("files under the tree root read by each figure script, traced by tools/trace_reads.py during "
               "reproduce_figures.sh --record; paths relative to the directory that holds PRFull/")
json.dump(old, open(dst, "w"), indent=1, sort_keys=True)
print("recorded the inputs of %d traced scripts in %s" % (len(old["steps"]), dst))
EOF
fi

"$PY" - "$SUM" "$LOG/summary.txt" "$MODE" "$(( $(date +%s) - T0 ))" <<'EOF'
import json, os, sys
recs = [json.loads(l) for l in open(sys.argv[1])] if os.path.exists(sys.argv[1]) else []
last = {r["name"]: i for i, r in enumerate(recs)}                 # --only reruns: the latest record of each step
recs = [r for i, r in enumerate(recs) if last[r["name"]] == i]
out = ["reproduce_figures.sh --%s: %d steps in %s s" % (sys.argv[3], len(recs), sys.argv[4])]
out.append("%-10s %-7s %8s %9s %7s %6s %7s %9s" % ("step", "status", "seconds", "peak MiB", "written", "same", "differ", "unchecked"))
for r in recs:
    out.append("%-10s %-7s %8s %9s %7s %6s %7s %9s" % (r["name"], r["status"], r.get("seconds", "-"), r.get("peak_rss_mib", "-"),
               r.get("written", "-"), r.get("same", "-"), r.get("differ", "-"), r.get("new", "-")))
nf = sum(r["status"] == "failed" for r in recs); nd = sum(r.get("differ", 0) for r in recs)
ns = sum(r["status"] == "skipped" for r in recs)
out.append("steps failed: %d; skipped: %d; files that differ from the deposit: %d; largest peak memory %.0f MiB"
           % (nf, ns, nd, max([r.get("peak_rss_mib", 0) for r in recs] or [0])))
for r in recs:
    if r["status"] == "skipped":
        out.append("  skipped: %s: %s" % (r["name"], r.get("why", "")))
        continue
    for w in r.get("why", []):
        out.append("  failed: %s: %s" % (r["name"], w))
    for d in r.get("differing", []):
        out.append("  differs: %s: %s" % (r["name"], d))
open(sys.argv[2], "w").write("\n".join(out) + "\n"); print("\n" + "\n".join(out))
sys.exit(1 if nf else (2 if nd else 0))
EOF
