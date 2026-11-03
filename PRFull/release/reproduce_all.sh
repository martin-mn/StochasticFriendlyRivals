#!/usr/bin/env bash
# reproduce_all.sh -- redraw every figure of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity", write or
# recompute the tables that come from deposited data, rerun the exact certificates, and rebuild the paper and its Computational Supplement, all on
# working copies; compare every output with the deposited file.  Release 1.0.0 does not contain the paper's sources
# (they are added in release 1.1.0, when the article is published): then the LaTeX step is skipped, and the script says so.
#
#   bash PRFull/release/reproduce_all.sh            reproduce_figures.sh, reproduce_certificates.sh, the LaTeX build
#                                                   (about 2 h 30 min on 3 cores, nearly all of it the certificates)
#   bash PRFull/release/reproduce_all.sh --quick    both scripts with --quick, and the LaTeX builds (about 13 min, 2.5 GB)
#   options:  --out DIR     working directory (default PRFull/release/work)
#             --threads N   threads of the certificate programs (default 3)
#             --no-latex    skip the LaTeX builds (ms.tex and cs.tex; skipped anyway when the sources are absent)
#             --latex-only  only the LaTeX builds, in the working copy of an earlier run of reproduce_figures.sh
#   PYTHON=..., CC=... as for the two scripts.
#
# The LaTeX build runs in the figure script's working copy (DIR/figures/PRFull/Text), so that it includes the
# figures and SI Table 4 just redrawn; in quick mode every figure is redrawn, and a figure that was not (--only
# runs of reproduce_figures.sh) is taken from the deposit and said so.  It runs Text/reorder_bib.py, which must leave
# ms.tex unchanged (the reference lists in the deposit are the ones the pool gives), then pdflatex three times, and
# fails on a LaTeX error or an undefined reference.  It then does the same for the Computational Supplement:
# Text/cs_bib.py must leave cs.tex unchanged, and cs.tex is built by pdflatex three times.  Neither PDF is compared
# with the deposited one (pdfTeX writes the date and the included figures' dates).
#
# Exit status: 0 all passed and every rewritten file is the deposited one; 2 all passed but some file differs;
# 1 a step failed.  DIR/logs/summary.txt collects the three summaries.

set -u
REL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRF="$(dirname "$REL")"
PY="${PYTHON:-python3}"
MODE=""; OUT="$REL/work"; NTH=3; LATEX=1; ONLYTEX=0
while [ $# -gt 0 ]; do
  case "$1" in
    --quick) MODE=--quick ;;
    --out) OUT="$2"; shift ;;
    --threads) NTH="$2"; shift ;;
    --no-latex) LATEX=0 ;;
    --latex-only) ONLYTEX=1 ;;
    -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
  shift
done
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
T0=$(date +%s)
echo "reproduce_all.sh ${MODE:---full} in $(dirname "$PRF"), working directory $OUT"

SF=0; SC=0
if [ "$ONLYTEX" = 0 ]; then
  bash "$REL/reproduce_figures.sh" $MODE --out "$OUT"; SF=$?
  bash "$REL/reproduce_certificates.sh" $MODE --out "$OUT" --threads "$NTH"; SC=$?
elif [ ! -d "$OUT/figures/PRFull/Text" ]; then
  echo "--latex-only: no working copy in $OUT/figures (run reproduce_figures.sh first)"; exit 1
fi

SL=0; NOTEX=""
LOG="$OUT/logs/latex"; mkdir -p "$LOG"; rm -f "$LOG/summary.jsonl"
if [ "$LATEX" = 1 ] && { [ ! -f "$PRF/Text/ms.tex" ] || [ ! -f "$PRF/Text/cs.tex" ]; }; then
  NOTEX="LaTeX: skipped, the paper's sources are not in this release (PRFull/Text/ms.tex, cs.tex; they are added in release 1.1.0, when the article is published)"
  echo; echo "$NOTEX"
  LATEX=0
  [ "$ONLYTEX" = 1 ] && exit 1
fi
if [ "$LATEX" = 1 ]; then
  WT="$OUT/figures/PRFull/Text"; WF="$OUT/figures/PRFull/Figures"
  for f in Fig1 Fig2 Fig3 Fig4 Fig5 SIFig1 SIFig2 SIFig3; do
    [ -f "$WF/$f.pdf" ] || { cp -p "$PRF/Figures/$f.pdf" "$WF/$f.pdf"; echo "  LaTeX: $f.pdf was not redrawn here; the deposited one is used"; }
  done
  [ -f "$WT/si_table4.tex" ] || cp -p "$PRF/Text/si_table4.tex" "$WT/si_table4.tex"
  "$PY" "$REL/tools/runstep.py" --name latex --label "the paper and the Computational Supplement: reorder_bib.py and cs_bib.py, then pdflatex three times each" --cwd "$WT" \
      --logdir "$LOG" --summary "$LOG/summary.jsonl" --need pdflatex --scan "$WT=$PRF/Text" \
      --ignore '^(ms|cs)\.(pdf|aux|log|out|fls|fdb_latexmk|synctex\.gz)$' --ignore '\.csv$|^si_table4\.tex$' --nocounts \
      --expect 'Output written on ms\.pdf \([0-9]+ pages' --must-write ms.pdf \
      --expect 'Output written on cs\.pdf \([0-9]+ pages' --must-write cs.pdf -- bash -c "
      set -e; cp -p ms.tex ms.tex.before
      $PY reorder_bib.py
      cmp -s ms.tex ms.tex.before && echo 'reorder_bib.py left ms.tex unchanged' || { echo 'FAIL: reorder_bib.py changed ms.tex'; exit 1; }
      rm -f ms.tex.before
      for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error ms.tex > pdflatex_\$i.log 2>&1 || { tail -30 pdflatex_\$i.log; exit 1; }; done
      grep -E '^Output written on' ms.log | tail -1
      if grep -Eq 'There were undefined references|Citation .* undefined|Reference .* undefined' ms.log; then echo 'FAIL: undefined references'; exit 1; fi
      echo \"ms.tex: LaTeX warnings in the last pass: \$(grep -c 'LaTeX Warning' ms.log || true); overfull boxes: \$(grep -c '^Overfull' ms.log || true)\"
      cp -p cs.tex cs.tex.before
      $PY cs_bib.py
      cmp -s cs.tex cs.tex.before && echo 'cs_bib.py left cs.tex unchanged' || { echo 'FAIL: cs_bib.py changed cs.tex'; exit 1; }
      rm -f cs.tex.before
      for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error cs.tex > pdflatex_cs_\$i.log 2>&1 || { tail -30 pdflatex_cs_\$i.log; exit 1; }; done
      grep -E '^Output written on' cs.log | tail -1
      if grep -Eq 'There were undefined references|Citation .* undefined|Reference .* undefined' cs.log; then echo 'FAIL: undefined references in cs.tex'; exit 1; fi
      echo \"cs.tex: LaTeX warnings in the last pass: \$(grep -c 'LaTeX Warning' cs.log || true); overfull boxes: \$(grep -c '^Overfull' cs.log || true)\""
  SL=$?
  [ $SL -eq 3 ] && SL=0          # pdflatex not installed: skipped, said so
fi

{ echo "reproduce_all.sh ${MODE:---full}: $(( $(date +%s) - T0 )) s in all"
  if [ "$ONLYTEX" = 0 ]; then
    echo; cat "$OUT/logs/figures/summary.txt" 2>/dev/null
    echo; cat "$OUT/logs/certificates/summary.txt" 2>/dev/null
  fi
  [ -n "$NOTEX" ] && { echo; echo "$NOTEX"; }
  if [ "$LATEX" = 1 ]; then echo; "$PY" -c "
import json, sys
for l in open(sys.argv[1]):
    r = json.loads(l)
    print('latex: %s in %s s, peak memory %s MiB%s' % (r['status'], r.get('seconds', '-'), r.get('peak_rss_mib', '-'),
          ''.join('; ' + w for w in r.get('why', []))))
" "$LOG/summary.jsonl" 2>/dev/null; grep -E '^(Output written on|ms\.tex: LaTeX|cs\.tex: LaTeX)' "$LOG/latex.out" 2>/dev/null; fi
} > "$OUT/logs/summary.txt"
echo; echo "=================================================================== reproduce_all.sh summary"
cat "$OUT/logs/summary.txt"
st=0
for s in $SF $SC $SL; do
  if [ "$s" = 1 ]; then st=1; elif [ "$s" = 2 ] && [ $st = 0 ]; then st=2; fi
done
exit $st
