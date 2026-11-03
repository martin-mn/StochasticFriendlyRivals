#!/usr/bin/env bash
# reproduce_certificates.sh -- rerun the exact certificates of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity"
# (SI section "Certificates", SI Box 2, SI Table 2; their implementation and checks, and the witness archive, are in
# the Computational Supplement, Sections 3 and 5) on a working copy, and compare everything they write with the
# deposited files.  It also reruns the exact checks of check/ref3 (SI Box 3 and the reduction to Eq. (3) on the
# boundary strata, the other memory-one stability statements, the family (0,p_CD,p_DC,1), the offer range of the
# Eq. (4) example) and the replay of the efficiency test of the rarity sweep (cannon/sweep13/validate/rcheck.c).
#
#   bash PRFull/release/reproduce_certificates.sh            everything that runs on a laptop (about 2 h 20 min on 3 cores)
#   bash PRFull/release/reproduce_certificates.sh --quick    a representative subset (about 11 min on 3 cores, 2.5 GB)
#   options:  --out DIR        working directory (default PRFull/release/work); the logs go to DIR/logs/certificates
#             --threads N      threads and processes of the certificate library (default 3; CERT_MAX_THREADS)
#             --no-bingap      skip the census of binary pairs (full mode runs all 64 tasks, about 30 min on 3 cores)
#             --only A,B,...   run only these steps of the chosen mode (the names --list prints, e.g. box_W,
#                              family_cert_pt, m1_box, rcheck), in the existing working copy (a fresh one if there is
#                              none; the C code is built only by the step build); a name that is no step of the mode
#                              is reported and makes the exit status 1; the summary counts this run's steps only and
#                              lists the earlier records of the working copy apart
#             --list           print the steps, run nothing
#   PYTHON=/path/to/python3 and CC=... select the interpreter and the C compiler (also inside check/ubound/run_all.sh
#   and check/wsep/run_all.sh); RCHECK_TASKS='0 1 2' replays only these tasks of the rarity sweep in the step rcheck.
#
# Nothing in the deposit is changed, apart from the absolute paths of its programs, which relocate.py (run first; a
# no-op in the author's tree) points at the tree and relocate.py --undo restores.  check/, cannon/bingap,
# cannon/stratcert and cannon/evoaudit (the last two without their logs), cannon/sweep13/validate and the ALPHA lines of
# cannon/sweep13/out_13 are copied to DIR/certificates/PRFull (with the paper, also Text/si_box3_legend.tex), and the
# sources of libm2f with its wrapper m2f.py to DIR/certificates/PartnersRivalsM2f; the absolute paths of the copy's
# programs are pointed at the copy (relocate.py --map), the C code is built there
# with the flags the outward rounding needs (-O2 -ffp-contract=off -fno-fast-math), and every step runs there under
# nice -n 10.  After each step every file it wrote is compared with the deposited file at the same place
# (tools/compare.py: exact, apart from run times, dates and paths); the printed output of the steps whose output was
# deposited as a log is compared with that log.  A step fails if it exits non-zero or prints FAIL, NOT CERTIFIED,
# a non-zero failure count or a traceback.
#
# Steps (quick mode runs those marked *):
#   * build        the interval engine (check/cert/compile.sh), ub (check/ubound), wsep and cover (check/wsep), and
#                  libm2f.dylib (PartnersRivalsM2f/cannon/m2f.c), which tests (b) compares with, in the working copy
#                  (no compiled code is deposited)
#   * tests        check/cert/tests.py: quick b c e (engine against exact rationals, the corner lemma L3, the
#                  certificate primitives); full a b c d e t
#   * box_<b>      check/cert/results/run_<b>.py, b = W N E S Snew11 Snew15 Swit W_large cube5 (SI Box 2, SI Table 2);
#                  quick: box_W (the W box) and box_cube5 (the whole 5-cube of wedge W)
#   * family_<f>   check/families/<f>.py, f = cert_pt (quick: only this) cert_ab boundary_S
#   * witnesses    check/cert/witnesses/make_witnesses.py: all 25 (quick: W_SI2, N_SI2, pt_t1_2, MB8_000,
#                  fig3a_right_centre, ab_sliver); each margins file against the deposited one
#   * witness_sha  each regenerated margins file against the SHA-256 recorded in its deposited witness
#     certify8N    check/certify8N.py (the N witness, exact Fractions, against the deposited log)
#   * exclusions   check/ubound/run_all.sh all (quick: the exact classes, the 12-co-player certificate, its exact rows,
#                  the census consistency)
#   * separation   check/wsep/run_all.sh (quick: the negative control, class c9a with its re-check, the lemma check);
#                  the full step needs the Zenodo archive PRFull-1.0.0_certificates.tar.gz and is skipped without it
#   * bingap       cannon/bingap: all 64 tasks of the census of binary pairs and collect.py (quick: task 0 of the
#                  Cannon probe, the focal genomes p < 8192, against the probe's deposited output)
#   * m1_boundary  check/ref3/eq3_boundary/P1/prove_strata.py: the partner conditions on the 13 boundary strata
#     m1_cross     P1/crosscheck.py: its limits against PartnersRivalsM1f/check/exact_partners.py (832 comparisons)
#     m1_p2        P2 (as run_all.sh): the second proof, 245 exact limits, the root search; its three output files
#                  against the deposited ones (its printed factor lists follow the hash seed)
#   * m1_compare   referee/referee_check.py: the two proofs compared, Eq. (3) at 38 514 games
#   * m1_box       box/mk_strata_table.py: SI Box 3 from P1 (without the paper only its body, si_box3_strata_body.tex),
#                  read back against exact_partners.py (check_box_vs_exact_partners.py)
#     m1_dim       stable_dim: facets_cert.py, zd_check.py, f2_display_check.py, zd_derivation.py
#     m1_sets      check/ref3/stable_sets: certify.py (24 jobs; it prints FAIL for D1, the one it leaves open) and
#                  d1_switchline.py (D1 proved)
#     m1_audit     stable_sets: certify_edges.py, normlead.py, audit.py (audit_summary.txt) and audit3.py
#   * m1_cycle     check/ref3/cycle_family/prove.py: the family (0,p_CD,p_DC,1) on its nine strata
#     m1_cycle_x   cycle_family: crosscheck.py (405 769 exact tests) and counterexample.py
#     offer_screen check/ref3/eq4_example/offer/screen_all.py at x = 1/100, 1/4, 2/5, 5/13 (the margins_*.txt)
#   * offer        offer/interval.py and check_eight.py: the offer range 0 < x <= 5/13 from the exact margins, the
#                  eight binding co-players; full mode also renewal, lengths, binders, two_class, selfplay, xcheck
#   * rcheck       cannon/sweep13/validate/rcheck.sh: the replay of the efficiency test of the rarity sweep (rcheck.c)
#                  against rc.tsv, flags.tsv and the Cannon ALPHA lines (quick: tasks 1 and 3; full: all 512, about
#                  17 min on 3 cores: 5.8 s of CPU per task, measured on 2026-11-03 with the x86_64 binary the step
#                  builds when the toolchain runs under Rosetta)
#   * stratcert    cannon/stratcert/collect.py on the 64 deposited task outputs of the certified points of the
#                  stratum counts (Cannon 48975432 and 48981870): completeness, the audit of every line, the exact
#                  self-play of all 12 748 witnesses recomputed; verdict True; against the deposited collect.out
#   * evoaudit     the exact audit of the evolutionary atoms (Cannon 48975385 and 48981982): the dumps it samples
#                  taken from PartnersRivalsM2f/cannon/archive (tools/evoaudit_dumps.py), then src/analyze.py on the
#                  deposited chunks; 27 261 audited, 2 changed atoms, 0 unresolved; against the deposited results/
#                  (the same two steps in both modes: they check the deposited cluster output, which is not rerun)
#
# Exit status (of the steps this run ran): 0 every step passed and every rewritten file is the deposited one; 2 every
# step passed but some file differs; 1 a step failed, or an --only name is no step of the mode.  The summary is
# DIR/logs/certificates/summary.txt.

set -u
REL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRF="$(dirname "$REL")"
ROOT="$(dirname "$PRF")"
PY="${PYTHON:-python3}"
CC="${CC:-cc}"
MODE=full; OUT="$REL/work"; NTH=3; BINGAP=1; LIST=0; ONLY=""
while [ $# -gt 0 ]; do
  case "$1" in
    --quick) MODE=quick ;;
    --out) OUT="$2"; shift ;;
    --threads) NTH="$2"; shift ;;
    --no-bingap) BINGAP=0 ;;
    --only) ONLY=",$2,"; shift ;;
    --list) LIST=1 ;;
    -h|--help) sed -n '2,84p' "$0"; exit 0 ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
  shift
done
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
W="$OUT/certificates"; WP="$W/PRFull"; LOG="$OUT/logs/certificates"; SUM="$LOG/summary.jsonl"
export CERT_MAX_THREADS="$NTH" PYTHONDONTWRITEBYTECODE=1 CC

if [ "$LIST" = 1 ]; then grep -E '^#   [ *] [a-z][A-Za-z0-9_<>-]* ' "$0" | sed 's/^#   //'; exit 0; fi

command -v "$PY" >/dev/null 2>&1 || { echo "no Python interpreter '$PY'"; exit 1; }
"$PY" -c "import numpy, gmpy2, sympy" 2>/dev/null || { echo "$PY cannot import numpy, gmpy2 and sympy (see PRFull/release/RELEASE.md)"; exit 1; }
command -v "$CC" >/dev/null 2>&1 || { echo "no C compiler '$CC'"; exit 1; }
T0=$(date +%s)
echo "reproduce_certificates.sh ($MODE${ONLY:+, only ${ONLY//,/ }}) in $ROOT"
echo "  $("$PY" -c 'import sys, numpy, gmpy2, sympy; print("Python " + sys.version.split()[0] + ", numpy " + numpy.__version__ + ", gmpy2 " + gmpy2.version() + ", sympy " + sympy.__version__)'), $("$CC" --version 2>&1 | head -1)"
echo "  working copy $W, logs $LOG, $NTH threads"

# the tree's own absolute paths (a no-op in the author's tree), then a fresh working copy (--only reuses the last one)
"$PY" "$REL/relocate.py" | tail -1
if [ -z "$ONLY" ] || [ ! -f "$W/.stamp" ]; then
  rm -rf "$W" "$LOG"; mkdir -p "$WP/cannon" "$LOG"
  rsync -a --exclude __pycache__ --exclude .DS_Store --exclude '*.dylib' --exclude '*.so' \
        --exclude '/ubound/ub' --exclude '/wsep/wsep' --exclude '/wsep/cover' "$PRF/check/" "$WP/check/"
  rsync -a --exclude /bingap --exclude logs --exclude probe_out "$PRF/cannon/bingap/" "$WP/cannon/bingap/" 2>/dev/null || true
  KX=(--exclude __pycache__ --exclude .DS_Store --exclude '*.dylib' --exclude '*.so' --exclude /logs --exclude /venv --exclude /pylib)
  rsync -a "${KX[@]}" --exclude /prep --exclude /collect.out --exclude /out/stratcert_all.tsv \
        "$PRF/cannon/stratcert/" "$WP/cannon/stratcert/"
  rsync -a "${KX[@]}" --exclude /data --exclude /results --exclude /results_probe --exclude /src/review \
        "$PRF/cannon/evoaudit/" "$WP/cannon/evoaudit/"
  # the replay of the efficiency test (step rcheck): the kit's validate/ and the ALPHA lines of the Cannon run
  mkdir -p "$WP/cannon/sweep13/out_13"
  rsync -a "${KX[@]}" --exclude /build "$PRF/cannon/sweep13/validate/" "$WP/cannon/sweep13/validate/"
  rsync -a --include 'a_1o3_*.tsv' --exclude '*' "$PRF/cannon/sweep13/out_13/" "$WP/cannon/sweep13/out_13/"
  if [ -f "$PRF/Text/si_box3_legend.tex" ]; then     # the title and legend of SI Box 3 (part of the paper: release 1.1.0)
    mkdir -p "$WP/Text"; cp -p "$PRF/Text/si_box3_legend.tex" "$WP/Text/"
  fi
  "$PY" "$REL/relocate.py" --map "$PRF/check=$WP/check" --map "$PRF/cannon/bingap=$WP/cannon/bingap" \
        --map "$PRF/cannon/stratcert=$WP/cannon/stratcert" --map "$PRF/cannon/evoaudit=$WP/cannon/evoaudit" \
        --map "$PRF/cannon/sweep13=$WP/cannon/sweep13" "$WP" | sed 's/^/  /'
  touch "$W/.stamp"
fi
# libm2f, which tests (b) compares the engine with, is built by the step build in the working copy, beside a copy of its
# wrapper m2f.py (which finds it at ../cannon/libm2f.dylib); tests.py is pointed at that copy, not at the deposit
if [ ! -f "$W/PartnersRivalsM2f/check/m2f.py" ]; then
  mkdir -p "$W/PartnersRivalsM2f/check" "$W/PartnersRivalsM2f/cannon"
  cp -p "$ROOT/PartnersRivalsM2f/cannon/m2f.c" "$ROOT/PartnersRivalsM2f/cannon/m2f_py.c" "$W/PartnersRivalsM2f/cannon/"
  "$PY" "$REL/relocate.py" --map "$ROOT/PartnersRivalsM2f/check=$W/PartnersRivalsM2f/check" "$WP/check/cert/tests.py" \
        | sed 's/^/  /'
  cp -p "$ROOT/PartnersRivalsM2f/check/m2f.py" "$W/PartnersRivalsM2f/check/"
fi
mkdir -p "$LOG"

NFAIL=0
# the deposited outputs were written under the author's tree; paths under it in them are compared as paths, so that the
# release verifies wherever it is unpacked
AUTHOR_ROOT=/Users/martin/Documents/CL1
SEEN=","; RUNID="$MODE.$T0.$$"    # the steps of this mode (for --only); this run, recorded with each of its steps
step() {   # step NAME LABEL DIR [runstep options] -- COMMAND ...
  local name=$1 label=$2 dir=$3; shift 3
  SEEN="$SEEN$name,"
  if [ -n "$ONLY" ]; then case "$ONLY" in *",$name,"*) ;; *) return 0 ;; esac; fi
  "$PY" "$REL/tools/runstep.py" --name "$name" --label "$label" --cwd "$dir" --logdir "$LOG" --summary "$SUM" \
      --run-id "$RUNID" --strip "$WP" --strip "$PRF" --strip "$ROOT" --strip "$AUTHOR_ROOT/PRFull" --strip "$AUTHOR_ROOT" "$@"
  local st=$?
  [ $st -eq 1 ] && NFAIL=$((NFAIL + 1))
  return 0
}
SC=(--scan "$WP/check=$PRF/check")
C="$WP/check"

# ---------------------------------------------------------------------------------------------------------- build
M2C="$W/PartnersRivalsM2f/cannon"; ARCH=""; [ "$(uname)" = Darwin ] && ARCH="-arch arm64 -arch x86_64"
step build "build the interval engine, ub, wsep and cover, and libm2f in the working copy (tests (b) compares with it)" "$C" "${SC[@]}" -- bash -c "
  set -e; ./cert/compile.sh
  (cd $M2C && $CC -O3 $ARCH -fPIC -shared -o libm2f.dylib m2f.c m2f_py.c -lm && echo built libm2f.dylib)
  $CC -O2 -ffp-contract=off -fno-fast-math -o ubound/ub ubound/ub.c -lm -lpthread && cmp ubound/certeng.c cert/certeng.c
  ./wsep/compile.sh && $CC -O2 -ffp-contract=off -fno-fast-math -o wsep/cover wsep/cover.c -lm -lpthread
  echo built"

# ---------------------------------------------------------------------------------------------------------- tests
if [ "$MODE" = quick ]; then
  step tests "the certificate library: engine against exact rationals (b), the corner lemma L3 (c), primitives (e)" \
       "$C" "${SC[@]}" --expect '==== [0-9]+ checks, 0 failed' -- "$PY" cert/tests.py b c e
else
  step tests "the certificate library, all sections (a b c d e t)" "$C" "${SC[@]}" --expect '==== [0-9]+ checks, 0 failed' \
       -- "$PY" cert/tests.py
fi

# ---------------------------------------------------------------------------------------------------------- boxes
if [ "$MODE" = quick ]; then BOXES="W cube5"; else BOXES="W N E S Snew11 Snew15 Swit W_large cube5"; fi
for b in $BOXES; do
  step "box_$b" "certified box or stratum: cert/results/run_$b.py" "$C" "${SC[@]}" -- "$PY" "cert/results/run_$b.py"
done

# ---------------------------------------------------------------------------------------------------------- families
if [ "$MODE" = quick ]; then FAM="cert_pt"; else FAM="cert_pt cert_ab boundary_S"; fi
for f in $FAM; do
  step "family_$f" "families/$f.py" "$C" "${SC[@]}" --expect 'CERTIFIED|not exhibited' --stdout-vs "$PRF/check/families/$f.out" \
       -- "$PY" "families/$f.py"
done

# ---------------------------------------------------------------------------------------------------------- witnesses
if [ "$MODE" = quick ]; then
  step witnesses "six of the 25 exact witnesses (cert/witnesses/make_witnesses.py)" "$C" "${SC[@]}" --expect '^done: 6 witnesses' \
       -- "$PY" cert/witnesses/make_witnesses.py W_SI2 N_SI2 pt_t1_2 MB8_000 fig3a_right_centre ab_sliver
else
  step witnesses "all 25 exact witnesses (cert/witnesses/make_witnesses.py)" "$C" "${SC[@]}" --expect '^done: 25 witnesses' \
       -- "$PY" cert/witnesses/make_witnesses.py
fi
step witness_sha "every regenerated margins file against the SHA-256 recorded in its deposited witness" "$C" -- "$PY" - <<EOF
import gzip, hashlib, json, os, sys
new, dep = "$C/cert/witnesses", "$PRF/check/cert/witnesses"
bad = n = 0
for f in sorted(os.listdir(new)):
    if not f.endswith(".json") or f == "index.json": continue
    m = os.path.join(new, "margins", f[:-5] + ".tsv.gz")
    if not os.path.exists(m) or os.path.getmtime(m) < os.path.getmtime("$W/.stamp"): continue      # not regenerated
    J = json.load(open(os.path.join(dep, f)))
    sha = hashlib.sha256(gzip.open(m).read()).hexdigest(); n += 1
    ok = sha == J["margins_sha256"] and json.load(open(os.path.join(new, f)))["recomputed"] == J["recomputed"]
    bad += not ok
    print("%-20s margins sha256 %s  %s" % (f[:-5], sha[:16], "matches the deposit" if ok else "DIFFERS"))
print("%d regenerated witnesses, %d differ from the deposit" % (n, bad)); sys.exit(1 if bad or not n else 0)
EOF

# ---------------------------------------------------------------------------------------------------------- certify8N
if [ "$MODE" = full ]; then
  step certify8N "the N witness in exact Fractions (check/certify8N.py)" "$C" "${SC[@]}" --expect 'FRIENDLY RIVAL in wedge N: True' \
       --stdout-vs "$PRF/check/certify8N.out" -- "$PY" certify8N.py
fi

# ---------------------------------------------------------------------------------------------------------- exclusions
if [ "$MODE" = quick ]; then
  step exclusions "exclusions: exact classes, the 12-co-player certificate, its exact rows, census consistency" "$C/ubound" "${SC[@]}" \
       --expect 'classes: none 31580692 C 6206077 A 2788386 Z 2471566' --expect 'class mismatches with the exact classes here: 0' -- bash -c "
    set -e; mkdir -p runs results
    ./ub classes cls.bin $NTH | tee results/classes.out
    $PY counts.py | tee results/counts.out
    rm -f runs/x.idx runs/x.qp runs/x.qm runs/x.log
    ./ub cover cls.bin qstar.txt runs/x $NTH | tee results/cover12.out
    $PY status.py runs/x | tee results/status_cover12.out
    $PY exact_rows.py qstar.txt | tee results/exact_rows.out
    $PY census_check.py runs/x | tee results/census_check_cover12.out"
else
  step exclusions "exclusions: check/ubound/run_all.sh all (the certificate, its re-check, the by-product)" "$C/ubound" "${SC[@]}" \
       --expect 'class mismatches with the exact classes here: 0' -- env PYTHON="$PY" bash run_all.sh all
fi

# ---------------------------------------------------------------------------------------------------------- separation
# wsep names, for each excluded stratum, the first co-player that works; with 3 threads the pool order and hence the
# co-player in a row depend on the schedule, so its tables are compared as certificates (tools/wsep_same.py: the
# same strata, counts and verdicts), not byte for byte.  The negative control must leave all 364 strata standing.
WS=(--nocounts --ignore '^wsep/out/'
    --expect '^NOT_EXCLUDED 364$' --expect '"rows_failing": 0' --expect '"multilinear_bad": \[\]' --expect '^finite eps ok: True'
    --expect 'control binary 19079 .*certified FR of W: True' --expect 'control eq:certW .*certified FR of W: True'
    --expect 'wsep classes compared: [0-9]+, differing: 0')
if [ "$MODE" = quick ]; then
  step separation "separation in W: negative control, class c9a and its re-check, the lemma check" "$C/wsep" "${SC[@]}" "${WS[@]}" -- bash -c "
    set -e; mkdir -p out/extra out/cover
    ./wsep file:ctrl_W364.txt $NTH ctrl_W364.out > ctrl_W364.summary; cat ctrl_W364.summary
    ./wsep c9a $NTH out/c9a.tsv > out/c9a.summary; cat out/c9a.summary
    $PY verify.py $NTH c9a
    $PY lemma_check.py
    $PY $REL/tools/wsep_same.py out $PRF/check/wsep/out c9a; true"   # not a lone command last (peak memory, see evoaudit)
else
  # the deposited runs of the eight extra classes, which wsep_same.py compares with, are in the Zenodo record's archive
  # PRFull-1.0.0_certificates.tar.gz, not in the GitHub repository; without them the step is skipped and says so
  NEEDW=(--need-what "they are in the archive PRFull-1.0.0_certificates.tar.gz of the Zenodo record; bash PRFull/release/fetch_zenodo.sh certificates fetches it")
  for x in x4 x13 x5 x15 f1p11 f1p14 f2p11 f2p14; do
    NEEDW+=(--need-file "$PRF/check/wsep/out/extra/$x.summary" --need-file "$PRF/check/wsep/out/extra/$x.tsv")
  done
  step separation "separation in W: check/wsep/run_all.sh (with the greedy cover)" "$C/wsep" "${NEEDW[@]}" "${SC[@]}" "${WS[@]}" -- bash -c "
    set -e; PYTHON=$PY bash run_all.sh; cat ctrl_W364.summary out/*.summary
    $PY $REL/tools/wsep_same.py out $PRF/check/wsep/out c6 c8 c68 c9a c9b c0 extra/x4 extra/x13 extra/x5 extra/x15 \
        extra/f1p11 extra/f1p14 extra/f2p11 extra/f2p14; true"
fi

# ---------------------------------------------------------------------------------------------------------- bingap
if [ "$BINGAP" = 1 ]; then
  B="$WP/cannon/bingap"
  BS=(--scan "$B=$PRF/cannon/bingap")
  if [ "$MODE" = quick ]; then
    step bingap "census of binary pairs: the probe's task 0 (7.9 M pairs, p < 8192, as Cannon's probe 48916871)" "$B" "${BS[@]}" \
         --must-write probe_out/task_0.txt -- bash -c "
      set -e; $CC -O2 -ffp-contract=off -fno-fast-math -o bingap bingap.c certeng.c -lm -lpthread
      mkdir -p probe_out; ./bingap 0 64 --checksym --pmax 8192 > probe_out/task_0.txt 2> probe_out/task_0.err
      cat probe_out/task_0.txt; test ! -s probe_out/task_0.err"
  else
    step bingap "census of binary pairs: all 64 tasks ($NTH at a time) and collect.py" "$B" "${BS[@]}" \
         --must-write results.txt --expect 'smallest nonzero .* 1/4140 ' -- bash -c "
      set -e; $CC -O2 -ffp-contract=off -fno-fast-math -o bingap bingap.c certeng.c -lm -lpthread
      mkdir -p out; seq 0 63 | xargs -P $NTH -I{} env SLURM_ARRAY_TASK_ID={} NTASKS=64 bash task.sh
      for f in out/task_*.err; do test ! -s \$f || { echo FAIL: \$f not empty; exit 1; }; done
      $PY collect.py > results.txt; cat results.txt"
  fi
fi

# ---------------------------------------------------------------------------------------------------------- ref3
# The exact checks of check/ref3 (sympy over Q, or Fractions); each script runs in its own directory.  The printed
# output of a script whose log is deposited is written to that log in the working copy (tee), so that the comparison
# after the step checks it; the other files they write are compared as usual.
R3="$C/ref3"; E3="$R3/eq3_boundary"
step m1_boundary "SI Box 3: the partner conditions on the 13 boundary strata (ref3/eq3_boundary/P1/prove_strata.py)" "$E3/P1" \
     "${SC[@]}" --must-write prove_strata.json --must-write prove_strata.log --expect '^  OVERALL PROVED: True ' \
     -- "$PY" prove_strata.py
if [ "$MODE" = full ]; then
  step m1_cross "the limits of P1 against PartnersRivalsM1f/check/exact_partners.py (P1/crosscheck.py)" "$E3/P1" "${SC[@]}" \
       --stdout-vs "$PRF/check/ref3/eq3_boundary/P1/crosscheck.log" \
       --expect '^cross-check passed: 832 Gaussian eps=1e-40 comparisons, 832 exact comparisons with exact_partners\.limitW$' \
       -- "$PY" crosscheck.py
  # as run_all.sh, with $PY; its printed lists of sympy factors come in an order that depends on the hash seed, so its
  # three output files are compared with the deposited ones, not its printed output with run_all.log
  step m1_p2 "the second proof of the 13 strata (P2: p2_proof, p2_numcheck, p2_rootsearch, p2_tables)" "$E3/P2" "${SC[@]}" \
       --must-write p2_summary.json --must-write p2_proof_report.txt --must-write p2_polys.pkl \
       --expect '^ALL CERTIFICATES PASSED; ' --expect '^points 245, failures 0, ' \
       --expect '^TOTAL factors with a zero in the open stratum: 0 ' -- bash -c "
    set -e; $PY p2_proof.py; $PY p2_numcheck.py; $PY p2_rootsearch.py; $PY p2_tables.py; true"
fi
step m1_compare "the two proofs of SI Box 3 compared; Eq. (3) at 38 514 games (referee/referee_check.py)" "$E3/referee" \
     "${SC[@]}" --stdout-vs "$PRF/check/ref3/eq3_boundary/referee/referee_check.log" --expect '^problems: 0$' \
     -- "$PY" referee_check.py
# without the paper (release 1.0.0) the title and legend of SI Box 3 are absent, and the box is written without them
if [ -f "$WP/Text/si_box3_legend.tex" ]; then
  BOX=(--must-write si_box3_strata.tex --expect '^wrote si_box3_strata\.tex: the same with its title and legend$')
else
  BOX=(--expect '^si_box3_strata\.tex not written: its title and legend \(Text/si_box3_legend\.tex\) are part of the paper')
fi
step m1_box "SI Box 3 written from P1 and read back against exact_partners.py (box/mk_strata_table.py, check_box_vs_exact_partners.py)" \
     "$E3/box" "${SC[@]}" --must-write si_box3_strata_body.tex --must-write check_box.log "${BOX[@]}" \
     --expect '^wrote si_box3_strata_body\.tex: 13 strata, 34 gains that differ from the display, 65 Bernstein factors$' \
     --expect '^typeset box agrees with exact_partners\.limitW: 204 listed gains, 1044 display gains, 78 L>0 checks; 0 mismatches$' \
     -- bash -c "
  set -eo pipefail; $PY mk_strata_table.py --out si_box3_strata.tex
  $PY check_box_vs_exact_partners.py si_box3_strata_body.tex | tee check_box.log; true"
if [ "$MODE" = full ]; then
  step m1_dim "the open facets and the zero-determinant plane (eq3_boundary/stable_dim, four scripts)" "$E3/stable_dim" "${SC[@]}" \
       --must-write facets_cert.log --must-write zd_check.log --must-write f2_display_check.log --must-write zd_derivation.log \
       --expect '^SUMMARY  F2 reduction: True;  F3-F8 all states recur .*: \[True, True, True, True, True, True\];  alternator reduction: True ' \
       --expect '^interior points of E: 200, all 16 deviators tie exactly; 200 nearby points off E: not all tie$' \
       --expect '^proposed p_DD = 0 display agrees with exact_partners\.limitW at 960 \(point, game, deviator\) triples$' \
       --expect '^alpha = beta = 0 at all four states  <=>  p on the plane E .*: confirmed$' -- bash -c "
    set -eo pipefail; for s in facets_cert zd_check f2_display_check zd_derivation; do $PY \$s.py | tee \$s.log; done; true"
  # certify.py prints FAIL for the one job it leaves open (D1, (1,1,0,y) on the switch line), which d1_switchline.py
  # then proves: the fail words are not read here; the expected lines and the comparison of certify.log decide
  step m1_sets "the other memory-one stability statements: ref3/stable_sets/certify.py (24 jobs) and d1_switchline.py" \
       "$R3/stable_sets" "${SC[@]}" --no-fail-words --must-write certify.log --must-write d1_switchline.log \
       --expect '^FAIL D1 110o ' --expect '^jobs 24 with problems 1$' --expect '^D1 proved: .*: True$' -- bash -c "
    set -eo pipefail; $PY certify.py | tee certify.log; $PY d1_switchline.py | tee d1_switchline.log; true"
  # the claim table of audit.py (audit_summary.txt) counts failures in a column that no fail word reads: its sum is printed
  step m1_audit "the exact sweep over all 81 strata (stable_sets: certify_edges, normlead, audit, audit3)" "$R3/stable_sets" \
       "${SC[@]}" --must-write certify_edges.log --must-write normlead.log --must-write audit_summary.txt --must-write audit3.log \
       --expect '^1o10 V \[1, 0, 0, 0\] binding constraints .* Farkas failures \[\] ' \
       --expect '^11o0 V \[1, 0, 0, 0\] binding constraints .* Farkas failures \[\] ' \
       --expect '^1110 V \[1, 0, 0, 0\] binding constraints .* Farkas failures \[\] ' --expect '^strata with a lead of uncertain sign: 0 of 81$' \
       --expect '^== 3\. efficiency closed form mismatches: 0$' --expect '^== 4\. games without a partner among 457 candidates: 0 of 1385 \[\]$' \
       --expect '^audit\.py claims: 9, failures: 0$' -- bash -c "
    set -eo pipefail; $PY certify_edges.py | tee certify_edges.log; $PY normlead.py | tee normlead.log
    $PY audit.py | tee audit_summary.txt; $PY audit3.py | tee audit3.log
    awk 'NR == 1 && \$1 == \"claim\" {t = 1; next} t && NF == 0 {exit} t {r++; n += \$3}
         END {print \"audit.py claims: \" r \", failures: \" n + 0}' audit_summary.txt; true"
fi
step m1_cycle "the family (0,p_CD,p_DC,1) on its nine strata (ref3/cycle_family/prove.py)" "$R3/cycle_family" "${SC[@]}" \
     --must-write prove.json --stdout-vs "$PRF/check/ref3/cycle_family/prove.log" \
     --expect '^positive-dimensional strata all match the claim: True$' -- "$PY" prove.py
if [ "$MODE" = full ]; then
  step m1_cycle_x "the family (0,p_CD,p_DC,1): 405 769 exact tests and the counterexample (cycle_family)" "$R3/cycle_family" \
       "${SC[@]}" --must-write crosscheck.log --must-write counterexample.log \
       --expect '^A\. 405769 exact \(strategy, game\) tests: corrected predicate disagrees 0 times; ' \
       --expect '^self-play uniform for every eps: True$' -- bash -c "
    set -eo pipefail; $PY crosscheck.py | tee crosscheck.log; $PY counterexample.py | tee counterexample.log; true"
  step offer_screen "Eq. (4): the exact margins of all 65536 binary co-players at x = 1/100, 1/4, 2/5, 5/13 (offer/screen_all.py)" \
       "$R3/eq4_example/offer" "${SC[@]}" --must-write margins_1_100.txt --must-write margins_1_4.txt \
       --must-write margins_2_5.txt --must-write margins_5_13.txt --must-write screen_all.log \
       --expect '^x = 2/5: 65536 co-players, min margin -1/198 attained by 8; negative: 400; zero: 38068; ' \
       --expect '^x = 5/13: 65536 co-players, min margin 0 attained by 38396; negative: 0; zero: 38396; ' -- bash -c "
    set -eo pipefail; for x in '1 100' '1 4' '2 5' '5 13'; do $PY screen_all.py \$x; done | tee screen_all.log; true"
fi
OFFER="interval check_eight"
OFX=(--expect '^identically zero: 37996 ; constant non-zero: 8172 \(negative: 0 \)$'   # a constant negative margin: an empty range
     --expect '^numerator negative as x -> 0\+: 0 \[\]$' --expect "^roots of numerators inside \(0,1\): \{'5/13': 400, "
     --expect '^reconstruction verified at x = 7/10 ' --expect '^all 21 x-values match \(5-13x\)/\(2\(25-13x\)\) for all eight; ')
if [ "$MODE" = full ]; then
  OFFER="$OFFER renewal lengths binders two_class selfplay xcheck"
  OFX+=(--expect 'renewal margin 7/174 = exact limit 7/174$' --expect '^zero at 5/13 but not at 1/4: 400 ; negative at 2/5: 400 ; same set: True$'
        --expect '^400 binders; per-branch gains \+1 and -5/8 for all; min offer length 3 min no-offer length 25/4$'   # "3 min" reads as a time
        --expect '^   exact: 1 \+ \(1 - x/2\)/\(x/2\) = 8 ; offers before recovery = x \* 2/x = 2$'
        --expect '^sympy exact: co-player 33153 margin\(x\) = -\(13\*x - 5\)/\(11\*x \+ 201\)  root: \[5/13\]$')
fi
step offer "Eq. (4): the offer range 0 < x <= 5/13 from the exact margins; the eight binding co-players (ref3/eq4_example/offer)" \
     "$R3/eq4_example/offer" "${SC[@]}" "${OFX[@]}" -- bash -c "set -eo pipefail; for s in $OFFER; do $PY \$s.py | tee \$s.log; done; true"

# ---------------------------------------------------------------------------------------------------------- rcheck
# The replay of the efficiency test of the rarity sweep at alpha = 1/3 (SI Box 1; Computational Supplement): each task's
# draws regenerated from its seed and their self-play class read twice (cannon/sweep13/validate/rcheck.sh, rcheck.c),
# against the deposited per-task counts rc.tsv, the misread draws flags.tsv and the Cannon run's ALPHA lines
if [ "$MODE" = quick ]; then RCT="${RCHECK_TASKS:-1 3}"; else RCT="${RCHECK_TASKS:-$(seq 0 511 | tr '\n' ' ')}"; fi
SV="$WP/cannon/sweep13/validate"
step rcheck "the replay of the efficiency test (cannon/sweep13/validate/rcheck.c), $(echo $RCT | wc -w | tr -d ' ') task(s), $NTH at a time" \
     "$SV" --scan "$SV=$PRF/cannon/sweep13/validate" \
     --expect "^replayed $(echo $RCT | wc -w | tr -d ' ') tasks: RC lines differ 0, FLAG lines differ 0, ALPHA counts differ 0\$" \
     -- env JOBS="$NTH" PYTHON="$PY" bash rcheck.sh $RCT

# ---------------------------------------------------------------------------------------------------------- stratcert
# The cluster run is not repeated (30.8 core-hours); its 64 task outputs are checked again as deposited: every retained
# stratum has one line, every line is audited, and the exact self-play limit of every witness is recomputed.
SK="$WP/cannon/stratcert"
step stratcert "certified points of the stratum counts: cannon/stratcert/collect.py on the 64 deposited task outputs" "$SK" \
     --scan "$SK=$PRF/cannon/stratcert" --must-write out/stratcert_all.tsv --stdout-vs "$PRF/cannon/stratcert/collect.out" \
     --expect '^units 12748, lines 12748, duplicated 0, malformed 0, missing 0 \[\]$' --expect '^tasks without \.done \(of 64\): 0 \[\]$' \
     --expect '^audit of 12748 certified lines \(exact self-play recomputed for 12748\): no problem$' \
     --expect '^verdict: every retained stratum of W, N, E and ray\+ has a certified exact friendly rival: True$' \
     -- "$PY" collect.py out 64

# ---------------------------------------------------------------------------------------------------------- evoaudit
# The audit itself (10.7 core-hours) is not repeated; its 128 deposited chunks are analysed again.  analyze.py redraws
# the design from the dumps and asserts that it gives the deposited sample; the dumps are taken from the deposited
# archives (0.67 GB in the working copy while the step runs, removed after it).
EK="$WP/cannon/evoaudit"
# its inputs that the GitHub repository does not carry are in the Zenodo record's archive
# PRFull-1.0.0_m2-archives-and-audit.tar.gz (fetch_zenodo.sh); without them the step is skipped and says so
NEEDE=(--need-file "$ROOT/PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz"
       --need-file "$ROOT/PartnersRivalsM2f/cannon/archive/cannon_out_f2_2026-09-26.tgz"
       --need-file "$PRF/cannon/evoaudit/out" --need-file "$PRF/cannon/evoaudit/out_probe"
       --need-what "they are in the archive PRFull-1.0.0_m2-archives-and-audit.tar.gz of the Zenodo record; bash PRFull/release/fetch_zenodo.sh fetches it")
step evoaudit "exact audit of the evolutionary atoms: the dumps from the archives, then src/analyze.py on the deposited chunks" "$EK" "${NEEDE[@]}" \
     --scan "$EK=$PRF/cannon/evoaudit" --ignore '^data/' --must-write results/summary.json --must-write results/confusion.txt \
     --must-write results/records.tsv.gz --must-write results/analyze.out \
     --expect '^dumps: 1024 files, 4772732 records' --expect '^sampled strategies equal their records in the dumps: 27261 of 27261$' \
     --expect '^audit: 27261 of 27261 sampled strategies audited, 0 errors, 0 missing, 2 changed atoms, 0 not reproduced$' \
     --expect '^audited 27261, errors 0, missing 0, unresolved 0; changed atoms 2 \(010->000 2\); changed bits: E 0, C 0, S 2; ' \
     -- bash -c "
    set -e; rm -rf data results
    $PY $REL/tools/evoaudit_dumps.py extract $ROOT/PartnersRivalsM2f/cannon/archive .
    mkdir -p results; $PY src/analyze.py out_probe,out results > results/analyze.out; head -3 results/analyze.out
    $PY $REL/tools/evoaudit_dumps.py summary results
    rm -rf data; echo 'the extracted dumps are removed from the working copy'"   # a builtin last: bash would exec a last
                                                                             # external command, and its peak memory is lost

# ---------------------------------------------------------------------------------------------------------- summary
UNKNOWN=""                       # --only names that are no step of this mode
for n in ${ONLY//,/ }; do case "$SEEN" in *",$n,"*) ;; *) UNKNOWN="$UNKNOWN $n" ;; esac; done
"$PY" - "$SUM" "$LOG/summary.txt" "$MODE" "$(( $(date +%s) - T0 ))" "$RUNID" "$UNKNOWN" <<'EOF'
import json, os, sys, time
recs = [json.loads(l) for l in open(sys.argv[1])] if os.path.exists(sys.argv[1]) else []
run, unknown = sys.argv[5], sys.argv[6].split()
# the steps of this run; with --only, the latest record of every other step of the working copy is listed apart
this = [r for r in recs if r.get("run") == run]
last = {r["name"]: i for i, r in enumerate(recs)}
old = [r for i, r in enumerate(recs) if last[r["name"]] == i and r.get("run") != run]
def row(r):
    return "%-14s %-7s %8s %9s %7s %6s %7s %9s" % (r["name"], r["status"], r.get("seconds", "-"), r.get("peak_rss_mib", "-"),
           r.get("written", "-"), r.get("same", "-"), r.get("differ", "-"), r.get("new", "-"))
def when(r):
    m = r.get("run", "").split(".")
    return "--%s run of %s" % (m[0], time.strftime("%Y-%m-%d %H:%M", time.localtime(int(m[1])))) if len(m) == 3 else "an earlier run"
recs = this
out = ["reproduce_certificates.sh --%s: %d steps in %s s" % (sys.argv[3], len(recs), sys.argv[4])]
out.append("%-14s %-7s %8s %9s %7s %6s %7s %9s" % ("step", "status", "seconds", "peak MiB", "written", "same", "differ", "unchecked"))
out += [row(r) for r in recs]
nf = sum(r["status"] == "failed" for r in recs); nd = sum(r.get("differ", 0) for r in recs)
ns = sum(r["status"] == "skipped" for r in recs)
out.append("steps failed: %d; skipped: %d; files that differ from the deposit: %d; largest peak memory %.0f MiB"
           % (nf, ns, nd, max([r.get("peak_rss_mib", 0) for r in recs] or [0])))
for r in recs:
    if r["status"] == "skipped":
        out.append("  skipped: %s: %s" % (r["name"], r.get("why", "")))
    for d in r.get("differing", []):
        out.append("  differs: %s: %s" % (r["name"], d))
if old:
    out.append("earlier records of this working copy, not counted above: %d steps" % len(old))
    out += ["  %s  (%s)" % (row(r), when(r)) for r in old]
if unknown:
    out.append("--only: no step named %s in --%s mode (see --list); nothing ran under that name"
               % (", ".join(unknown), sys.argv[3]))
open(sys.argv[2], "w").write("\n".join(out) + "\n"); print("\n" + "\n".join(out))
sys.exit(1 if nf or unknown else (2 if nd else 0))
EOF
