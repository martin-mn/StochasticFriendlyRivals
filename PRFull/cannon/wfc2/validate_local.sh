#!/bin/bash
# validate_local.sh -- the switch checks of the wfc2 kit, run on one machine (the Mac or a Cannon login/test node).
#   usage: [JOBS=2] validate_local.sh <workdir> [itend=100000]     (at most JOBS processes at a time, nice 10)
# It copies the kit into <workdir>, builds every binary there (compile.sh) plus the PRODUCTION source ref/wfc_prod.c
# as wfc_prod / wfc2_prod, and checks:
#   D1  default F1 build == production build, REP lines (timings aside) and strategy dumps (cmp), 4 games x 2 replicates
#   D2  default F2 build == production F2 build, the same, 2 games x 2 replicates (itend x 10)
#   D3  neutral mode (kernel draws + classification): NEU lines identical, default vs production
#   T1  TSERIES build: REP lines and strategy dumps identical to production, TS/TB/TSR lines present and consistent
#   A   every arm binary runs (A1 alpha=1/4, A2 alpha=2/5, K0, B0, B0K0, E3, E5, E3F2, E5F2) and prints its REP lines
# Exit status 0 only if every check passes.  Nothing is written outside <workdir>.
set -uo pipefail
KIT="$(cd "$(dirname "$0")" && pwd)"
W=${1:?usage: validate_local.sh <workdir> [itend]}; IT=${2:-100000}
mkdir -p "$W" && rsync -a --exclude /out --exclude '/wfc' --exclude '/wfc2' --exclude '/wfc_*' --exclude '/wfc2_*' "$KIT/" "$W/kit/"
cd "$W/kit" || exit 1
CC=${CC:-cc}; CFLAGS=${CFLAGS:--O3 -march=native}; F2="-DNPOP=100 -DBETAV=3.0 -DMUV=1e-4"
nice -n 10 bash compile.sh > build.log 2>&1 || { cat build.log; exit 1; }
$CC $CFLAGS -I. -o wfc_prod  ref/wfc_prod.c m2f.c -lm && $CC $CFLAGS -I. $F2 -o wfc2_prod ref/wfc_prod.c m2f.c -lm || exit 1
$CC $CFLAGS -DLOCALMUT -DK0VERIFY -o wfc_k0v wfc.c m2f.c -lm || exit 1          # K0 with production's W-111 test alongside
mkdir -p v; fail=0
ok()  { echo "PASS $*"; }
bad() { echo "FAIL $*"; fail=1; }
strip() { awk '$1=="REP"{ $25="t"; $26="t"; print }' "$1"; }       # REP lines without t_sim and t_class
stripn() { awk '$1=="NEU"{ $NF="t"; print }' "$1"; }
# ---- D1, T1, and the F1 arms: 4 games (W 30, S 47, N 58, E 165), 2 replicates each, production seeds
jobs=()
for g in 30 47 58 165; do
  jobs+=("wfc_prod p_$g $g 1/3 $IT 1000 4 2 9000000 1" "wfc d_$g $g 1/3 $IT 1000 4 2 9000000 1")
  jobs+=("wfc_ts t_$g $g 1/3 $IT 1000 4 2 9000000 1")
done
for g in 30 47; do
  jobs+=("wfc a1_$g $g 1/4 $IT 1000 4 1 31000000 1" "wfc a2_$g $g 2/5 $IT 1000 4 1 32000000 1")
  jobs+=("wfc_k0 k0_$g $g 1/3 $IT 1000 4 1 33000000 1" "wfc_b0 b0_$g $g 1/3 $IT 1000 4 1 34000000 1")
  jobs+=("wfc_b0k0 b0k0_$g $g 1/3 $IT 1000 4 1 34500000 1" "wfc_e3 e3_$g $g 1/3 $IT 1000 4 1 35000000 1")
  jobs+=("wfc_e5 e5_$g $g 1/3 $IT 1000 4 1 36000000 1")
done
jobs+=("wfc_k0v k0v_467 467 1/3 $IT 1000 4 1 33000000 1")      # W game where K0 grows friendly-rival clouds
# ---- D2 and the F2 arms: N = 100 at 10 x itend
for g in 30 172; do
  jobs+=("wfc2_prod p2_$g $g 1/3 $((10*IT)) 4000 1 2 9100000 1" "wfc2 d2_$g $g 1/3 $((10*IT)) 4000 1 2 9100000 1")
  jobs+=("wfc2_e3 e3f2_$g $g 1/3 $((10*IT)) 4000 1 1 35500000 1" "wfc2_e5 e5f2_$g $g 1/3 $((10*IT)) 4000 1 1 36500000 1")
done
export TS_NT=200 TS_NI=10 TS_NB=40
printf '%s\n' "${jobs[@]}" | while read -r bin tag rest; do
  echo "$bin $tag $rest"; done > v/jobs.txt
# at most 2 at a time
xargs -P "${JOBS:-2}" -L 1 bash -c 'bin=$0; tag=$1; shift 1; [ "$bin" = wfc_ts ] && export TS_DUMP=v/T_$tag.bin; nice -n 10 ./$bin run "$@" v/S_$tag.bin > v/r_$tag.txt 2> v/e_$tag.txt; echo $? > v/x_$tag' < v/jobs.txt
# ---- neutral mode (D3)
for g in 30 47; do
  nice -n 10 ./wfc_prod neutral $g 1/3 300 > v/n_p_$g.txt; nice -n 10 ./wfc neutral $g 1/3 300 > v/n_d_$g.txt
  [ "$(stripn v/n_p_$g.txt)" = "$(stripn v/n_d_$g.txt)" ] && ok "D3 neutral game $g" || bad "D3 neutral game $g"
done
for x in v/x_*; do [ "$(cat $x)" = 0 ] || bad "exit status ${x#v/x_} = $(cat $x): $(cat v/e_${x#v/x_}.txt)"; done
for g in 30 47 58 165; do
  a=$(strip v/r_p_$g.txt); [ $(echo "$a" | grep -c REP) = 2 ] || bad "D1 production gave no REP lines at game $g"
  [ "$a" = "$(strip v/r_d_$g.txt)" ] && ok "D1 REP game $g" || bad "D1 REP game $g"
  cmp -s v/S_p_$g.bin v/S_d_$g.bin && ok "D1 dump game $g" || bad "D1 dump game $g"
  [ "$a" = "$(strip v/r_t_$g.txt)" ] && ok "T1 REP game $g (TSERIES on)" || bad "T1 REP game $g (TSERIES on)"
  cmp -s v/S_p_$g.bin v/S_t_$g.bin && ok "T1 dump game $g (TSERIES on)" || bad "T1 dump game $g"
  nts=$(grep -c '^TS ' v/r_t_$g.txt); ntb=$(grep -c '^TB ' v/r_t_$g.txt)
  [ "$nts" = 400 ] && [ "$ntb" = 80 ] && ok "T1 TS/TB lines game $g ($nts, $ntb)" || bad "T1 TS/TB line counts game $g ($nts, $ntb)"
  # the TB payoffs over the second half reproduce the REP line's pay (blocks are 1/40 of the run: 20 of them)
  python3 - "v/r_t_$g.txt" <<'EOF' && ok "T1 TB second-half payoff = REP pay, TS counts = TS_NI, W-111 blocks = REP ex111, game $g" || bad "T1 TB/TS consistency game $g"
import sys
rep = {}; tb = {}; ts = {}
for l in open(sys.argv[1]):
    f = l.split()
    if not f: continue
    if f[0] == "REP": rep[int(f[9])] = f
    elif f[0] == "TB": tb.setdefault(int(f[1]), []).append(f)
    elif f[0] == "TS": ts.setdefault(int(f[1]), []).append(f)
ok = True
for r, f in rep.items():
    b = tb[r]; nb = len(b); half = [x for x in b[nb // 2:]]
    pay = sum(float(x[5]) * float(x[8]) for x in half) / sum(float(x[8]) for x in half)
    ok &= abs(pay - float(f[6])) < 1e-9 * max(1.0, abs(pay))
    ok &= all(sum(int(v) for v in x[5:13]) == int(x[4]) for x in ts[r])
    if int(f[26]):                        # wedge W: exact per-block counts; the second half must equal ex111
        ok &= sum(int(x[7]) for x in half) == int(f[27])
    else:
        ok &= all(int(x[7]) == -1 for x in b)
sys.exit(0 if ok else 1)
EOF
done
for g in 30 172; do
  a=$(strip v/r_p2_$g.txt); [ $(echo "$a" | grep -c REP) = 2 ] || bad "D2 production F2 gave no REP lines at game $g"
  [ "$a" = "$(strip v/r_d2_$g.txt)" ] && ok "D2 F2 REP game $g" || bad "D2 F2 REP game $g"
  cmp -s v/S_p2_$g.bin v/S_d2_$g.bin && ok "D2 F2 dump game $g" || bad "D2 F2 dump game $g"
done
for t in a1_30 a1_47 a2_30 a2_47 k0_30 k0_47 b0_30 b0_47 b0k0_30 b0k0_47 e3_30 e3_47 e5_30 e5_47 e3f2_30 e3f2_172 e5f2_30 e5f2_172; do
  [ "$(grep -c '^REP' v/r_$t.txt)" = 1 ] && ok "A arm $t runs: $(grep '^#' v/r_$t.txt | head -1 | cut -c1-90)" || bad "A arm $t"
done
for t in k0_30 k0_47 b0k0_30 b0k0_47; do grep -q '^MUT' v/r_$t.txt && ok "A $t MUT: $(grep '^MUT' v/r_$t.txt)" || bad "A $t no MUT line"; done
grep -q 'verdicts identical' v/e_k0v_467.txt && [ "$(cat v/x_k0v_467)" = 0 ] && ok "K0 W-111 test = production's, $(grep -o '[0-9]* W-111 verdicts' v/e_k0v_467.txt) at game 467" || bad "K0VERIFY at game 467: $(cat v/e_k0v_467.txt)"
[ $fail = 0 ] && echo "ALL CHECKS PASSED (itend $IT)" || echo "SOME CHECKS FAILED"
exit $fail
