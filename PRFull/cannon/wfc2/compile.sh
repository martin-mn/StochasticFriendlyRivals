#!/bin/bash
# compile.sh -- every binary of the wfc2 kit from the one source wfc.c (+ the production engine m2f.c).
#   wfc       F1 (N = 1000, beta = 100, mu = 1e-2, eps = 1e-4)          production; arms A1, A2 (alpha is runtime)
#   wfc2      F2 (N = 100,  beta = 3,   mu = 1e-4, eps = 1e-4)          production
#   wfc_k0    F1 with the local kernel K0                        -DLOCALMUT          arm K0
#   wfc_b0    F1 at beta = 0 (neutral), global kernel            -DBETAV=0.0         arm B0
#   wfc_b0k0  F1 at beta = 0 (neutral), local kernel K0          -DBETAV=0.0 -DLOCALMUT
#   wfc_e3    F1 at eps = 1e-3                                   -DEPSV=1e-3         arm E3
#   wfc_e5    F1 at eps = 1e-5                                   -DEPSV=1e-5         arm E5
#   wfc2_e3   F2 at eps = 1e-3                                   F2 flags -DEPSV=1e-3  arm E3F2
#   wfc2_e5   F2 at eps = 1e-5                                   F2 flags -DEPSV=1e-5  arm E5F2
#   wfc_ts    F1 with the instrumented time series               -DTSERIES           arm TS
# Before building, it checks that with no switch the preprocessed program (comments and whitespace aside) is token for
# token ref/wfc_prod.c, the production source (PartnersRivalsM2f/cannon/wfc/wfc.c, sha256 below), at F1 and at F2:
# the same compiler then emits the same code, so the default builds reproduce production bit for bit.
set -euo pipefail
cd "$(dirname "$0")"
CC=${CC:-cc}
CFLAGS=${CFLAGS:--O3 -march=native}
F2="-DNPOP=100 -DBETAV=3.0 -DMUV=1e-4"
PROD_SHA=06d09ef9ab6e420b2621a61df4b5ee40bc467e57d11d2ebdaa0684d129e6eb8f
sha() { if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
[ "$(sha ref/wfc_prod.c)" = "$PROD_SHA" ] || { echo "ref/wfc_prod.c is not the production source" >&2; exit 1; }
tok() { $CC -E -P -I. "$@" | tr -s ' \t\n' ' '; }
# (and the switches' own defaults spelled out, -DEPSV=1e-4 and -DBETAV=100.0, give the production program too: the
#  EPSV switch touches EPS and nothing else)
for flags in "" "$F2" "-DEPSV=1e-4" "$F2 -DEPSV=1e-4" "-DBETAV=100.0 -DEPSV=1e-4"; do
  a=$(tok $flags wfc.c | cksum); b=$(tok $flags ref/wfc_prod.c | cksum)
  [ "$a" = "$b" ] || { echo "build '${flags:-F1}' is NOT token-identical to production" >&2; exit 1; }
done
echo "default F1 and F2 builds (and -DEPSV=1e-4, -DBETAV=100.0) are token-identical to the production source"
mkdir -p out
build() { local out=$1; shift; $CC $CFLAGS "$@" -o "$out" wfc.c m2f.c -lm; echo "built $out  $*"; }
build wfc
build wfc2     $F2
build wfc_k0   -DLOCALMUT
build wfc_b0   -DBETAV=0.0
build wfc_b0k0 -DBETAV=0.0 -DLOCALMUT
build wfc_e3   -DEPSV=1e-3
build wfc_e5   -DEPSV=1e-5
build wfc2_e3  $F2 -DEPSV=1e-3
build wfc2_e5  $F2 -DEPSV=1e-5
build wfc_ts   -DTSERIES
