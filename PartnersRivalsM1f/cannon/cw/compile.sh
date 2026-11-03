#!/bin/bash
# cw: the WF process over the CONTINUOUS memory-one strategy space (cw.f),
# its BINMODE audit build, and the reference sf1 for the audit.
#
#   bash compile.sh                      cw.x cwb.x kappa.x   (the Cannon set)
#   bash compile.sh probe                cwp.x                (the rate-probe build)
#   bash compile.sh cw cwb cwp sf1t cwt cwbt cwc cwr kappa   any subset, by name
#   bash compile.sh test                 = sf1t cwt cwbt  (reduced itend)
#
#   run F2 (N = 100, beta = 3, u = 1e-4, itend = 1e8; cw2.f = mkcw2.py(cw.f)):
#   bash compile.sh f2                   cw2.x cwb2.x         (the F2 Cannon set; kappa.x is shared)
#   bash compile.sh probe2               cwp2.x               (the F2 rate-probe build, itend 1e7)
#   bash compile.sh test2                = sf2t cwt2 cwbt2 (reduced itend)
#   NEVER plain `bash compile.sh` (or cw/cwb/cwp) on Cannon while the F1 array 48366968 is live:
#   the compiler rewrites cw.x in place under the running tasks.  The F2 targets write only *2.x.
#
# targets
#   cw.x    continuum production build           cw.f        + classifier
#   cwb.x   -DBINMODE: sf1's 16-corner process through the lineage/cache
#           machinery (audit build; bit-identical to sf1 at the same itend)
#   cwt.x   cw.x  with -DITEND=$ITS                (test, default ITS=20000)
#   cwbt.x  cwb.x with -DITEND=$ITS                (test)
#   sf1t.x  ref_sf1.f with ONLY its itend line changed to $ITS (test; the
#           copied source is build/sf1t.f and the diff is asserted to be
#           exactly that one line)
#   cwp.x   cw.x  with -DITEND=$ITP  (the Cannon rate PROBE, default
#           ITP=1000000 = itend/10; its own binary and its own directory
#           probe/, so a probe cell can never pass for a production cell)
#   cwc.x   cw.x  with -DITEND=2: runs in milliseconds; its h-file's
#           '# corner atoms' line is the 16-corner classification made by
#           the LINKED classifier through cw.f's own call (the V1 re-check)
#   cwr.x   cw.x  with -DITEND=$ITS -DCLSREP=$REP (default 11): timing build,
#           every m1cls call repeated REP times; (t_cwr - t_cwt)/(REP-1) is
#           the classifier's share of the wall time
#   kappa.x kernel-measure program kappa.c
#   --- run F2: the same builds of cw2.f (generated from cw.f by mkcw2.py; every F2 target first
#       runs `python3 mkcw2.py --check` and refuses a cw2.f that is not exactly its output) ---
#   cw2.x   continuum F2 production build        cw2.f       + classifier
#   cwb2.x  cw2.f -DBINMODE (bit-identical to sf2 at the same itend)
#   cwp2.x  cw2.x with -DITEND=$ITP2 (the F2 rate PROBE, default ITP2=10000000 = itend/10;
#           its own directory probe2/)
#   cwt2.x  cw2.x  with -DITEND=$ITS              (test)
#   cwbt2.x cwb2.x with -DITEND=$ITS              (test)
#   sf2t.x  ref_sf2.f (byte copy of DiskM1WF/Cannon/mw/sf2.f) with ONLY its itend line changed to
#           $ITS (build/sf2t.f; the one-line diff is asserted)
#
# classifier: m1cls.c (the real one).  CLS=<file.c> overrides it for tests
# only; the old test stub is in quarantine/m1cls_stub.c, and a build that
# names a stub is REFUSED for the production targets cw, cwb, cwp, kappa.
#
# compilers
#   Cannon (ifx found after module load):  ifx -O3 -xHost -fp-model precise
#                                          icx -O2 -xHost -fp-model precise -ffp-contract=off
#   -xHost on the login node (Gold 6426Y, Sapphire Rapids) targets SPR: run the binaries on
#   sapphire or test ONLY, never on shared (Cascade Lake).
#   local  (no ifx):   gfortran -O2 -ffp-contract=off -fno-fast-math
#                      gcc      -O2 -ffp-contract=off
#   -ffp-contract=off forbids fused multiply-add contraction (the one thing
#   that makes gfortran's results depend on the target), so sf1t.x and
#   cwbt.x built here are comparable bit for bit.  The SAME flags are used
#   for every Fortran target, which is what the BINMODE audit requires.
#
# BUILD ONLY WHAT YOU NEED WHILE AN ARRAY IS RUNNING (the mw lesson: the
# compiler writes the executable in place).
set -e
cd "$(dirname "$0")"
ITS=${ITS:-20000}

if ! type module >/dev/null 2>&1; then
  for f in /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh; do
    [ -f $f ] && . $f && break
  done
fi
type module >/dev/null 2>&1 && module load intel/25.2.1-fasrc01
if command -v ifx >/dev/null 2>&1; then
  FC=ifx;  FFLAGS="-O3 -xHost -fp-model precise"; FPP=-fpp
  CC=icx;  CFLAGS="-O2 -xHost -fp-model precise -ffp-contract=off"   # no FMA: the arithmetic V1/V2 validated
else
  FC=gfortran; FFLAGS="-O2 -ffp-contract=off -fno-fast-math"; FPP=-cpp
  CC=gcc;      CFLAGS="-O2 -ffp-contract=off"
fi
LIBF="payf3.f payf2.f mtb.f mt.f"

ITP=${ITP:-1000000}
ITP2=${ITP2:-10000000}
REP=${REP:-11}
CLS=${CLS:-m1cls.c}
[ -f "$CLS" ] || { echo "classifier $CLS not found"; exit 1; }

T="$*"
[ -n "$T" ] || T="cw cwb kappa"
[ "$T" = "test" ] && T="sf1t cwt cwbt"
[ "$T" = "probe" ] && T="cwp"
[ "$T" = "test2" ] && T="sf2t cwt2 cwbt2"
[ "$T" = "probe2" ] && T="cwp2"
[ "$T" = "f2" ] && T="cw2 cwb2"
case "$CLS" in *stub*)
  for t in $T; do case $t in cw|cwb|cwp|kappa|cw2|cwb2|cwp2)
    echo "REFUSED: target $t with the stub classifier $CLS"; exit 1 ;; esac; done ;;
esac

# F2: cw2.f must be exactly mkcw2.py's output of the present cw.f (diff cw.f cw2.f is the audit)
for t in $T; do case $t in cw2|cwb2|cwp2|cwt2|cwbt2)
  python3 -B mkcw2.py --check || { echo "REFUSED: cw2.f is not mkcw2.py's output"; exit 1; }; break ;;
esac; done

mkdir -p build
echo "FC  = $FC $FFLAGS"
echo "CC  = $CC $CFLAGS"
echo "CLS = $CLS"
case "$CLS" in *stub*) echo "*** WARNING: linking the TEST STUB classifier (atom = 0): NOT for production ***" ;; esac
$CC $CFLAGS -c "$CLS" -o build/m1cls.o

fbuild () {   # fbuild <exe> <source> [defines...]
  local exe=$1 src=$2; shift 2
  (cd build && $FC $FFLAGS $FPP "$@" -o "../$exe" "../$src" \
       $(for f in $LIBF; do echo ../$f; done) m1cls.o)
  ls -la "$exe"
}

for t in $T; do
  case $t in
    cw)    fbuild cw.x   cw.f ;;
    cwb)   fbuild cwb.x  cw.f -DBINMODE ;;
    cwt)   fbuild cwt.x  cw.f -DITEND=$ITS ;;
    cwbt)  fbuild cwbt.x cw.f -DBINMODE -DITEND=$ITS ;;
    cwp)   fbuild cwp.x  cw.f -DITEND=$ITP ;;
    cwc)   fbuild cwc.x  cw.f -DITEND=2 ;;
    cwr)   fbuild cwr.x  cw.f -DITEND=$ITS -DCLSREP=$REP ;;
    sf1t)
      sed "s|^      data itv  /10000000/\$|      data itv  /$ITS/|" ref_sf1.f > build/sf1t.f
      nd=$(diff ref_sf1.f build/sf1t.f | grep -c '^[<>]' || true)
      [ "$nd" -eq 2 ] || { echo "sf1t.f: expected exactly one changed line, got diff of $nd"; exit 1; }
      diff ref_sf1.f build/sf1t.f | grep '^[<>]'
      (cd build && $FC $FFLAGS $FPP -o ../sf1t.x sf1t.f $(for f in $LIBF; do echo ../$f; done))
      ls -la sf1t.x ;;
    cw2)   fbuild cw2.x   cw2.f ;;
    cwb2)  fbuild cwb2.x  cw2.f -DBINMODE ;;
    cwt2)  fbuild cwt2.x  cw2.f -DITEND=$ITS ;;
    cwbt2) fbuild cwbt2.x cw2.f -DBINMODE -DITEND=$ITS ;;
    cwp2)  fbuild cwp2.x  cw2.f -DITEND=$ITP2 ;;
    sf2t)
      sed "s|^      data itv  /100000000/\$|      data itv  /$ITS/|" ref_sf2.f > build/sf2t.f
      nd=$(diff ref_sf2.f build/sf2t.f | grep -c '^[<>]' || true)
      [ "$nd" -eq 2 ] || { echo "sf2t.f: expected exactly one changed line, got diff of $nd"; exit 1; }
      diff ref_sf2.f build/sf2t.f | grep '^[<>]'
      (cd build && $FC $FFLAGS $FPP -o ../sf2t.x sf2t.f $(for f in $LIBF; do echo ../$f; done))
      ls -la sf2t.x ;;
    kappa)
      $CC $CFLAGS -o kappa.x kappa.c build/m1cls.o -lm
      ls -la kappa.x ;;
    *) echo "unknown target $t"; exit 1 ;;
  esac
done
