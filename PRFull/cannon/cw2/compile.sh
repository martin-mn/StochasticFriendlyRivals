#!/bin/bash
# PRFull cannon/cw2: builds of cwx.f (= PartnersRivalsM1f cw.f + the switches [cw2 1-3]) and kappa.c.
#
#   bash compile.sh              the Cannon run set: cw cwe5 cwts kappa kappae5
#   bash compile.sh probe        the probe set:      cwp cwpe5 cwpts cwts1
#   bash compile.sh test         the local test set: refcwt cwt cwte5 cwtts   (ITS generations, default 20000)
#   bash compile.sh <target...>  any subset by name
#
# targets (every Fortran target is cwx.f + mt.f mtb.f payf2.f payf3.f + m1cls.o, the same flags):
#   cw       cw.x      cwx.f, no switch: THE PRODUCTION PROGRAM (cw.f's itend 1e7).  Arm E3 runs it at
#                      isl = 3*ipt+2 (eps slot 2 = 1e-3, no code change); the validation task VAL runs it
#                      at isl = 3*ipt+3, replicate 1, production seed, against the archived run1 line.
#   cwe5     cw_e5.x   -DEPS1=1.0d-5: eps slot 1 = 1e-5 (default 1e-2).  Arm E5, isl = 3*ipt+1.
#   cwts     cw_ts.x   -DTSERIES -DITEND=40000000: arm TS, 4x production length, 2000 blocks of 2e4.
#   kappa    kappa.x   kappa.c default (= ref/kappa.c): kernel masses at eps slots 1e-2/1e-3/1e-4.
#   kappae5  kappa_e5.x kappa.c -DEPS1=1e-5.
#   cwp      cwp.x     cw.x    at -DITEND=$ITP (default 1e6)      probe: E3 short rate
#   cwpe5    cwp_e5.x  cw_e5.x at -DITEND=$ITP                   probe: E5 short rate
#   cwpts    cwp_ts.x  -DTSERIES -DITEND=$ITPTS (default 4e6)    probe: TS short rate
#   cwts1    cw_ts1.x  -DTSERIES at production length 1e7        probe: TS validation vs the archive
#   refcwt   refcw_t.x ref/cw.f (the byte copy of the production source) at -DITEND=$ITS   (local test)
#   cwt      cw_t.x    cwx.f at -DITEND=$ITS                                                (local test)
#   cwte5    cw_te5.x  cwx.f -DEPS1=1.0d-5 -DITEND=$ITS                                     (local test)
#   cwtts    cw_tts.x  cwx.f -DTSERIES -DITEND=$ITS  (ITS must be a multiple of 2*NBLK=4000) (local test)
#   refcw    refcw.x   ref/cw.f, production itend                                           (local test)
#
# compilers: exactly PartnersRivalsM1f/cannon/cw/compile.sh's
#   Cannon:  module load intel/25.2.1-fasrc01; ifx -O3 -xHost -fp-model precise -fpp
#                                               icx -O2 -xHost -fp-model precise -ffp-contract=off
#            -xHost on a login node (Sapphire Rapids): run the binaries on sapphire or test ONLY.
#   local:   gfortran -O2 -ffp-contract=off -fno-fast-math -cpp;  gcc -O2 -ffp-contract=off
# Never build while an array that uses the binary is live: the compiler rewrites it in place.
set -e
cd "$(dirname "$0")"
ITS=${ITS:-20000}
ITP=${ITP:-1000000}
ITPTS=${ITPTS:-4000000}

if ! type module >/dev/null 2>&1; then
  for f in /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh; do
    [ -f $f ] && . $f && break
  done
fi
type module >/dev/null 2>&1 && module load intel/25.2.1-fasrc01
if command -v ifx >/dev/null 2>&1; then
  FC=ifx;  FFLAGS="-O3 -xHost -fp-model precise"; FPP=-fpp
  CC=icx;  CFLAGS="-O2 -xHost -fp-model precise -ffp-contract=off"
else
  FC=gfortran; FFLAGS="-O2 -ffp-contract=off -fno-fast-math"; FPP=-cpp
  CC=gcc;      CFLAGS="-O2 -ffp-contract=off"
fi
LIBF="payf3.f payf2.f mtb.f mt.f"

T="$*"
[ -n "$T" ] || T="cw cwe5 cwts kappa kappae5"
[ "$T" = "probe" ] && T="cwp cwpe5 cwpts cwts1"
[ "$T" = "test" ] && T="refcwt cwt cwte5 cwtts"

# the copied sources must still be the production files (ref/SHA256SUMS; cwx.f and kappa.c are the edited ones)
for f in mt.f mtb.f payf2.f payf3.f m1cls.c m1cls.h games.dat; do
  cmp -s "$f" "ref/$f" || { echo "REFUSED: $f differs from ref/$f"; exit 1; }
done
if command -v sha256sum >/dev/null 2>&1; then (cd ref && sha256sum -c --quiet SHA256SUMS)
else (cd ref && shasum -a 256 -c --quiet SHA256SUMS); fi || { echo "REFUSED: ref/ does not match SHA256SUMS"; exit 1; }

mkdir -p build
echo "FC = $FC $FFLAGS $FPP"
echo "CC = $CC $CFLAGS"
$CC $CFLAGS -c m1cls.c -o build/m1cls.o

fbuild () {   # fbuild <exe> <source> [defines...]
  local exe=$1 src=$2; shift 2
  (cd build && $FC $FFLAGS $FPP "$@" -o "../$exe" "../$src" \
       $(for f in $LIBF; do echo ../$f; done) m1cls.o)
  echo "built $exe  ($src $*)"
}

for t in $T; do
  case $t in
    cw)      fbuild cw.x      cwx.f ;;
    cwe5)    fbuild cw_e5.x   cwx.f -DEPS1=1.0d-5 ;;
    cwts)    fbuild cw_ts.x   cwx.f -DTSERIES -DITEND=40000000 ;;
    cwp)     fbuild cwp.x     cwx.f -DITEND=$ITP ;;
    cwpe5)   fbuild cwp_e5.x  cwx.f -DEPS1=1.0d-5 -DITEND=$ITP ;;
    cwpts)   fbuild cwp_ts.x  cwx.f -DTSERIES -DITEND=$ITPTS ;;
    cwts1)   fbuild cw_ts1.x  cwx.f -DTSERIES ;;
    refcwt)  fbuild refcw_t.x ref/cw.f -DITEND=$ITS ;;
    refcw)   fbuild refcw.x   ref/cw.f ;;
    cwt)     fbuild cw_t.x    cwx.f -DITEND=$ITS ;;
    cwte5)   fbuild cw_te5.x  cwx.f -DEPS1=1.0d-5 -DITEND=$ITS ;;
    cwtts)   fbuild cw_tts.x  cwx.f -DTSERIES -DITEND=$ITS ;;
    kappa)   $CC $CFLAGS -o kappa.x kappa.c build/m1cls.o -lm; echo "built kappa.x" ;;
    kappae5) $CC $CFLAGS -DEPS1=1e-5 -o kappa_e5.x kappa.c build/m1cls.o -lm; echo "built kappa_e5.x" ;;
    refkappa) $CC $CFLAGS -I. -o refkappa.x ref/kappa.c build/m1cls.o -lm; echo "built refkappa.x" ;;
    *) echo "unknown target $t"; exit 1 ;;
  esac
done
