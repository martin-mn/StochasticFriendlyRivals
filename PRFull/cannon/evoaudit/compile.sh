#!/bin/bash
# Build the two libraries of the audit.
#   lib/cert/libcert.{so,dylib} : the certificate engine (outward-rounded intervals; -O2, no FP contraction, no fast-math)
#   src/libprod.{so,dylib}      : production's classifier (src/m2f.c = PartnersRivalsM2f/cannon/wfc/m2f.c unchanged,
#                                 src/prodcls.c = wfc.c's classify with full scans), built with production's flags
#                                 on Cannon: GCC 8.5.0 (/usr/bin/cc), -O3 -march=native.
set -e
cd "$(dirname "$0")"
bash lib/cert/compile.sh
if [ "$(uname)" = "Darwin" ]; then
  cc -O3 -arch arm64 -arch x86_64 -fPIC -shared -o src/libprod.dylib src/prodcls.c src/m2f.c -lm
  echo "built $(pwd)/src/libprod.dylib (macOS: not production's compiler; the reproduction check is only indicative here)"
else
  /usr/bin/cc --version | head -1
  /usr/bin/cc -O3 -march=native -fPIC -shared -o src/libprod.so src/prodcls.c src/m2f.c -lm
  echo "built $(pwd)/src/libprod.so"
fi
