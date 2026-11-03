#!/bin/bash
# Build the stratcert kit: libcert.so (the certificate engine, a verbatim copy of PRFull/check/cert) and, on
# Cannon, the Python dependencies numpy + gmpy2 in ./pylib (user-level pip install --target; no Cannon python
# module ships gmpy2: tested 2026-09-27 on python/3.10.9, 3.12.x and intelpython/3.9.16).
set -e
cd "$(dirname "$(readlink -f "$0")")"
if command -v sha256sum >/dev/null; then sha256sum -c --quiet SHA256SUMS; else shasum -a 256 -c --quiet SHA256SUMS; fi
echo "library copy verified against SHA256SUMS"
(cd cert && ./compile.sh)
if [ "$(uname)" != "Darwin" ]; then
  PY=${PY:-/n/sw/Mambaforge-22.11.1-4/bin/python3}      # = module python/3.10.9-fasrc01 (Mambaforge 22.11.1)
  if ! PYTHONPATH=$PWD/pylib $PY -c "import numpy, gmpy2" 2>/dev/null; then
    $PY -m pip install --no-cache-dir --target=$PWD/pylib "numpy<2" gmpy2
  fi
  PYTHONPATH=$PWD/pylib $PY -c "import sys, numpy, gmpy2; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'gmpy2', gmpy2.version())"
fi
