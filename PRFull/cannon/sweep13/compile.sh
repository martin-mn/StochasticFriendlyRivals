#!/bin/bash
# compile.sh -- build census (modes psweep, pdraw, strata) for the sweep13 kit.  On Cannon: gcc from the default
# environment (16.2 as of 2026-09); the archived alpha = 0.333 run (45831015) and its probe (45830542) were built
# with the same flags by the system gcc 8.5.0 (the ELF .comment of the archived binary).  Writes BUILD.txt.
#   CC=cc bash compile.sh      builds with another compiler (the Mac: Apple clang, for the local checks)
set -euo pipefail
cd "$(dirname "$0")"
CC=${CC:-gcc}
CFLAGS=${CFLAGS:--O3 -march=native}
$CC $CFLAGS -o census census.c m2f.c -lm
mkdir -p out out_probe
sha() { if command -v sha256sum >/dev/null; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
{ echo "built $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname)"; echo "$CC $CFLAGS -o census census.c m2f.c -lm"
  $CC --version | head -1; sha census.c m2f.c witnesses.h census; } > BUILD.txt
cat BUILD.txt
