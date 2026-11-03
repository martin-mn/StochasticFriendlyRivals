#!/usr/bin/env bash
# pack_release.sh -- pack an assembled release tree into the archives of the Zenodo record.
#
#   bash PRFull/release/pack_release.sh TREE OUTDIR [ARCHIVE ...]
#
# TREE is a tree written by assemble_release.py --dest.  The assignment of files to archives is computed afresh by
# hosting/hosting_plan.py from TREE/CONTENTS.tsv (into OUTDIR/plan); ARCHIVE names select some of them (e.g.
# binary-paper-dumps); without names all six are packed.  Each archive holds its files under the prefix PRFull-1.0.0/,
# so that unpacking all of them in one directory rebuilds the tree.  The archives are made reproducibly and without
# anything of the packing machine: files in sorted order, owner and group root (uid and gid 0), no extended attributes,
# no macOS metadata or ._ files (COPYFILE_DISABLE=1), no ACLs, and gzip -n (no name or time in the gzip header).
# OUTDIR/ARCHIVES.sha256 gets the SHA-256 of each archive written.
set -euo pipefail
[ $# -ge 2 ] || { sed -n '2,12p' "$0"; exit 1; }
TREE="$(cd "$1" && pwd)"; OUT="$2"; shift 2
REL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX=PRFull-1.0.0
mkdir -p "$OUT/plan"; OUT="$(cd "$OUT" && pwd)"
python3 "$REL/hosting/hosting_plan.py" "$TREE" "$OUT/plan" > /dev/null
ARCS=("$@"); [ ${#ARCS[@]} -gt 0 ] || ARCS=(code binary-paper-dumps m2-archives-and-audit census evolution-runs certificates)
export COPYFILE_DISABLE=1
if tar --version 2>/dev/null | grep -q bsdtar; then
  TARF=(--no-xattrs --no-mac-metadata --no-acls --uid 0 --gid 0 --uname root --gname root -s ",^,$PREFIX/,")
else                                                                    # GNU tar
  TARF=(--no-xattrs --no-acls --owner=0 --group=0 --numeric-owner --transform "s,^,$PREFIX/,")
fi
for a in "${ARCS[@]}"; do
  name="${PREFIX}_$a.tar.gz"
  awk -F'\t' -v n="$name" 'NR > 1 && $3 == n {print $1}' "$OUT/plan/zenodo_archive_of_file.tsv" | LC_ALL=C sort > "$OUT/plan/$a.list"
  nfiles=$(wc -l < "$OUT/plan/$a.list" | tr -d ' ')
  [ "$nfiles" -gt 0 ] || { echo "no files for archive $a"; exit 1; }
  if [ "$a" = code ]; then        # CONTENTS.tsv and MANIFEST.sha256 are not listed in CONTENTS.tsv; they go with the code
    printf 'CONTENTS.tsv\nMANIFEST.sha256\n' >> "$OUT/plan/$a.list"
  fi
  (cd "$TREE" && tar -c -f - "${TARF[@]}" -T "$OUT/plan/$a.list") | gzip -n -9 > "$OUT/$name"
  echo "$name: $nfiles files, $(du -h "$OUT/$name" | cut -f1)"
done
(cd "$OUT" && shasum -a 256 "${PREFIX}"_*.tar.gz > ARCHIVES.sha256)
echo "wrote $OUT/ARCHIVES.sha256"
