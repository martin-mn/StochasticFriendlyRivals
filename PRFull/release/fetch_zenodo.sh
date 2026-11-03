#!/usr/bin/env bash
# fetch_zenodo.sh -- fetch the archives of the Zenodo record that the GitHub repository does not carry, unpack them
# into the release root, and check every file they bring against MANIFEST.sha256.
#
#   bash PRFull/release/fetch_zenodo.sh            the two archives the quick verification needs beyond the repository:
#                                                  binary-paper-dumps (Figures 4 and 5, SI Table 4) and
#                                                  m2-archives-and-audit (the audit of the atoms)
#   bash PRFull/release/fetch_zenodo.sh --all      all five data archives (everything reproduce_all.sh reads)
#   bash PRFull/release/fetch_zenodo.sh NAME ...   these archives (binary-paper-dumps, m2-archives-and-audit, census,
#                                                  evolution-runs, certificates)
#
# ZENODO_RECORD   the record id of release 1.0.0 (default 23923779, https://doi.org/10.5281/zenodo.23923779)
# ZENODO_BASE     where the archives are, default https://zenodo.org/records/$ZENODO_RECORD/files; a local directory
#                 (file:///path or /path) can stand in for it, e.g. to test with archives made by pack_release.sh
# The downloads go to PRFull/release/downloads and are kept; delete them after a successful run if space is short.
set -euo pipefail
REL="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$(dirname "$REL")")"
PREFIX=PRFull-1.0.0
RECORD="${ZENODO_RECORD:-23923779}"
BASE="${ZENODO_BASE:-https://zenodo.org/records/$RECORD/files}"
MAP="$REL/hosting/zenodo_archive_of_file.tsv"
case "${1:-}" in
  -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
  --all) set -- binary-paper-dumps m2-archives-and-audit census evolution-runs certificates ;;
  "") set -- binary-paper-dumps m2-archives-and-audit ;;
esac
DL="$REL/downloads"; mkdir -p "$DL"
for a in "$@"; do
  name="${PREFIX}_$a.tar.gz"
  case "$BASE" in
    file://*) cp "${BASE#file://}/$name" "$DL/$name" ;;
    /*)       cp "$BASE/$name" "$DL/$name" ;;
    *)        curl -fL --retry 3 -o "$DL/$name" "$BASE/$name?download=1" ;;
  esac
  echo "unpacking $name into $ROOT"
  tar -xzf "$DL/$name" -C "$ROOT" --strip-components 1
  # check every file of this archive against the manifest
  awk -F'\t' -v n="$name" 'NR > 1 && $3 == n {print $1}' "$MAP" > "$DL/$a.list"
  bad=$(cd "$ROOT" && awk 'NR == FNR {want[$0] = 1; next} ($2 in want)' "$DL/$a.list" MANIFEST.sha256 \
        | shasum -a 256 -c 2>/dev/null | grep -vc ': OK$' || true)
  n=$(wc -l < "$DL/$a.list" | tr -d ' ')
  if [ "$bad" != 0 ]; then echo "  $name: $bad of $n files do not match MANIFEST.sha256"; exit 1; fi
  echo "  $name: $n files, every one matches MANIFEST.sha256"
done
