#!/usr/bin/env bash
# the release root: run PRFull/release/reproduce_all.sh (see README.md)
exec bash "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/PRFull/release/reproduce_all.sh" "$@"
