#!/usr/bin/env python3
"""relocate.py -- point the absolute paths in the programs of this release at the place where it was unpacked.

The programs were run on the author's machine, where the projects sat side by side in one directory,
/Users/martin/Documents/CL1 (PRFull, PartnersRivalsM1f, PartnersRivalsM2f, PartnersRivals1, ...).  Many scripts name
their inputs by absolute paths under it.  The release keeps that layout, so the only change needed is the prefix.

    python3 PRFull/release/relocate.py            rewrite the prefix in every program file of the tree to the
                                                  tree's own location (the directory that holds PRFull/)
    python3 PRFull/release/relocate.py --check    list the files that carry the prefix, change nothing
    python3 PRFull/release/relocate.py --undo     restore the original prefix (the files are then byte-identical
                                                  to MANIFEST.sha256 again) and remove the record of the relocation
    python3 PRFull/release/relocate.py --map OLD=NEW [--map ...] PATH ...
                                                  rewrite given prefixes in the files under PATH only, or in PATH if
                                                  it is a program file (used by the reproduce_*.sh scripts on their
                                                  working copies)

Only program files are touched (.py .sh .sbatch .slurm .c .h .f; data, logs and documentation keep the author's
paths), and only whole path prefixes (the prefix followed by '/', a quote, white space or the end of the string).
In the author's own tree the default call does nothing.  While the tree is relocated, the state is kept in
PRFull/release/relocate.json: the original prefix, the current one, and every file changed; --undo restores them byte
for byte and then removes relocate.json, after which `shasum -a 256 -c MANIFEST.sha256` passes again and the tree
holds no file that the manifest does not list.  The release's own tools (PRFull/release/) name no absolute paths and
are never rewritten.
"""
import argparse
import json
import os
import re
import sys

ORIGINAL = "/Users/martin/Documents/CL1"
EXT = (".py", ".sh", ".sbatch", ".slurm", ".c", ".h", ".f")
SKIP_DIRS = {"__pycache__", ".git", "out", "out_13", "out_probe", "run1", "run2", "log", "logs"}
MAXSIZE = 4 << 20
HERE = os.path.dirname(os.path.abspath(__file__))
STATE = os.path.join(HERE, "relocate.json")


def tree_root():
    return os.path.dirname(os.path.dirname(HERE))               # .../<tree>/PRFull/release -> <tree>


def program_files(top):
    if os.path.isfile(top):                                       # one program file (--map on a single file)
        if top.endswith(EXT) and not os.path.islink(top):
            yield top
        return
    for d, dirs, files in os.walk(top):
        if os.path.abspath(d) == HERE:                            # PRFull/release: this script and its tools
            dirs[:] = []
            continue
        dirs[:] = sorted(x for x in dirs if x not in SKIP_DIRS and not os.path.islink(os.path.join(d, x)))
        for f in sorted(files):
            p = os.path.join(d, f)
            if f.endswith(EXT) and not os.path.islink(p) and os.path.getsize(p) <= MAXSIZE:
                yield p


def pattern(prefixes):
    alts = "|".join(re.escape(p) for p in sorted(prefixes, key=len, reverse=True))
    return re.compile("(" + alts + ")(?=[/'\"\\s]|$)", re.M)


def rewrite(files, mapping, dry):
    rx = pattern(mapping.keys())
    changed = []
    for p in files:
        try:
            s = open(p, encoding="utf-8").read()
        except (UnicodeDecodeError, OSError):
            continue
        n = len(rx.findall(s))
        if not n:
            continue
        changed.append((p, n))
        if not dry:
            t = rx.sub(lambda m: mapping[m.group(1)], s)
            st = os.stat(p)
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(t)
            os.utime(p, (st.st_atime, st.st_mtime))             # keep the time stamps of the deposit
    return changed


def load_state():
    try:
        return json.load(open(STATE))
    except (OSError, ValueError):
        return {"original": ORIGINAL, "current": ORIGINAL, "files": []}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--undo", action="store_true")
    ap.add_argument("--root", default=None, help="the directory that holds PRFull/ (default: this tree's)")
    ap.add_argument("--map", action="append", default=[])
    ap.add_argument("paths", nargs="*")
    a = ap.parse_args()
    if a.map:                                                   # working copies: explicit prefixes, no state
        mapping = dict(m.split("=", 1) for m in a.map)
        files = [f for top in (a.paths or ["."]) for f in program_files(top)]
        ch = rewrite(files, mapping, a.check)
        print("relocate: %d file(s) %s under %s" % (len(ch), "would change" if a.check else "changed", " ".join(a.paths)))
        return 0
    st = load_state()
    root = os.path.abspath(a.root or tree_root())
    target = st["original"] if a.undo else root
    cur = st["current"]
    if cur == target and not a.check:
        if a.undo and os.path.exists(STATE):                    # a record left with every file restored
            os.remove(STATE)
            print("relocate: the programs already point at %s; removed %s" % (target, os.path.relpath(STATE, root)))
            return 0
        print("relocate: the programs already point at %s; nothing to do" % target)
        return 0
    files = list(program_files(root))
    ch = rewrite(files, {cur: target}, a.check)
    for p, n in ch:
        print("  %-80s %3d" % (os.path.relpath(p, root), n))
    if a.check:
        print("relocate --check: %d program file(s) name %s" % (len(ch), cur))
        return 0
    if target == st["original"]:
        # --undo: every file is restored (rewrite() stops at the first file it cannot write, and the record stays), so
        # the tree is the deposited one again and the record goes too (MANIFEST.sha256 and CONTENTS.tsv do not list it)
        if os.path.exists(STATE):
            os.remove(STATE)
        print("relocate: %d file(s) point at %s again; removed %s" % (len(ch), target, os.path.relpath(STATE, root)))
        return 0
    rel = sorted(set(st.get("files", [])) | {os.path.relpath(p, root) for p, _ in ch})
    st.update(current=target, files=rel)
    with open(STATE, "w") as fh:
        json.dump(st, fh, indent=1)
    print("relocate: %d file(s) now point at %s (state in %s)" % (len(ch), target, os.path.relpath(STATE, root)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
