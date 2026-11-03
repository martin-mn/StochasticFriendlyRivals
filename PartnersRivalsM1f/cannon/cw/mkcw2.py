#!/usr/bin/env python3
"""
mkcw2.py -- build cw2.f (run F2) from cw.f (run F1) by an asserted list of substitutions, each of which must
match EXACTLY ONCE (the discipline of DiskM1WF/Opt/mksf.py, which built mw's sf2.f from the same source as sf1.f
the same way).  `diff cw.f cw2.f` is therefore the audit of F2 against F1, and this script asserts its size.

    F1  cw.f    N = 1000, m = 1, beta = 100, u = 1e-2, itend = 1e7   continuum seed 40000000, BINMODE 10000000
    F2  cw2.f   N = 100,  m = 1, beta = 3,   u = 1e-4, itend = 1e8   continuum seed 42000000, BINMODE 11000000

These are exactly the five parameter lines by which DiskM1WF/Cannon/mw/sf2.f differs from sf1.f (n, iseedb,
uval, bval, itv), with cw's two seed bases (the BINMODE one is sf2's own 11000000, so that BINMODE-F2 can be
compared with sf2 bit for bit; the continuum one is 42000000, disjoint from F1's 40M, kappa's 41M, and 8M 9M
10M 11M 20-25M 30-32M 39M) and itv's default, which in cw.f is the ITEND macro.  The remaining edits are comments
(and the one w-file header string) that name F1's parameters or the F1 class-table size 'maxc 2002' (-> 202,
since maxc = 2n+2).  NOTHING ELSE CHANGES: eps, nrep = nper = 10, games.dat, the kernel, the classifier, the
task decode, the loop, the RNG discipline, the outputs' formats.

cw.f is the source of the LIVE F1 array (Cannon 48366968) and is never written; its sha256 is asserted first.

usage:  python3 mkcw2.py           writes cw2.f next to cw.f and prints the non-comment diff lines
        python3 mkcw2.py --diff    also prints the full diff
        python3 mkcw2.py --check   writes nothing; exits 1 unless cw2.f is exactly this script's output (compile.sh)
"""
import difflib
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "cw.f")
DST = os.path.join(HERE, "cw2.f")
SRC_SHA = "37ec579b248ffb57716aebfa4d96f871d3852eaf5810874874baba4c729a2ebb"

# (tag, old, new): each `old` must occur exactly once in cw.f
EDITS = [
    # ---- the parameters (code) -------------------------------------------------------------------------------
    ("n", "      parameter (n=1000,m=1)\n",
          "      parameter (n=100,m=1)\n"),
    ("seed BINMODE", "      parameter (iseedb=10000000)\n",
                     "      parameter (iseedb=11000000)\n"),
    ("seed continuum", "      parameter (iseedb=40000000)\n",
                       "      parameter (iseedb=42000000)\n"),
    ("u", "      data uval /1.d-02/\n",
          "      data uval /1.d-04/\n"),
    ("beta", "      data bval /100.d+00/\n",
             "      data bval /3.d+00/\n"),
    ("itend default", "#define ITEND 10000000\n",
                      "#define ITEND 100000000\n"),
    # ---- the w-file header string that states the F1 class-table size (maxc = 2n+2) --------------------------
    ("maxc string", "     & //', nhw = class-table high water (maxc 2002); neartie = m1cls'\n",
                    "     & //', nhw = class-table high water (maxc 202); neartie = m1cls'\n"),
    # ---- comments that name F1's parameters or F1's comparison program ---------------------------------------
    ("c: BINMODE proc", "*   [cw 7] -DBINMODE: sf1's 16-corner process THROUGH the lineage/cache\n",
     "*   [cw 7] -DBINMODE: sf2's 16-corner process THROUGH the lineage/cache\n"),
    ("c: BINMODE ref", "*          machinery; bit-identical to sf1.x (the audit of 1, 2, 4)\n",
                       "*          machinery; bit-identical to sf2.x (the audit of 1, 2, 4)\n"),
    ("c: itend", "*   [cw 8] -DITEND=<n>: itend override for tests (default sf1's 1e7)\n",
                 "*   [cw 8] -DITEND=<n>: itend override for tests (default sf2's 1e8)\n"),
    ("c: seeds", "* for their history.  Continuum seeds 40000000 + (isl-1)*nper + irep;\n"
                 "* BINMODE keeps sf1's 10000000 so that it can be compared with sf1.x.\n",
                 "* for their history.  Continuum seeds 42000000 + (isl-1)*nper + irep;\n"
                 "* BINMODE keeps sf2's 11000000 so that it can be compared with sf2.x.\n"),
    ("c: BINMODE out", "* output (BINMODE): sf1.x's three files exactly; diagnostics on stdout.\n",
                       "* output (BINMODE): sf2.x's three files exactly; diagnostics on stdout.\n"),
    ("c: seed base", "C     [cw 7] seed base: sf1's in BINMODE, a disjoint one in the continuum\n",
                     "C     [cw 7] seed base: sf2's in BINMODE, a disjoint one in the continuum\n"),
    ("c: ITEND", "C     [cw 8] the literal 10000000 is now the macro ITEND (default 10000000)\n",
                 "C     [cw 8] the literal 100000000 is now the macro ITEND (default 1e8)\n"),
]
N_CODE, N_COMMENT = 7, 8          # changed lines on each side of the diff: 7 code, 8 comment (15 in all)


def is_comment(line):
    """fixed-form Fortran: a line whose column 1 is C, c, * or ! is a comment"""
    return line[:1] in ("C", "c", "*", "!")


def main():
    raw = open(SRC, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    assert sha == SRC_SHA, "cw.f sha256 %s is not the F1 production source %s" % (sha, SRC_SHA)
    s = raw.decode("ascii")
    for tag, old, new in EDITS:
        n = s.count(old)
        assert n == 1, "edit %r matched %d times, expected exactly 1" % (tag, n)
        s = s.replace(old, new)
    a, b = raw.decode("ascii").splitlines(True), s.splitlines(True)
    assert len(a) == len(b), "line count changed (%d -> %d)" % (len(a), len(b))
    changed = [(i + 1, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]
    code = [c for c in changed if not is_comment(c[1])]
    comm = [c for c in changed if is_comment(c[1])]
    assert all(is_comment(y) == is_comment(x) for _, x, y in changed), "an edit turned code into comment or back"
    assert (len(code), len(comm)) == (N_CODE, N_COMMENT), \
        "diff has %d code + %d comment lines, expected %d + %d" % (len(code), len(comm), N_CODE, N_COMMENT)
    if "--check" in sys.argv:          # compile.sh: refuse a cw2.f that is not this script's output of cw.f
        cur = open(DST).read() if os.path.exists(DST) else None
        if cur != s:
            sys.exit("cw2.f is %s -- run: python3 mkcw2.py" % ("missing" if cur is None else "NOT mkcw2.py's output"))
        print("cw2.f = mkcw2.py(cw.f): checked (%d edits, %d code + %d comment lines)"
              % (len(EDITS), len(code), len(comm)))
        return
    tmp = DST + ".new"
    with open(tmp, "w") as fh:
        fh.write(s)
    os.replace(tmp, DST)
    print("cw2.f written: %d edits, each matched exactly once; %d lines changed (%d code, %d comment), "
          "line count %d unchanged" % (len(EDITS), len(changed), len(code), len(comm), len(b)))
    print("sha256 cw.f  %s\nsha256 cw2.f %s" % (sha, hashlib.sha256(s.encode()).hexdigest()))
    print("--- diff cw.f cw2.f, non-comment lines ---")
    for ln, x, y in code:
        print("%4d < %s%4d > %s" % (ln, x, ln, y), end="")
    if "--diff" in sys.argv:
        print("--- full diff ---")
        sys.stdout.writelines(difflib.unified_diff(a, b, "cw.f", "cw2.f", n=0))


if __name__ == "__main__":
    main()
