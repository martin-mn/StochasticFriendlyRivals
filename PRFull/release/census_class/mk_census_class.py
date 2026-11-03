#!/usr/bin/env python3
"""mk_census_class.py -- the support-based efficiency class of every one of the 3^16 memory-two strata, in the
census's own stratum code, and joined to the archived census output (referee 2, additional comment 3).

    nice -n 10 python3 mk_census_class.py [--threads 3]        (about 30 s; writes only in this directory)

What the census archived.  The census (PartnersRivalsM2f/cannon/census.c, Cannon array 45751662) wrote one line
"g cls nP nM K" for each stratum in which at least one draw was confirmed competitive, where cls is the class read
from the weights at the stratum's first draw with the tolerance 1e-12 ('C' w_CC > 1-1e-12, 'A' w_CD+w_DC > 1-1e-12,
'Z' w_DD < 1e-12, '-' otherwise).  PartnersRivalsM2f/check/st_eff.txt keeps the 1 524 791 of these lines whose
class is not '-'.  It is therefore the tolerance-based class at one sampled point, and only for those strata.

What this adds.  The class decided by the support alone (lemma L1 of SI section "Certificates"): the exact integer
exponents of the leading-order state reduction give the support of the self-play limit, which is constant on the
open stratum, and the class follows without a tolerance -- C if the support lies in the states whose last outcome
is CC, A if it lies in CD and DC, Z if it avoids DD, none otherwise.  This is PRFull/check/ubound/ub.c, mode
"classes" (compiled here into a temporary directory, from the unchanged source), for all 43 046 721 strata.

Writes (here):
  census_support_class.u8.gz    43 046 721 bytes after gunzip, one per stratum, indexed by the census code
                                g = sum_j d_j 3^j (d_j = 0: p_j = 0, 1: p_j free in (0,1), 2: p_j = 1, as census.c);
                                0 none, 1 C, 2 A, 3 Z
  st_eff_support_class.txt.gz   st_eff.txt with a sixth column: the support-based class (-, C, A or Z)
  census_class.json             counts, the checks below, and SHA-256 of both uncompressed files
Checks (all must pass; exit status 1 otherwise): the class counts are those of the Methods (C 6 206 077,
A 2 788 386, Z 2 471 566); every C stratum has p_0 = 1; every st_eff.txt line has the same class by support
as the census read from the weights (0 mismatches).
"""
import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PRF = os.path.dirname(os.path.dirname(HERE))                        # .../PRFull
CL1 = os.path.dirname(PRF)                                          # the directory that holds the projects
UB = os.path.join(PRF, "check", "ubound")
STEFF = os.path.join(CL1, "PartnersRivalsM2f", "check", "st_eff.txt")
N3 = 3 ** 16
CODE = {0: "-", 1: "C", 2: "A", 3: "Z"}
EXPECT = {"-": 31580692, "C": 6206077, "A": 2788386, "Z": 2471566}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=3)
    ap.add_argument("--cc", default=os.environ.get("CC", "cc"))
    a = ap.parse_args()
    t0 = time.time()
    tmp = tempfile.mkdtemp(prefix="census_class_")
    try:
        exe = os.path.join(tmp, "ub")
        subprocess.check_call([a.cc, "-O2", "-ffp-contract=off", "-fno-fast-math", "-I", UB, "-o", exe,
                               os.path.join(UB, "ub.c"), "-lm", "-lpthread"])
        out = subprocess.check_output([exe, "classes", os.path.join(tmp, "cls.bin"), str(a.threads)], text=True)
        print(out.strip())
        ub = np.fromfile(os.path.join(tmp, "cls.bin"), dtype=np.uint8)
    finally:
        shutil.rmtree(tmp, True)
    assert ub.size == N3
    # ub.c indexes a stratum by sum_j u_j 3^j with u = 0 (p_j = 0), 1 (p_j = 1), 2 (free); the census by
    # sum_j d_j 3^j with d = 0 (0), 1 (free), 2 (1).  Relabel every digit: new[d] = old[u(d)], u = (0, 2, 1)[d].
    arr = ub.reshape((3,) * 16)                                    # C order: axis k is the digit j = 15 - k
    for k in range(16):
        arr = np.take(arr, [0, 2, 1], axis=k)
    cen = np.ascontiguousarray(arr).reshape(-1)
    # checks: the counts, p_0 = 1 on class C (census digit 0 of stratum g is g % 3; p_0 = 1 is d_0 = 2)
    counts = {CODE[c]: int((cen == c).sum()) for c in range(4)}
    ok_counts = counts == EXPECT
    g = np.arange(N3, dtype=np.int64)
    ok_p0 = bool(np.all(g[cen == 1] % 3 == 2))
    del g
    # the join with the archived census lines
    lines = open(STEFF).read().split("\n")
    rows = [l.split() for l in lines if l.strip()]
    gi = np.array([int(r[0]) for r in rows], dtype=np.int64)
    got = cen[gi]
    want = np.array([{"-": 0, "C": 1, "A": 2, "Z": 3}[r[1]] for r in rows], dtype=np.uint8)
    mism = int((got != want).sum())
    by = {c: int((want == v).sum()) for c, v in (("C", 1), ("A", 2), ("Z", 3))}
    # write
    raw = cen.tobytes()
    with gzip.GzipFile(os.path.join(HERE, "census_support_class.u8.gz"), "wb", mtime=0) as fh:
        fh.write(raw)
    joined = "".join("%s %s\n" % (" ".join(r), CODE[int(c)]) for r, c in zip(rows, got)).encode()
    with gzip.GzipFile(os.path.join(HERE, "st_eff_support_class.txt.gz"), "wb", mtime=0) as fh:
        fh.write(joined)
    rep = dict(
        strata=N3, counts_by_support=counts, counts_expected_from_the_Methods=EXPECT, counts_ok=ok_counts,
        class_C_has_p0_equal_1=ok_p0,
        st_eff_lines=len(rows), st_eff_lines_by_census_class=by, st_eff_class_mismatches_with_support=mism,
        source_ub=os.path.relpath(os.path.join(UB, "ub.c"), CL1), source_census=os.path.relpath(STEFF, CL1),
        sha256_st_eff_txt=hashlib.sha256("\n".join(lines).encode()).hexdigest(),
        sha256_census_support_class_u8=hashlib.sha256(raw).hexdigest(),
        sha256_st_eff_support_class_txt=hashlib.sha256(joined).hexdigest(),
        encoding=("census_support_class.u8: byte g = class of the stratum with census code g = sum_j d_j 3^j, "
                  "d_j = 0: p_j = 0, 1: p_j free, 2: p_j = 1; 0 none, 1 C (w_CC = 1), 2 A (w_CD + w_DC = 1), "
                  "3 Z (w_DD = 0 only)"),
        seconds=round(time.time() - t0, 1))
    with open(os.path.join(HERE, "census_class.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    print("support-based classes of all %d strata: %s (the Methods: %s) %s" % (N3, counts, EXPECT, "ok" if ok_counts else "DIFFERS"))
    print("class C only with p_0 = 1: %s" % ("ok" if ok_p0 else "FAILED"))
    print("st_eff.txt: %d lines %s, class by support differs from the census's class at %d" % (len(rows), by, mism))
    print("wrote census_support_class.u8.gz, st_eff_support_class.txt.gz, census_class.json in %.1f s" % (time.time() - t0))
    sys.exit(0 if (ok_counts and ok_p0 and mism == 0) else 1)


if __name__ == "__main__":
    main()
