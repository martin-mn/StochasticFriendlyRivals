#!/usr/bin/env python3
"""frmass.py -- the kernel mass of the friendly rivals (atom 111) at any alpha, per wedge, for the enrichments of
the alpha arms A1 (1/4) and A2 (2/5).  Local, read-only, a few seconds.

The friendly-rival set of each wedge is fixed (111 depends only on the wedge: every friendly rival is Nash), so its
alpha-mass is a sum over the support strata s of  alpha^npin(s) (1-2alpha)^nfree(s) x (friendly-rival fraction of s).
The fractions are the exhaustive stratum census's hit fractions nP/K (wedges W and N) and nM/K (S and E),
PartnersRivalsM2f/check/st_eff.txt, exactly as PRFull/notes/m2_verify/verify_m2.py section D reads it.

Since 2026-09-27 the base is the rarity sweep at alpha = 1/3 exactly (PRFull cannon/sweep13, array 48928780,
1 633 280 000 draws; FR13 below), the same rates as Figures/fig45.py and SI Data 1 use:
    mass(alpha) = FR13 x census(alpha) / census(1/3)
so the arms and production share one convention.  (Before, the base was the sweep at alpha = 0.333, FR0333, carried by
census(alpha) / census(0.333); that version is frmass_0333.py.bak / frmass_0333.json.)  The raw census values are printed alongside (they agree with
Referee-triage report 5: 1/4 -> W 3.21e-7, S 9.53e-3, N, E 1.31e-5; 2/5 -> W 1.955e-5, S 1.676e-2, N, E 3.35e-4).

    python3 frmass.py            prints the table and writes frmass.json (read by analyze.py)
"""
import json
import os
import sys

import numpy as np

ST = "/Users/martin/Documents/CL1/PartnersRivalsM2f/check/st_eff.txt"
FR0333 = {"W": 3.9111e-06, "S": 1.3426e-02, "N": 9.9431e-05, "E": 9.9162e-05}     # figwfc.py, Table 4 (alpha = 0.333)
NDRAW13 = 1633280000                                                                   # sweep13: counts at alpha = 1/3
FR13 = {"W": 6508 / NDRAW13, "S": 21949101 / NDRAW13, "N": 164451 / NDRAW13, "E": 163245 / NDRAW13}
ALPHAS = {"1o3": 1 / 3, "1o10": 1 / 10, "1o5": 1 / 5, "1o4": 1 / 4, "2o5": 2 / 5, "1o2": 1 / 2}


def census(a, d):
    g, cls, nP, nM, K, npin, nfree = d
    w = a ** npin * (1 - 2 * a) ** nfree
    return {lab: float((w[cls == c] * n[cls == c] / K[cls == c]).sum())
            for lab, c, n in (("W", "C", nP), ("S", "C", nM), ("N", "A", nP), ("E", "A", nM))}


def main():
    raw = np.loadtxt(ST, dtype=str)
    g = raw[:, 0].astype(np.int64); cls = raw[:, 1]; nP = raw[:, 2].astype(float); nM = raw[:, 3].astype(float)
    K = raw[:, 4].astype(float)
    nfree = np.zeros(len(g), int); x = g.copy()
    for _ in range(16):
        nfree += (x % 3 == 1); x //= 3
    d = (g, cls, nP, nM, K, 16 - nfree, nfree)
    base = census(1 / 3, d)
    b333 = census(0.333, d)
    out = {}
    print("alpha   wedge  census(alpha)   census/census(1/3)   mass used (FR13 x ratio)   old (FR0333 route)")
    for tag, a in ALPHAS.items():
        c = census(a, d)
        out[tag] = {w: FR13[w] * c[w] / base[w] for w in "WSNE"}
        for w in "WSNE":
            print("%-6s  %s      %.4e      %.5f             %.4e                %.4e" %
                  (tag, w, c[w], c[w] / base[w], out[tag][w], FR0333[w] * c[w] / b333[w]))
    for w in "WSNE":
        assert abs(out["1o3"][w] / FR13[w] - 1) < 1e-12
    print("alpha = 1/3: the sweep's own rates")
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "frmass.json")
    json.dump(out, open(fn, "w"), indent=1)
    print("wrote", fn)


if __name__ == "__main__":
    sys.exit(main())
