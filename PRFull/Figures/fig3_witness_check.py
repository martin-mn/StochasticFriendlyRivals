#!/usr/bin/env python3
"""Exact check of the one-probability witness of Figure 3c (not needed to draw the figure).

    p(t) = (0,0,0,t, 0,0,0,0, 0,1,0,1, 1,0,1,0)

For each t in a small set of rationals it computes, in exact Fractions with the eps->0 engine
../../PartnersRivalsM2f/check/exact.py (imported read-only, no bytecode written there):
  * the self-play limit (w_CC, w_CD, w_DC, w_DD);
  * min over all 65536 binary memory-two co-players of w_DC - w_CD (the competitive margin for
    T > S; by the box-vertex lemma of PartnersRivalsM2f this is the worst case over [0,1]^16).

    PYTHONDONTWRITEBYTECODE=1 nice -n 19 python3 fig3_witness_check.py  > fig3_witness_check.out
"""
import os
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivalsM2f/check")
from fractions import Fraction as F
from multiprocessing import Pool

import exact                                                    # noqa: E402

NPROC = 6


def P(t):
    return [F(0), F(0), F(0), F(t), F(0), F(0), F(0), F(0),
            F(0), F(1), F(0), F(1), F(1), F(0), F(1), F(0)]


def job(a):
    lo, hi, t = a
    p = P(t)
    worst, warg, nzero = None, None, 0
    for g in range(lo, hi):
        W = exact.freqs(p, exact.g2p(g))
        d = W[2] - W[1]
        if d == 0:
            nzero += 1
        if worst is None or d < worst:
            worst, warg = d, g
    return worst, warg, nzero


if __name__ == "__main__":
    TS = [F(0), F(1, 10), F(1, 4), F(1, 2), F(3, 4), F(9, 10), F(1)]
    chunks = [(i, min(i + 2048, 65536)) for i in range(0, 65536, 2048)]
    with Pool(NPROC) as pool:
        for t in TS:
            W = exact.freqs(P(t), P(t))
            worst, warg, nz = None, None, 0
            for w, a, z in pool.imap_unordered(job, [(lo, hi, t) for lo, hi in chunks]):
                nz += z
                if worst is None or w < worst:
                    worst, warg = w, a
            print("t = %-5s self-play %-28s min_q (w_DC - w_CD) = %-6s at q = %-6d  [%d co-players at 0]"
                  % (t, "(" + ", ".join(str(x) for x in W) + ")", worst, warg, nz), flush=True)
