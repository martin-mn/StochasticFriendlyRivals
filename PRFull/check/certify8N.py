#!/usr/bin/env python3
"""Exact rational certificate for wedge N (u+v>1, T>S), dimension 8.  The point
      p = (2/3, 1, 0, 1/3, 1/3, 1/3, 0, 0, 2/3, 1, 2/3, 1/3, 1/3, 0, 1, 0)
lies in the 8-dimensional support stratum *10***00*1***010, is EFFICIENT ABOVE the switch line
(its self-play limit is perfect alternation, (W_CC, W_CD, W_DC, W_DD) = (0, 1/2, 1/2, 0) exactly)
and is a RIVAL for T > S: W_DC - W_CD >= 0 against every one of the 65536 binary memory-two
co-players, computed in exact rational arithmetic with no floating point.  By the box-vertex lemma
that is the exact worst case over the whole memory-two co-player continuum, so p is a FRIENDLY
RIVAL of the wedge N.  The minimum margin is exactly 0, attained at 31024 co-players (SI Table 2).

Corrected copy of PartnersRivalsM2f/check/certify8N.py (which is frozen and was not edited):
  * the docstring of the original was a copy of certify7.py's (the wedge-W point
    (1,1/2,1,1/2,0,0,1,1/2,1,1/4,1,1/4,1/4,0,1/4,0), "7-dimensional stratum 1*1*001*1*1**0*0",
    "EFFICIENT below the switch line (self-play W_CC = 1 exactly)"); the code always certified the
    N point above, and this docstring now says so;
  * the efficiency test checks W_CD = W_DC = 1/2 explicitly (the original tested W_CD + W_DC = 1,
    which is equivalent because self-play is symmetric under exchanging the players);
  * the process pool is capped at 3 (Martin's local-compute rule) instead of 12.
The computation itself is unchanged: PartnersRivalsM2f/check/exact.py (read-only, imported without
writing bytecode), all 65536 co-players, exact Fractions.

Run:   cd PRFull/check && nice -n 10 python3 certify8N.py        (about 30 s on 3 processes)
"""
import sys
sys.dont_write_bytecode = True                   # nothing may be written into the read-only M2f tree
sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivalsM2f/check")
from fractions import Fraction as F
from multiprocessing import Pool
import exact                                     # PartnersRivalsM2f/check/exact.py (read-only)

NPROC = 3
P = [F(2,3), F(1), F(0), F(1,3), F(1,3), F(1,3), F(0), F(0),
     F(2,3), F(1), F(2,3), F(1,3), F(1,3), F(0), F(1), F(0)]

def job(args):
    lo, hi = args
    worst = None; warg = None; nzero = 0
    for t in range(lo, hi):
        W = exact.freqs(P, exact.g2p(t))
        d = W[2] - W[1]
        if d == 0: nzero += 1
        if worst is None or d < worst or (d == worst and t < warg): worst = d; warg = t
    return worst, warg, nzero

if __name__ == "__main__":
    W = exact.freqs(P, P)
    print("self-play limit (exact):  W_CC=%s  W_CD=%s  W_DC=%s  W_DD=%s" % tuple(str(x) for x in W))
    eff = (W[0] == 0 and W[3] == 0 and W[1] == F(1, 2) and W[2] == F(1, 2))
    print("efficient ABOVE the switch line (W_CD = W_DC = 1/2):", eff)
    chunks = [(i, min(i+2048, 65536)) for i in range(0, 65536, 2048)]
    worst = None; warg = None; nzero = 0
    with Pool(NPROC) as pool:
        for w, a, z in pool.imap_unordered(job, chunks):
            nzero += z
            if worst is None or w < worst or (w == worst and a < warg): worst, warg = w, a
    print("min over all 65536 binary co-players of W_DC - W_CD  =  %s   (attained first at genome %d)"
          % (worst, warg))
    print("co-players attaining equality:", nzero)
    print("RIVAL for T > S:", worst >= 0)
    print("FRIENDLY RIVAL in wedge N:", eff and (worst >= 0))
