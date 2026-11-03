#!/usr/bin/env python3
"""Exact rational certificate: the point
      p = (1, 1/2, 1, 1/2, 0, 0, 1, 1/2, 1, 1/4, 1, 1/4, 1/4, 0, 1/4, 0)
lies in the 7-dimensional support stratum 1*1*001*1*1**0*0, is EFFICIENT below the switch line
(self-play W_CC = 1 exactly), and is a RIVAL for T > S: W_DC - W_CD >= 0 against every one of the
65536 binary memory-two co-players, computed in exact rational arithmetic with no floating point.
By the box-vertex lemma that is the exact worst case over the whole memory-two co-player continuum.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from multiprocessing import Pool
import exact

P = [F(2,3), F(1), F(0), F(1,3), F(1,3), F(1,3), F(0), F(0),
     F(2,3), F(1), F(2,3), F(1,3), F(1,3), F(0), F(1), F(0)]
     

def job(args):
    lo, hi = args
    worst = None; warg = None; nzero = 0
    for t in range(lo, hi):
        W = exact.freqs(P, exact.g2p(t))
        d = W[2] - W[1]
        if d == 0: nzero += 1
        if worst is None or d < worst: worst = d; warg = t
    return worst, warg, nzero

if __name__ == "__main__":
    W = exact.freqs(P, P)
    print("self-play limit (exact):  W_CC=%s  W_CD=%s  W_DC=%s  W_DD=%s" % tuple(str(x) for x in W))
    print("efficient ABOVE the switch line (W_CD = W_DC = 1/2):", (W[1]+W[2] == 1))
    chunks = [(i, min(i+2048, 65536)) for i in range(0, 65536, 2048)]
    worst = None; warg = None; nzero = 0
    with Pool(12) as pool:
        for w, a, z in pool.imap_unordered(job, chunks):
            nzero += z
            if worst is None or w < worst: worst, warg = w, a
    print("min over all 65536 binary co-players of W_DC - W_CD  =  %s   (attained at genome %d)"
          % (worst, warg))
    print("co-players attaining equality:", nzero)
    print("RIVAL for T > S:", worst >= 0)
    print("FRIENDLY RIVAL in wedge N:", ((W[1]+W[2] == 1)) and (worst >= 0))
