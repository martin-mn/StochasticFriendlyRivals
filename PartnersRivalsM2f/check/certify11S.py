#!/usr/bin/env python3
"""Exact rational certificate for wedge S (u+v<1, T<S), dimension 11.  The point
      p = (1, 1/2, 1/2, 1/2, 1, 1/2, 1/2, 1/2, 1, 1/2, 1/2, 1, 0, 1/2, 1/2, 1/2)
lies in the support stratum 1***1***1**10***, is EFFICIENT below the switch line (self-play
W_CC = 1 exactly) and is a RIVAL for T < S: W_CD - W_DC >= 0 against every one of the 65536 binary
memory-two co-players, in exact rational arithmetic with no floating point."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from multiprocessing import Pool
import exact
H = F(1,2)
P = [F(1), H, H, H, F(1), H, H, H, F(1), H, H, F(1), F(0), H, H, H]
def job(a):
    lo, hi = a; worst = None; warg = None; nzero = 0
    for t in range(lo, hi):
        W = exact.freqs(P, exact.g2p(t)); d = W[1] - W[2]
        if d == 0: nzero += 1
        if worst is None or d < worst: worst, warg = d, t
    return worst, warg, nzero
if __name__ == "__main__":
    W = exact.freqs(P, P)
    print("self-play limit (exact):  W_CC=%s W_CD=%s W_DC=%s W_DD=%s" % tuple(str(x) for x in W))
    print("efficient below the switch line (W_CC = 1):", W[0] == 1)
    worst = None; warg = None; nzero = 0
    with Pool(12) as pool:
        for w, a, z in pool.imap_unordered(job, [(i, min(i+2048,65536)) for i in range(0,65536,2048)]):
            nzero += z
            if worst is None or w < worst: worst, warg = w, a
    print("min over all 65536 binary co-players of W_CD - W_DC = %s  (at genome %d)" % (worst, warg))
    print("co-players attaining equality:", nzero)
    print("FRIENDLY RIVAL in wedge S:", (W[0] == 1) and (worst >= 0))
