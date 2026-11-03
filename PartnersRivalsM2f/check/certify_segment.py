#!/usr/bin/env python3
"""The minimal binary-invisible family: an OPEN one-parameter segment of friendly rivals in
wedge N whose two binary endpoints are not friendly rivals at all.

    p(t) = (0,0,0,t, 0,0,0,0, 0,1,0,1, 1,0,1,0),   t in (0,1)

Everything below is exact rational arithmetic (../check/exact.py), no floating point."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from multiprocessing import Pool
import exact

def P(t): return [F(0),F(0),F(0),F(t), F(0),F(0),F(0),F(0),
                  F(0),F(1),F(0),F(1), F(1),F(0),F(1),F(0)]

def job(a):
    lo, hi, t = a; worst=None; warg=None
    p = P(t)
    for g in range(lo, hi):
        W = exact.freqs(p, exact.g2p(g)); d = W[2]-W[1]
        if worst is None or d < worst: worst, warg = d, g
    return worst, warg

if __name__ == "__main__":
    print("endpoints (binary):")
    for t in (F(0), F(1)):
        W = exact.freqs(P(t), P(t))
        print("   t = %s  self-play limit %s   perfect alternation: %s"
              % (t, tuple(str(x) for x in W), W[1]+W[2] == 1))
    print("\ninterior:")
    for t in (F(1,4), F(1,2), F(3,4), F(1,10), F(9,10)):
        W = exact.freqs(P(t), P(t))
        print("   t = %-4s self-play limit %s   perfect alternation: %s"
              % (t, tuple(str(x) for x in W), W[1]+W[2] == 1))
    t = F(1,2)
    worst=None; warg=None
    with Pool(12) as pool:
        for w,a in pool.imap_unordered(job, [(i, min(i+2048,65536), t) for i in range(0,65536,2048)]):
            if worst is None or w < worst: worst, warg = w, a
    print("\nat t = 1/2: min over all 65536 binary co-players of W_DC - W_CD = %s (at genome %d)"
          % (worst, warg))
    print("FRIENDLY RIVAL in wedge N at t = 1/2:", worst >= 0)
