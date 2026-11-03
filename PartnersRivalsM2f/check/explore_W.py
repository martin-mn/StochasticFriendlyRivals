#!/usr/bin/env python3
"""Wedge W: how much of the continuum does the single binary friendly-rival family open into?

W = {u+v<1, u-v<1} contains the classical Prisoner's Dilemma.  Its eight binary friendly rivals
are the family 111*000101*1*010 (position j = state j = 4*recent+before), ids
19079 19087 20103 20111 23175 23183 24199 24207.  Below the switch line efficiency is
self-play W_CC = 1; rivalry for T > S is min over co-players of (W_DC - W_CD) >= 0.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m2f

EIGHT = [19079, 19087, 20103, 20111, 23175, 23183, 24199, 24207]
TOL = -1e-9

def eff_CC(p):   return m2f.selfplay(p)[0] > 1 - 1e-12
def rivP(p):     return m2f.rival_margin(p, +1, earlyquit=TOL)[0] >= TOL
def frW(p):      return eff_CC(p) and rivP(p)

def label(p, j, t):
    return "j=%2d t=%.3f" % (j, t)

print("=== coordinate-wise relaxation of each of the eight, below the switch line ===")
grid = np.array([0.0, 1e-6, 0.01, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.95, 0.99, 1-1e-6, 1.0])
summary = {}
for g in EIGHT:
    p0 = m2f.g2p(g)
    assert frW(p0), g
    row = []
    for j in range(16):
        keep = []
        for t in grid:
            p = p0.copy(); p[j] = t
            e = eff_CC(p); r = rivP(p)
            keep.append("F" if (e and r) else ("e" if e else ("r" if r else ".")))
        row.append("".join(keep))
    summary[g] = row
    print("\n-- %d  %s" % (g, "".join(str((g >> j) & 1) for j in range(16))))
    print("   grid:", " ".join("%.3g" % t for t in grid))
    for j in range(16):
        free = all(c == "F" for c in summary[g][j])
        print("   j=%2d (recent=%s,before=%s) p0=%d  %s %s"
              % (j, "CC CD DC DD".split()[j >> 2], "CC CD DC DD".split()[j & 3],
                 int(p0[j]), summary[g][j], "  <-- FULLY FREE" if free else ""))
