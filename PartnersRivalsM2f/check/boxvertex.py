#!/usr/bin/env python3
"""Numerical test of the BOX-VERTEX LEMMA: for a fixed p, is the worst memory-two co-player
always a BINARY one?  Take strategies certified as rivals by minimising over the 65536 binary
co-players and hunt for a CONTINUUM co-player q in [0,1]^16 that beats them."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m2f
from scipy.optimize import minimize

rng = np.random.default_rng(20260909)
# a bag of certified rivals: the 8 binary friendly rivals of W, the exact 7-dim continuum point,
# and random continuum rivals found by rejection
P7 = np.array([1,.5,1,.5,0,0,1,.5,1,.25,1,.25,.25,0,.25,0], dtype=float)
bag = [m2f.g2p(g) for g in [19079,19087,20103,20111,23175,23183,24199,24207]] + [P7]
PAT = '1*1*001*1*1**0*0'; free = [j for j,c in enumerate(PAT) if c=='*']
def build(v):
    p = np.zeros(16); k = 0
    for j,c in enumerate(PAT):
        p[j] = 0. if c=='0' else (1. if c=='1' else v[k]); k += (c=='*')
    return p
while len(bag) < 24:
    v = rng.uniform(0,1,len(free)); p = build(v)
    if m2f.selfplay(p)[0] > 1-1e-12 and m2f.rival_margin(p,+1,earlyquit=-1e-9)[0] >= -1e-9:
        bag.append(p)

worst_overall = 0.0
for i, p in enumerate(bag):
    binmin, arg = m2f.rival_margin(p, +1)
    def f(z):
        q = 1/(1+np.exp(-z))
        return m2f.freqs(p, q)[2] - m2f.freqs(p, q)[1]
    best = binmin; bq = None
    for trial in range(40):
        z0 = rng.normal(0, 3, 16)
        r = minimize(f, z0, method='Nelder-Mead',
                     options=dict(maxiter=4000, maxfev=4000, xatol=1e-10, fatol=1e-12))
        if r.fun < best: best = r.fun; bq = 1/(1+np.exp(-r.x))
    # also plain random search
    for trial in range(4000):
        q = rng.uniform(0,1,16)
        d = m2f.freqs(p,q)[2] - m2f.freqs(p,q)[1]
        if d < best: best = d; bq = q.copy()
    worst_overall = min(worst_overall, best - binmin)
    print("  p#%-2d  binary min %+.12f   best continuum co-player %+.12f   gap %+.3e %s"
          % (i, binmin, best, best - binmin, "  <-- VIOLATION" if best < binmin - 1e-9 else ""))
print("\nworst (continuum - binary) over the bag: %+.3e" % worst_overall)
print("box-vertex lemma survives:", worst_overall > -1e-9)
