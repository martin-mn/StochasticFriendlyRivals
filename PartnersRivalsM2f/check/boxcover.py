#!/usr/bin/env python3
"""Minimal 'box' cover of a set of support strata -- the 3-valued analogue of
../../PartnersRivals/check/families.py.

A stratum is a 16-tuple over {0, i, 1} (p_j = 0 / interior / 1).  A BOX is a product of subsets
S_j of {0,i,1}; the six subsets have readable meanings

   {0}       p_j = 0            {1}       p_j = 1          {i}     0 < p_j < 1
   {0,i}     p_j < 1            {i,1}     p_j > 0          {0,1}   p_j binary
   {0,i,1}   p_j free

A box is PRIME for the set if all its strata lie in the set and no coordinate's subset can be
enlarged.  A minimum cover is a smallest family of primes whose union is the set.
"""
import sys, itertools, collections
from scipy.optimize import milp, LinearConstraint, Bounds
import numpy as np

SYM = {0:'0', 1:'i', 2:'1'}
# readable name for each subset of {0,1,2}
NAME = {frozenset([0]):'0', frozenset([1]):'i', frozenset([2]):'1',
        frozenset([0,1]):'<', frozenset([1,2]):'>', frozenset([0,2]):'b',
        frozenset([0,1,2]):'*'}
MEAN = {'0':'= 0','1':'= 1','i':'in (0,1)','<':'< 1','>':'> 0','b':'in {0,1}','*':'free'}

def digits(g):
    d=[]
    for _ in range(16): d.append(g%3); g//=3
    return tuple(d)

def grow(strata):
    """all maximal boxes (primes) contained in the stratum set, by coordinate-wise expansion"""
    S = set(strata)
    primes = set()
    # start from singleton boxes and greedily expand each coordinate in every order-independent way
    frontier = set(tuple(frozenset([d]) for d in s) for s in S)
    seen = set(frontier)
    while frontier:
        nxt = set()
        for box in frontier:
            grew = False
            for j in range(16):
                for extra in (0,1,2):
                    if extra in box[j]: continue
                    nb = list(box); nb[j] = frozenset(box[j] | {extra}); nb = tuple(nb)
                    if nb in seen: continue
                    if all(t in S for t in itertools.product(*nb)):
                        seen.add(nb); nxt.add(nb); grew = True
            if not grew: primes.add(box)
        frontier = nxt
    # keep only maximal ones
    out = []
    for b in primes:
        if not any(b is not c and all(b[j] <= c[j] for j in range(16)) for c in primes):
            out.append(b)
    return out

def render(box):
    return "".join(NAME[box[j]] for j in range(16))

def mincover(strata, primes):
    S = sorted(set(strata)); idx = {s:i for i,s in enumerate(S)}
    cols = []
    for b in primes:
        col = np.zeros(len(S))
        for t in itertools.product(*b): col[idx[t]] = 1
        cols.append(col)
    A = np.array(cols).T
    res = milp(c=np.ones(len(primes)), constraints=LinearConstraint(A, lb=1, ub=np.inf),
               integrality=np.ones(len(primes)), bounds=Bounds(0,1))
    if not res.success: return None
    return [primes[k] for k in range(len(primes)) if res.x[k] > 0.5]

if __name__ == "__main__":
    import glob
    tag = sys.argv[1]; pat = sys.argv[2]
    S = set()
    for f in glob.glob(pat):
        for line in open(f): S.add(int(line.split()[0]))
    D = [digits(g) for g in sorted(S)]
    print("%s: %d strata" % (tag, len(D)))
    primes = grow(D)
    print("prime boxes: %d" % len(primes))
    cov = mincover(D, primes)
    print("minimum cover: %d boxes" % len(cov))
    rows = ["CC","CD","DC","DD"]
    for b in sorted(cov, key=lambda b: -sum(len(x) for x in b)):
        r = render(b)
        n = 1
        for x in b: n *= len(x)
        print("\n  %s   (%d strata)" % (r, n))
        for k in range(4):
            print("        %s | %s" % (rows[k], "  ".join(r[4*k+m] for m in range(4))))
    print("\n  legend:", ", ".join("%s : p %s" % (k, v) for k,v in MEAN.items()))
