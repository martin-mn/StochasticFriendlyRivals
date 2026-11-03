"""Exact eps->0 limit stationary distributions for the memory-one pair chain.

Polynomials in eps with Fraction coefficients (list, index = power of eps).
w is obtained by the Markov-chain-tree theorem (no division, no threshold), and the
eps -> 0 limit is the ratio of the lowest-order coefficients.  Conventions: ../m1lib.py.
"""
from fractions import Fraction as F
from itertools import product

SWAP = (0, 2, 1, 3)
DET = [tuple(F(x) for x in b) for b in product([0, 1], repeat=4)]   # genome LSB first


def pmul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return r


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else F(0)) + (b[i] if i < len(b) else F(0)) for i in range(n)]


def psub(a, b):
    return padd(a, [-x for x in b])


def clamp_poly(x):
    x = F(x)
    return [x, 1 - 2 * x]          # x + eps(1-2x)


def transition(p, q):
    """M[s][t], polynomials in eps.  p focal, q co-player; state = (focal, co)."""
    pc = [clamp_poly(x) for x in p]
    qc = [clamp_poly(x) for x in q]
    M = [[None] * 4 for _ in range(4)]
    for s in range(4):
        a = pc[s]
        b = qc[SWAP[s]]
        na = psub([F(1)], a)
        nb = psub([F(1)], b)
        M[s] = [pmul(a, b), pmul(a, nb), pmul(na, b), pmul(na, nb)]
    return M


_TREES = None


def _trees():
    """for each root i, the list of acyclic parent-maps on the other three nodes."""
    global _TREES
    if _TREES is None:
        _TREES = []
        for i in range(4):
            others = [j for j in range(4) if j != i]
            good = []
            for par in product(range(4), repeat=3):
                pm = dict(zip(others, par))
                ok = True
                for j in others:              # follow parents; must reach i
                    seen, cur = set(), j
                    while cur != i:
                        if cur in seen:
                            ok = False
                            break
                        seen.add(cur)
                        cur = pm[cur]
                    if not ok:
                        break
                if ok:
                    good.append(pm)
            _TREES.append(good)
    return _TREES


def limit_W(p, q):
    """lim_{eps->0} stationary distribution of the pair chain (p focal vs q), exact."""
    M = transition(p, q)
    T = []
    for i in range(4):
        tot = [F(0)]
        for pm in _trees()[i]:
            w = [F(1)]
            for j, k in pm.items():
                w = pmul(w, M[j][k])
            tot = padd(tot, w)
        T.append(tot)
    d = min(next(k for k, c in enumerate(t) if c != 0) for t in T)
    lead = [t[d] if d < len(t) else F(0) for t in T]
    s = sum(lead)
    return [x / s for x in lead]


def nash_rows(p):
    """(a,b,c) per deterministic deviator: p is limit-Nash at (u,v) iff a(u-1)+b v <= c
    for all 16 rows, with a=d_CD, b=d_DC, c=d_DD, d_X = W_X(sigma,p) - W_X(p,p)."""
    Wp = limit_W(p, p)
    rows = []
    for s in DET:
        Ws = limit_W(s, p)
        rows.append((Ws[1] - Wp[1], Ws[2] - Wp[2], Ws[3] - Wp[3]))
    return rows, Wp


def is_nash(p, u, v):
    rows, _ = nash_rows(p)
    u, v = F(u), F(v)
    return all(a * (u - 1) + b * v <= c for a, b, c in rows)
