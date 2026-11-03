"""ROUTE B (item 2) -- EFFICIENCY ABOVE THE SWITCH LINE: the ALTERNATORS.

Self-play only.  p = (p1,p2,p3,p4) in [0,1]^4, clamped p -> phat = eps + (1-2eps) p
(convention A).  The pair chain has states CC=0 CD=1 DC=2 DD=3, focal action first;
in SELF-PLAY the focal player uses phat_s and the co-player uses phat_{SWAP[s]} with
SWAP = (0,2,1,3).

A strategy is an ALTERNATOR iff its self-play limit is W = (0, 1/2, 1/2, 0).

The engine is the Markov-chain tree theorem with exact leading-order (order, coefficient)
arithmetic, exactly as ../check/strata_rivals.py does for the rival margin, but with the
co-player equal to the focal strategy instead of a vertex.  It is subtraction-free, so the
leading coefficient of every d_i is a polynomial with NON-NEGATIVE coefficients that is
STRICTLY POSITIVE on the open stratum; hence W_i = 0 iff ord(d_i) > min_j ord(d_j),
a purely combinatorial condition, CONSTANT ON EACH OPEN STRATUM.
"""
import itertools
import sys
from fractions import Fraction

import sympy as sp

sys.path.insert(0, '/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib  # noqa: E402

SWAP = (0, 2, 1, 3)
SNAME = ("CC", "CD", "DC", "DD")
CC, CD, DC, DD = 0, 1, 2, 3
X = sp.symbols('x1 x2 x3 x4', positive=True)


def arborescences(n=4):
    out = {}
    for r in range(n):
        others = [s for s in range(n) if s != r]
        trees = []
        for tgt in itertools.product(*[[t for t in range(n) if t != s] for s in others]):
            f = dict(zip(others, tgt))
            ok = True
            for s in others:
                cur, steps = s, 0
                while cur != r and steps <= n:
                    cur = f[cur]; steps += 1
                if cur != r:
                    ok = False; break
            if ok:
                trees.append([(s, f[s]) for s in others])
        out[r] = trees
    return out


ARB = arborescences()
assert all(len(v) == 16 for v in ARB.values())


def factors(slot, i):
    """((ord,coeff) of phat_i, (ord,coeff) of 1-phat_i) for slot in {0,'i',1}."""
    if slot == 0:
        return (1, sp.Integer(1)), (0, sp.Integer(1))
    if slot == 1:
        return (0, sp.Integer(1)), (1, sp.Integer(1))
    return (0, X[i]), (0, 1 - X[i])


def entries(stratum):
    """(order, leading coeff) of every self-play transition entry M[s][t]."""
    pf = [factors(stratum[i], i) for i in range(4)]
    ent = {}
    for s in range(4):
        for t in range(4):
            a = pf[s][0] if t in (CC, CD) else pf[s][1]              # focal, uses p_s
            b = pf[SWAP[s]][0] if t in (CC, DC) else pf[SWAP[s]][1]  # co-player, uses p_{SWAP s}
            ent[(s, t)] = (a[0] + b[0], sp.expand(a[1] * b[1]))
    return ent


def leading(stratum):
    ent = entries(stratum)
    res = []
    for r in range(4):
        best, acc = None, sp.Integer(0)
        for tree in ARB[r]:
            o = sum(ent[e][0] for e in tree)
            c = sp.Integer(1)
            for e in tree:
                c *= ent[e][1]
            if best is None or o < best:
                best, acc = o, sp.expand(c)
            elif o == best:
                acc = sp.expand(acc + c)
        res.append((best, sp.expand(acc)))
    return res


def limit_W(stratum):
    lead = leading(stratum)
    m = min(o for o, _ in lead)
    num = [L if o == m else sp.Integer(0) for o, L in lead]
    den = sp.expand(sum(num))
    return [sp.cancel(n / den) for n in num], lead, m


def is_alternator(stratum):
    lead = leading(stratum)
    m = min(o for o, _ in lead)
    return lead[CC][0] > m and lead[DD][0] > m, lead, m
