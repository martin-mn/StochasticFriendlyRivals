"""Affine dimension of the set of achievable limit stationary distributions {W(sigma,p)}_16.

Akin's lemma forces every achievable w into the plane {sum w = 1} cap {w.p* = 0}, p* =
(p1-1,p2-1,p3,p4), which is 2-dimensional.  If the 16 deterministic deviators affinely SPAN
that plane, then "all 16 deviators get the same payoff" forces the payoff vector to lie in
span{p*, 1}, i.e. the exact ZD equaliser form -- a two-parameter family in the cube."""
from fractions import Fraction as F
import random
from limlib import limit_W, DET

random.seed(5)
rnd = lambda: F(random.randint(1, 59), 60)


def affrank(Ws):
    base = Ws[0]
    rows = [[x - y for x, y in zip(W, base)] for W in Ws[1:]]
    r = 0
    cols = list(range(4))
    for c in cols:
        piv = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def akin_ok(p, Ws):
    ps = [p[0] - 1, p[1] - 1, p[2], p[3]]
    return all(sum(w * x for w, x in zip(W, ps)) == 0 for W in Ws)


strata = {'interior': lambda: (rnd(), rnd(), rnd(), rnd()),
          'p1=0': lambda: (F(0), rnd(), rnd(), rnd()),
          'p1=1': lambda: (F(1), rnd(), rnd(), rnd()),
          'p2=0': lambda: (rnd(), F(0), rnd(), rnd()),
          'p2=1': lambda: (rnd(), F(1), rnd(), rnd()),
          'p3=0': lambda: (rnd(), rnd(), F(0), rnd()),
          'p3=1': lambda: (rnd(), rnd(), F(1), rnd()),
          'p4=0': lambda: (rnd(), rnd(), rnd(), F(0)),
          'p4=1': lambda: (rnd(), rnd(), rnd(), F(1))}

for name, draw in strata.items():
    ranks, akins = set(), set()
    for _ in range(12):
        p = draw()
        Ws = [limit_W(s, p) for s in DET]
        ranks.add(affrank(Ws))
        akins.add(akin_ok(p, Ws))
    print(f'{name:9s} affine rank of the 16 limit-W  {sorted(ranks)}   Akin w.p*=0 holds {akins}')
