#!/usr/bin/env python3
"""Independent cross-check of the symbolic limits in prove_strata.json (not part of the proof).

At random rational interior points of every stratum and for all 16 deviators (and self-play):
 (A) the stationary distribution at eps = 10^-40, solved by exact Fraction Gaussian
     elimination (no cofactors), must agree with the symbolic limit W to within 10^-35;
 (B) the project's own exact engine exact_partners.limitW (cofactor/order based, separate
     code) must return exactly the same W.
    python3 crosscheck.py
"""
import json
import os
import random
import sys
from fractions import Fraction as F

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import exact_partners  # noqa: E402

M = (0, 2, 1, 3)
a_s, b_s, c_s = sp.symbols('a b c')


def stationary_gauss(sig, p, eps):
    P = []
    for s in range(4):
        x = eps + (1 - 2 * eps) * sig[s]
        y = eps + (1 - 2 * eps) * p[M[s]]
        P.append([x * y, x * (1 - y), (1 - x) * y, (1 - x) * (1 - y)])
    # solve w (P - I) = 0, sum w = 1: transpose system, replace last eq by normalisation
    Am = [[P[j][i] - (1 if i == j else 0) for j in range(4)] for i in range(4)]
    rhs = [F(0)] * 4
    Am[3] = [F(1)] * 4
    rhs[3] = F(1)
    n = 4
    Mx = [row[:] + [rhs[i]] for i, row in enumerate(Am)]
    for col in range(n):
        piv = next(r for r in range(col, n) if Mx[r][col] != 0)
        Mx[col], Mx[piv] = Mx[piv], Mx[col]
        for r in range(n):
            if r != col and Mx[r][col] != 0:
                f = Mx[r][col] / Mx[col][col]
                Mx[r] = [Mx[r][k] - f * Mx[col][k] for k in range(n + 1)]
    return [Mx[i][n] / Mx[i][i] for i in range(n)]


def main():
    data = json.load(open(os.path.join(HERE, 'prove_strata.json')))
    rng = random.Random(7)
    eps = F(1, 10 ** 40)
    nA = nB = 0
    for st in data['strata']:
        code = st['stratum'][3:-1].split(',')
        for trial in range(4):
            pt = {}
            vals = []
            for name, cd in zip('abc', code):
                if cd == '*':
                    q = F(rng.randint(1, 998), 999)
                    pt[sp.Symbol(name)] = sp.Rational(q.numerator, q.denominator)
                    vals.append(q)
                else:
                    vals.append(F(int(cd)))
            p = (F(1),) + tuple(vals)
            for s, row in st['rows'].items():
                sig = tuple(F(int(ch)) for ch in s)
                Wsym = [sp.Rational(sp.sympify(w).subs(pt)) for w in row['W']]
                Wsym = [F(int(w.p), int(w.q)) for w in Wsym]
                wg = stationary_gauss(sig, p, eps)
                assert all(abs(wg[i] - Wsym[i]) < F(1, 10 ** 35) for i in range(4)), (st['stratum'], s, pt)
                nA += 1
                Wp, m, akin = exact_partners.limitW(sig, p, one=1)
                assert list(Wp) == Wsym and akin, (st['stratum'], s, pt, Wp, Wsym)
                nB += 1
    print('cross-check passed: %d Gaussian eps=1e-40 comparisons, %d exact comparisons with exact_partners.limitW'
          % (nA, nB))


if __name__ == '__main__':
    main()
