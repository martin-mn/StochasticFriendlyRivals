"""Hidden exceptional sets inside open strata.

For every one of the 81 strata and each of the 17 chains that decide stability (self-play p vs p, and
the 16 binary deviators vs p) compute symbolically, with the free coordinates as symbols, the leading
eps-coefficient of the Markov-chain-tree normaliser Dtot(eps).  The eps->0 limit is a single rational
function on the open stratum exactly when that leading coefficient has no zero there.  (The corner
exceptions of the manuscript -- the segments (p_CD,p_DC) = (0,1),(1,0) on {p_CC=1}, the corners of the
cycle family -- are points where the lead of a stratum vanishes, so they sit in lower strata.)

Sign test: substitute x = t/(1+t) for every free coordinate, clear the denominators, and check that all
coefficients of the resulting polynomial in t > 0 have one sign.  If not, report the factorisation.
Run: python3 normlead.py   (12 processes)
"""
import itertools
import os
from multiprocessing import Pool

import sympy as sp

e = sp.symbols('e', positive=True)
M = [0, 2, 1, 3]
X = sp.symbols('x0:4', positive=True)
T = sp.symbols('t0:4', positive=True)


def clamp(x):
    return e + (1 - 2 * e) * x


def lead(sig, p):
    P = sp.zeros(4, 4)
    for s in range(4):
        x = clamp(sig[s]); y = clamp(p[M[s]])
        P[s, 0] = x * y; P[s, 1] = x * (1 - y); P[s, 2] = (1 - x) * y; P[s, 3] = (1 - x) * (1 - y)
    A = sp.eye(4) - P
    tot = 0
    for i in range(4):
        idx = [k for k in range(4) if k != i]
        tot += A.extract(idx, idx).det(method='berkowitz')
    pt = sp.Poly(sp.expand(tot), e)
    for k in range(pt.degree() + 1):
        c = sp.expand(pt.coeff_monomial(e ** k))
        if c != 0:
            return k, c
    return None, 0


def signcheck(c, free):
    if not free:
        return 'const', c > 0
    sub = {X[i]: T[i] / (1 + T[i]) for i in free}
    num, den = sp.fraction(sp.together(c.subs(sub)))
    poly = sp.Poly(sp.expand(num), *[T[i] for i in free])
    co = poly.coeffs()
    # den is a product of powers of (1+t), positive
    if all(x > 0 for x in co):
        return 'pos', True
    if all(x < 0 for x in co):
        return 'neg', True
    return 'mixed', False


def work(st):
    free = [i for i, ch in enumerate(st) if ch == 'o']
    p = tuple(sp.Integer(0) if ch == '0' else sp.Integer(1) if ch == '1' else X[i] for i, ch in enumerate(st))
    out = []
    chains = [('self', p)] + [(''.join(map(str, b)), tuple(sp.Integer(x) for x in b))
                              for b in itertools.product([0, 1], repeat=4)]
    for name, sig in chains:
        k, c = lead(sig, p)
        kind, ok = signcheck(c, free)
        if not ok:
            out.append((name, k, str(sp.factor(c))))
        else:
            out.append((name, k, kind))
    return ''.join(st), out


if __name__ == '__main__':
    strata = [''.join(s) for s in itertools.product('0o1', repeat=4)]
    with Pool(int(os.environ.get('CERT_MAX_THREADS') or 12)) as pool:
        res = pool.map(work, strata, chunksize=1)
    nbad = 0
    for st, out in res:
        bad = [o for o in out if o[2] not in ('pos', 'neg', 'const')]
        selfk = out[0][1]
        if bad:
            nbad += 1
            print(st, 'self-play order', selfk, 'NONCONSTANT-SIGN leads:', bad)
    print('strata with a lead of uncertain sign:', nbad, 'of', len(res))
