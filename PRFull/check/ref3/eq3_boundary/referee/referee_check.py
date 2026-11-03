"""Referee: independent checks of the two provers' results (not a proof by itself).
(1) P1 (cofactor route) and P2 (perturbation route) closed forms of W agree symbolically on all 13 strata.
(2) Pointwise: at rational points (random and near the stratum edges), the eps->0 limit computed
    by exact Fraction Gaussian elimination at eps1 = 1e-60, eps2 = 2e-60 and Richardson
    extrapolation equals P1's formula to < 1e-90 (limits are rationals with small denominators).
(3) End-to-end: max_sigma g_sigma <= 0  <=>  I <= 0 and II <= 0 on a (u,v) grid that includes
    points on and 1e-9 off the lines I = 0, II = 0; and l1, l2 >= 0, mu = 0 at each point.
"""
import json, itertools, random, sys
from fractions import Fraction as F
import sympy as sp

import os
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
P1 = json.load(open(os.path.join(ROOT, 'P1', 'prove_strata.json')))
P2 = json.load(open(os.path.join(ROOT, 'P2', 'p2_summary.json')))
a, b, c = sp.symbols('a b c')
SIG = [''.join(t) for t in itertools.product('01', repeat=4)]
SW = (0, 2, 1, 3)


def p1code(lbl):  # '(1,0,*,1)' -> '0*1'
    return ''.join(lbl.strip('()').split(',')[1:])


def stationary(M):
    n = 4
    A = [[(F(1) if i == j else F(0)) - M[j][i] for j in range(n)] + [F(0)] for i in range(n)]
    A[n - 1] = [F(1)] * n + [F(1)]
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col] / A[col][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [A[i][n] / A[i][i] for i in range(n)]


def pi_eps(sig, p, e):
    M = []
    for s in range(4):
        x = e + (1 - 2 * e) * sig[s]
        y = e + (1 - 2 * e) * p[SW[s]]
        M.append([x * y, x * (1 - y), (1 - x) * y, (1 - x) * (1 - y)])
    return stationary(M)


def limit(sig, p):
    e1, e2 = F(1, 10 ** 60), F(2, 10 ** 60)
    p1, p2 = pi_eps(sig, p, e1), pi_eps(sig, p, e2)
    return [(e2 * x1 - e1 * x2) / (e2 - e1) for x1, x2 in zip(p1, p2)]


rng = random.Random(7)
n_sym = n_pt = n_uv = 0
bad = []
for st in P1['strata']:
    code = p1code(st['stratum'])
    rows1, rows2 = st['rows'], P2[code]['rows']
    assert len(rows1) == 16 and len(rows2) == 16
    W1 = {s: [sp.sympify(w, locals={'a': a, 'b': b, 'c': c}) for w in rows1[s]['W']] for s in SIG}
    for s in SIG:
        W2 = [sp.sympify(w, locals={'a': a, 'b': b, 'c': c}) for w in rows2[s]['W']]
        for x, y in zip(W1[s], W2):
            n_sym += 1
            if sp.cancel(x - y) != 0:
                bad.append(('P1!=P2', code, s, x, y))
        if sp.sympify(rows1[s]['mu']) != 0 or sp.sympify(rows2[s]['mu']) != 0:
            bad.append(('mu!=0', code, s))
    free = [sym for sym, ch in zip((a, b, c), code) if ch == '*']
    # points: random interior + near the edges of the open stratum
    cand = [F(1, 1000), F(999, 1000), F(1, 2)] + [F(rng.randint(1, 96), 97) for _ in range(3)]
    pts = [dict(zip(free, vals)) for vals in itertools.product(cand, repeat=len(free))] if free else [{}]
    if len(pts) > 12:
        pts = rng.sample(pts, 12)
    for pt in pts:
        pv = []
        for sym, ch in zip((a, b, c), code):
            pv.append(F(int(ch)) if ch != '*' else pt[sym])
        A_, B_, C_ = pv
        p = (F(1), A_, B_, C_)
        # self-play
        Wpp = limit((1, A_, B_, C_), p)
        if max(abs(x - y) for x, y in zip(Wpp, [1, 0, 0, 0])) > F(1, 10 ** 90):
            bad.append(('selfplay', code, pt))
        Ws = {}
        for s in SIG:
            sig = tuple(int(ch) for ch in s)
            Wn = limit(sig, p)
            Wf = [F(str(sp.Rational(w.subs({k: sp.Rational(v.numerator, v.denominator) for k, v in pt.items()}))))
                  for w in W1[s]]
            n_pt += 1
            if max(abs(x - y) for x, y in zip(Wn, Wf)) > F(1, 10 ** 90):
                bad.append(('point', code, s, pt, [float(x) for x in Wn], [str(x) for x in Wf]))
            Ws[s] = Wf
        # end-to-end iff on a (u,v) grid, using exact limits
        uvs = [(F(i, 4), F(j, 4)) for i in range(-8, 13) for j in range(-8, 9)]
        # points on and next to I=0 (v = (1-a)/c) and II=0
        for u_ in [F(-3), F(0), F(1, 2), F(1), F(2), F(5)]:
            v0 = (1 - A_) / C_
            for d in (0, F(1, 10 ** 9), -F(1, 10 ** 9)):
                uvs.append((u_, v0 + d))
            if B_ != 0:
                v1 = -(1 - A_) * (u_ - 1) / B_
                for d in (0, F(1, 10 ** 9), -F(1, 10 ** 9)):
                    uvs.append((u_, v1 + d))
            else:
                for d in (0, F(1, 10 ** 9), -F(1, 10 ** 9)):
                    uvs.append((1 + d, u_))
        for u_, v_ in uvs:
            I = C_ * v_ - (1 - A_)
            II = (1 - A_) * (u_ - 1) + B_ * v_
            gmax = max((u_ - 1) * W[1] + v_ * W[2] - W[3] for W in Ws.values())
            n_uv += 1
            if (gmax <= 0) != (I <= 0 and II <= 0):
                bad.append(('iff', code, pt, u_, v_, gmax, I, II))
print('symbolic W comparisons P1 vs P2:', n_sym)
print('pointwise limit comparisons:', n_pt)
print('(u,v) end-to-end checks:', n_uv)
print('problems:', len(bad))
for x in bad[:20]:
    print(x)
