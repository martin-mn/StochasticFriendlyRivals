"""Competitive range of p(x) in x = p_<DC,CD>, other coordinates fixed, against all 65536 binary co-players.
On the open stratum 0 < x < 1 the margin is m_q(x) = (a + b x)/(g + d x), numerator and normaliser affine in x
(Markov chain tree theorem: x enters each tree at most once; leading order constant on the stratum).
Reconstruct (a,b,g,d) from exact values at x = 1/100, 1/4, 2/5; verify at x = 5/13; read off the sign of the
numerator on (0, 1) from a + b*0 and a + b*1."""
from fractions import Fraction as F
from lo import code2q, pstar, margin


def load(name):
    return {int(c): F(m) for c, m in (l.split() for l in open(name))}


X = [F(1, 100), F(1, 4), F(2, 5)]
V = [load('margins_1_100.txt'), load('margins_1_4.txt'), load('margins_2_5.txt')]
V513 = load('margins_5_13.txt')


def nullvec(rows):
    """1-dim null space of a 3x4 Fraction matrix (None if not 1-dim)"""
    A = [r[:] for r in rows]
    piv = []; r = 0
    for c in range(4):
        k = next((i for i in range(r, 3) if A[i][c] != 0), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        A[r] = [v / A[r][c] for v in A[r]]
        for i in range(3):
            if i != r and A[i][c] != 0:
                A[i] = [u - A[i][c] * w for u, w in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == 3:
            break
    free = [c for c in range(4) if c not in piv]
    if len(free) != 1:
        return None
    f = free[0]
    v = [F(0)] * 4; v[f] = F(1)
    for i, c in enumerate(piv):
        v[c] = -A[i][f]
    return v


zero = const = 0
bad = []
neg_at_0 = []
roots = {}
for q in range(65536):
    ms = [Vi[q] for Vi in V]
    if all(m == 0 for m in ms):
        assert V513[q] == 0
        zero += 1
        continue
    if ms[0] == ms[1] == ms[2]:
        assert V513[q] == ms[0]
        const += 1
        if ms[0] < 0:
            bad.append(q)
        continue
    v = nullvec([[F(1), x, -m, -m * x] for x, m in zip(X, ms)])
    assert v is not None, q
    a, b, g, d = v
    if g + d * F(1, 4) < 0:
        a, b, g, d = -a, -b, -g, -d
    assert g >= 0 and g + d >= 0 and (g > 0 or g + d > 0), (q, v)     # normaliser > 0 on (0,1)
    assert (a + b * F(5, 13)) / (g + d * F(5, 13)) == V513[q], q
    if a < 0:
        neg_at_0.append(q)
    if b != 0:
        r = -a / b
        if 0 < r < 1:
            roots[r] = roots.get(r, 0) + 1
print('identically zero:', zero, '; constant non-zero:', const, '(negative:', len(bad), ')')
print('numerator negative as x -> 0+:', len(neg_at_0), neg_at_0[:10])
print('roots of numerators inside (0,1):', {str(r): n for r, n in sorted(roots.items())})
# spot check the reconstruction at a fifth point for a sample of co-players
import random
random.seed(1)
x5 = F(7, 10)
for q in random.sample(range(65536), 200) + [36228, 33153, 33920]:
    ms = [Vi[q] for Vi in V]
    if all(m == 0 for m in ms) or ms[0] == ms[1] == ms[2]:
        continue
    a, b, g, d = nullvec([[F(1), x, -m, -m * x] for x, m in zip(X, ms)])
    assert (a + b * x5) / (g + d * x5) == margin(pstar(x5), code2q(q)), q
print('reconstruction verified at x = 7/10 on a sample of 203 co-players')
