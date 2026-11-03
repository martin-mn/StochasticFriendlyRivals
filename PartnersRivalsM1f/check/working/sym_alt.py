#!/usr/bin/env python3
"""ABOVE the switch line: the exact Nash polygon of an efficient alternator.

Efficiency above u+v = 1 forces (p2,p3) in {(0,1),(1,0)} with p1 < 1, p4 > 0 and not
(p1,p4) = (0,1)  (item 2, Eff_alt).  So the efficient set is TWO two-parameter sheets,
and the whole question "which games have a partner above the line" is a question about
the pair (p1,p4) only.  This script derives, in sympy with eps symbolic and p1, p4
symbolic, the eps -> 0 limit facet triples (a,b,c) of

    a (u-1) + b v <= c,        a = W_CD(s,p) - 1/2, b = W_DC(s,p) - 1/2, c = W_DD(s,p)

for all 16 deterministic deviators s and both corners.
"""
import sympy as sp

e, p1, p4, u, v = sp.symbols('e p1 p4 u v', real=True)
SWAP = (0, 2, 1, 3)


def clamp(x):
    return e + (1 - 2 * e) * x


def stat(P, Q):
    """Exact stationary row vector of the pair chain, focal P against co-player Q."""
    P = [clamp(x) for x in P]
    Q = [clamp(x) for x in Q]
    M = sp.zeros(4, 4)
    for s in range(4):
        a, b = P[s], Q[SWAP[s]]
        M[s, 0], M[s, 1], M[s, 2], M[s, 3] = a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)
    A = (M.T - sp.eye(4)).as_mutable()
    A[3, :] = sp.Matrix([[1, 1, 1, 1]])
    w = A.solve(sp.Matrix([0, 0, 0, 1]))
    return [sp.cancel(sp.together(x)) for x in w]


def limit(w):
    return [sp.simplify(sp.limit(x, e, 0, '+')) for x in w]


DET = [(a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]
# m1lib's genome order is LSB first over (p1,p2,p3,p4); reproduce it for readable names
DET = [(g0, g1, g2, g3) for g0 in (0, 1) for g1 in (0, 1) for g2 in (0, 1) for g3 in (0, 1)]


def run(corner, tag):
    p = (p1, corner[0], corner[1], p4)
    V = limit(stat(p, p))
    print("=" * 78)
    print("corner (p2,p3) = %s   self-play limit W = %s" % (corner, V))
    rows = []
    for s in DET:
        W = limit(stat([sp.Integer(x) for x in s], p))
        a = sp.simplify(W[1] - V[1])
        b = sp.simplify(W[2] - V[2])
        c = sp.simplify(W[3] - V[3])
        rows.append((s, a, b, c))
        gen = "".join(str(x) for x in s)
        print("  %s  a=%-28s b=%-28s c=%-28s" % (gen, a, b, c))
    return rows


if __name__ == "__main__":
    R01 = run((0, 1), "01")
    R10 = run((1, 0), "10")

    print("\n" + "=" * 78)
    print("The 16 inequalities a(u-1) + b v - c <= 0, simplified, corner (0,1):")
    seen = set()
    for s, a, b, c in R01:
        ex = sp.simplify(sp.together(a * (u - 1) + b * v - c))
        n, d = sp.fraction(sp.cancel(ex))
        # denominators are positive on 0 <= p1 < 1, 0 < p4 <= 1; report numerator
        key = sp.srepr(sp.simplify(n / sp.LC(sp.Poly(n, u, v)) if n != 0 else n))
        print("  %s : %s <= 0   [num %s ; den %s]"
              % ("".join(map(str, s)), sp.simplify(ex), sp.factor(n), sp.factor(d)))
