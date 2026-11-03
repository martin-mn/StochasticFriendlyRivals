#!/usr/bin/env python3
"""BELOW the switch line: the exact Nash polygon of an efficient mutual cooperator.

Efficiency below u+v = 1 forces p1 = 1 (item 2, Eff_CC), with (p2,p3) not an alternating
corner and (p4 > 0 or max(p2,p3) = 1).  Self-play is then W = (1,0,0,0), pi(p,p) = 1, and
the Nash condition for the deviator s is

    W_CD(s,p) (u-1) + W_DC(s,p) v  <=  W_DD(s,p).

All three coefficients are NON-NEGATIVE, which already proves that every efficient p is
Nash on the whole quadrant u <= 1, v <= 0.  This script computes the 16 triples exactly,
with eps symbolic and p2, p3, p4 symbolic.
"""
import sympy as sp

e, p2, p3, p4, u, v = sp.symbols('e p2 p3 p4 u v', nonnegative=True)
SWAP = (0, 2, 1, 3)


def clamp(x):
    return e + (1 - 2 * e) * x


def stat(P, Q):
    P = [clamp(x) for x in P]
    Q = [clamp(x) for x in Q]
    M = sp.zeros(4, 4)
    for s in range(4):
        a, b = P[s], Q[SWAP[s]]
        M[s, 0], M[s, 1], M[s, 2], M[s, 3] = a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)
    A = (M.T - sp.eye(4)).as_mutable()
    A[3, :] = sp.Matrix([[1, 1, 1, 1]])
    return list(A.LUsolve(sp.Matrix([0, 0, 0, 1])))


def lim(w):
    return [sp.simplify(sp.limit(sp.cancel(sp.together(x)), e, 0, '+')) for x in w]


DET = [(a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]

if __name__ == "__main__":
    p = (sp.Integer(1), p2, p3, p4)
    V = lim(stat(p, p))
    print("self-play limit V =", V)
    print()
    rows = []
    for s in DET:
        W = lim(stat([sp.Integer(x) for x in s], p))
        a, b, c = [sp.factor(sp.simplify(W[i] - V[i])) for i in (1, 2, 3)]
        rows.append((s, a, b, c))
        print("%s  a(CD)=%-34s b(DC)=%-34s c(DD)=%-34s"
              % ("".join(map(str, s)), a, b, c))
    print()
    print("the 16 constraints a(u-1) + b v - c <= 0:")
    for s, a, b, c in rows:
        ex = sp.cancel(sp.together(a * (u - 1) + b * v - c))
        n, d = sp.fraction(ex)
        print("  %s : num %s  <=0  [den %s]"
              % ("".join(map(str, s)), sp.factor(sp.expand(n)), sp.factor(d)))
