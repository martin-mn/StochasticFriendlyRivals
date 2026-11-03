#!/usr/bin/env python3
"""Decisive check of the reduced Nash forms of partners_forms.py, in sympy.

For each of the 16 deviators the raw facet expression  a(u-1) + b v - c  (with a, b, c the
EXACT eps -> 0 limits) is divided by the corresponding reduced expression LHS - RHS.  The
quotient must be a rational function of p alone that is strictly POSITIVE on the open
cube; then the two inequalities are the same inequality.  The six (below) / four (above)
deviators whose facet is identically zero are reported as ties.

Runs the two limit computations from scratch (~3 min) so that nothing is taken on trust.
"""
import numpy as np
import sympy as sp

e, u, v = sp.symbols('e u v', real=True)
p1, p2, p3, p4 = sp.symbols('p1 p2 p3 p4', real=True)
SWAP = (0, 2, 1, 3)
DET = [(a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]
RNG = np.random.default_rng(3)


def stat(P, Q):
    P = [e + (1 - 2 * e) * x for x in P]
    Q = [e + (1 - 2 * e) * x for x in Q]
    M = sp.zeros(4, 4)
    for s in range(4):
        a, b = P[s], Q[SWAP[s]]
        M[s, 0], M[s, 1], M[s, 2], M[s, 3] = a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)
    A = (M.T - sp.eye(4)).as_mutable()
    A[3, :] = sp.Matrix([[1, 1, 1, 1]])
    return list(A.LUsolve(sp.Matrix([0, 0, 0, 1])))


def lim(w):
    return [sp.simplify(sp.limit(sp.cancel(sp.together(x)), e, 0, '+')) for x in w]


def positive_on_cube(expr, syms, n=30000):
    """Sample the open cube densely; return (min, max) of the quotient."""
    f = sp.lambdify(syms, expr, "numpy")
    X = RNG.uniform(1e-4, 1 - 1e-4, (len(syms), n))
    val = np.asarray(f(*X), float) * np.ones(n)
    return float(np.nanmin(val)), float(np.nanmax(val))


def verify(p, reduced, syms, tag):
    """reduced[s] is the expression whose sign must match a(u-1)+bv-c."""
    V = lim(stat(p, p))
    print("\n%s : self-play limit V = %s" % (tag, V))
    ties, ok = [], 0
    for s in DET:
        gen = "".join(map(str, s))
        W = lim(stat([sp.Integer(x) for x in s], p))
        a, b, c = (sp.simplify(W[i] - V[i]) for i in (1, 2, 3))
        raw = sp.cancel(sp.together(a * (u - 1) + b * v - c))
        if raw == 0:
            ties.append(gen)
            continue
        red = reduced[gen]
        q = sp.cancel(sp.simplify(raw / red))
        assert not q.has(u) and not q.has(v), "%s: quotient still depends on the game: %s" % (gen, q)
        lo, hi = positive_on_cube(q, syms)
        flag = "OK " if lo > 0 else "BAD"
        print("   %s %s  quotient %-40s  range [%.4g, %.4g]"
              % (flag, gen, sp.factor(q), lo, hi))
        assert lo > 0, "%s: quotient is not positive on the cube" % gen
        ok += 1
    print("   ties (facet identically zero): %s" % " ".join(ties))
    print("   %d non-trivial facets, all matched" % ok)
    return ties


if __name__ == "__main__":
    A = 1 - p2
    K = p3 + p4 - p3 * p4
    below = {
        "0000": p4 * v - A, "0100": p4 * v - A, "1000": p4 * v - A, "1100": p4 * v - A,
        "0110": A * (u - 1) + p3 * v, "0111": A * (u - 1) + p3 * v,
        "0010": p4 * (A * (u - 1) + v) - A * (1 - p3),
        "0011": A * (u - 1) + K * v - A * (1 - p3),
        "0001": A * (1 - p4) * (u - 1) + K * v - A,
        "0101": A * (1 - p4) * (u - 1) + p3 * v - A * p3,
    }
    tb = verify((sp.Integer(1), p2, p3, p4), below, (p2, p3, p4), "BELOW  p = (1,p2,p3,p4)")

    ab = {
        "0000": p4 * (1 + v - u) - (1 + u + v),
        "0100": p4 * (1 + v - u) - (1 + u + v),
        "1000": p4 * (1 + v - u) - (1 + u + v),
        "1101": p1 * (1 - u + v) - 2 * v,
        "1110": p1 * (1 - u + v) - 2 * v,
        "1111": p1 * (1 - u + v) - 2 * v,
        "1100": -((1 + p4) * u + (2 - p1) * v + 1 - p1 - p4),
        "0001": -((p1 * p4 + p4 + 1) * u + (1 + p4 - p1 * p4) * v - (p1 * p4 + p4 - 1)),
        "0101": -((p1 + p4 + 1) * u + (3 - p1 - p4) * v - (p1 + p4 - 1)),
        "0111": -((1 + p4 - p1 * p4) * u + (3 - 2 * p1 - p4 + p1 * p4) * v
                  - (2 * p1 + p4 - p1 * p4 - 1)),
        "0110": -((1 - p1) * (u + v + 1) + p4 * (u + v - 1)),
        "1001": -((1 - p1) * (u + v + 1) + p4 * (u + v - 1)),
    }
    ta = verify((p1, sp.Integer(0), sp.Integer(1), p4), ab, (p1, p4),
                "ABOVE  corner (0,1),  p = (p1,0,1,p4)")

    # the (1,0) sheet: the whole polygon is  max(u, 1+v) <= (1+u+v)/2, i.e. u - v = 1
    mir = {}
    V10 = lim(stat((p1, sp.Integer(1), sp.Integer(0), p4), (p1, sp.Integer(1), sp.Integer(0), p4)))
    print("\nABOVE  corner (1,0) : self-play limit V = %s" % V10)
    worst = []
    for s in DET:
        W = lim(stat([sp.Integer(x) for x in s], (p1, sp.Integer(1), sp.Integer(0), p4)))
        a, b, c = (sp.simplify(W[i] - V10[i]) for i in (1, 2, 3))
        worst.append((("".join(map(str, s))), sp.simplify(a * (u - 1) + b * v - c)))
    lock_cd = [g for g, x in worst if sp.simplify(x - (u - sp.Rational(1, 2) * (1 + u + v))) == 0]
    lock_dc = [g for g, x in worst if sp.simplify(x - ((1 + v) - sp.Rational(1, 2) * (1 + u + v))) == 0]
    print("   deviators scoring exactly u    (lock into CD): %s" % " ".join(lock_cd))
    print("   deviators scoring exactly 1+v  (lock into DC): %s" % " ".join(lock_dc))
    assert lock_cd and lock_dc
    print("   so max deviation = max(u, 1+v) and the sheet is Nash only on u - v = 1")
