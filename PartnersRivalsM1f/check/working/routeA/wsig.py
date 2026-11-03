"""Symbolic W(sigma,p) for the 16 deterministic deviators sigma against a generic resident p.

Convention A: both players clamped, p -> eps + (1-2eps)p.  States CC=0 CD=1 DC=2 DD=3 with the
FOCAL (= deviator sigma) action first, so the co-player's response at state s is p[SWAP[s]].
"""
import sympy as sp
from itertools import product

SWAP = (0, 2, 1, 3)
e = sp.symbols('epsilon', positive=True)
p = sp.symbols('p1 p2 p3 p4', positive=True)


def clamp(x):
    return e + (1 - 2 * e) * x


def stat(sig, res):
    """Stationary row vector of the pair chain, focal sig (list of 4), co-player res (list of 4)."""
    a = [clamp(sig[s]) for s in range(4)]
    b = [clamp(res[SWAP[s]]) for s in range(4)]
    M = sp.Matrix(4, 4, lambda s, t: 0)
    for s in range(4):
        M[s, 0] = a[s] * b[s]
        M[s, 1] = a[s] * (1 - b[s])
        M[s, 2] = (1 - a[s]) * b[s]
        M[s, 3] = (1 - a[s]) * (1 - b[s])
    A = (M.T - sp.eye(4))
    A = A[:3, :].col_join(sp.Matrix([[1, 1, 1, 1]]))
    rhs = sp.Matrix([0, 0, 0, 1])
    w = A.LUsolve(rhs)
    return [sp.cancel(sp.together(x)) for x in w]


if __name__ == '__main__':
    out = {}
    for bits in product([0, 1], repeat=4):
        g = ''.join(str(b) for b in bits)
        w = stat(list(bits), list(p))
        lim = [sp.simplify(sp.limit(x, e, 0, '+')) for x in w]
        out[g] = lim
        print(g, [sp.factor(sp.cancel(x)) for x in lim])
    import pickle
    with open('/private/tmp/claude-501/-Users-martin-Documents-CL1-PartnersRivalsM1f/c1127d6a-179d-4b6f-be64-0786413868ce/scratchpad/A/wsig.pkl', 'wb') as f:
        pickle.dump({k: [sp.srepr(x) for x in v] for k, v in out.items()}, f)
