"""ROUTE C -- partners above the switch line.

Exact eps->0 limits of pi(sigma, p) for all 16 deterministic sigma against an
alternator p = (p1, c2, c3, p4) with (c2,c3) in {(0,1),(1,0)}, p1 and p4 symbolic.
"""
import sympy as sp
from itertools import product

a, d, e, u, v = sp.symbols('a d e u v', positive=True)
SWAP = (0, 2, 1, 3)
STATE = ['CC', 'CD', 'DC', 'DD']


def clamp(x):
    return e + (1 - 2*e)*x


def stat(pv, qv):
    """Stationary row vector of the pair chain, focal pv against co-player qv (unclamped in)."""
    P = [clamp(x) for x in pv]
    Q = [clamp(x) for x in qv]
    M = sp.zeros(4, 4)
    for s in range(4):
        A_, B_ = P[s], Q[SWAP[s]]
        M[s, 0] = A_*B_
        M[s, 1] = A_*(1 - B_)
        M[s, 2] = (1 - A_)*B_
        M[s, 3] = (1 - A_)*(1 - B_)
    L = (M.T - sp.eye(4))
    L = L[:3, :]
    L = L.col_join(sp.Matrix([[1, 1, 1, 1]]))
    rhs = sp.Matrix([0, 0, 0, 1])
    w = L.LUsolve(rhs)
    return [sp.cancel(sp.together(x)) for x in w]


def lim0(expr):
    return sp.simplify(sp.limit(sp.cancel(expr), e, 0, '+'))


def genome(s):
    return ''.join(str(x) for x in s)


rows = {}
for corner in [(0, 1), (1, 0)]:
    p = (a, sp.Integer(corner[0]), sp.Integer(corner[1]), d)
    # self play
    ws = stat(p, p)
    Ws = [lim0(x) for x in ws]
    print('=== corner (p2,p3) =', corner, ' self-play limit W =', Ws)
    for sig in product([0, 1], repeat=4):
        w = stat([sp.Integer(x) for x in sig], p)
        W = [lim0(x) for x in w]
        tot = sp.simplify(sum(W))
        pi = sp.simplify(W[0] + W[1]*u + W[2]*(1 + v))
        rows[(corner, sig)] = (W, pi)
        print(f'  sigma={genome(sig)}  W=({W[0]}, {W[1]}, {W[2]}, {W[3]})  sum={tot}')
        print(f'        pi = {sp.nsimplify(pi)}')

import pickle, sys
with open('/private/tmp/claude-501/-Users-martin-Documents-CL1-PartnersRivalsM1f/c1127d6a-179d-4b6f-be64-0786413868ce/scratchpad/routeC_rows.pkl','wb') as f:
    pickle.dump({k: (sp.srepr(sp.Matrix(vv[0])), sp.srepr(vv[1])) for k, vv in rows.items()}, f)
print('saved')
