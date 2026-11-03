"""Route D engine: exact eps->0 limits of self-play distribution and rival margin."""
import sys, itertools
from fractions import Fraction as F
import sympy as sp

sys.path.insert(0, '/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib

EPS = sp.Symbol('epsilon', positive=True)
SWAP = (0, 2, 1, 3)
VERT = [tuple(b) for b in itertools.product([0, 1], repeat=4)]   # LSB-first genome order


def cl(x, eps=EPS):
    return eps + (1 - 2*eps)*sp.sympify(x)


# ---------------------------------------------------------------- self-play (lumped, closed form)
def selfplay_nums(p, eps=EPS):
    """(X_num, Y_num, Z_num) for the lumped self-play chain, as exact expressions in eps.
    X=W_CC, Y=W_CD+W_DC, Z=W_DD, each = num / (X_num+Y_num+Z_num)."""
    a, b, c, d = [cl(x, eps) for x in p]
    beta = (1-b)*(1-c)
    gam = b*c
    Xn = d*(d*beta + gam*(2-d))
    Yn = 2*d*(1-a)*(1+a-d)
    Zn = (1-a)*((1+a)*beta + gam*(1-a))
    return sp.expand(Xn), sp.expand(Yn), sp.expand(Zn)


def selfplay_limit(p):
    Xn, Yn, Zn = selfplay_nums(p)
    D = sp.expand(Xn + Yn + Zn)
    return tuple(sp.limit(sp.cancel(n/D), EPS, 0, '+') for n in (Xn, Yn, Zn))


# ---------------------------------------------------------------- full pair chain, exact in eps
def pair_stationary(p, q, eps=EPS):
    """Stationary row vector of the 4-state pair chain, as exact rational functions of eps."""
    pc = [cl(x, eps) for x in p]
    qc = [cl(x, eps) for x in q]
    M = sp.zeros(4, 4)
    for s in range(4):
        a, b = pc[s], qc[SWAP[s]]
        M[s, 0] = a*b; M[s, 1] = a*(1-b); M[s, 2] = (1-a)*b; M[s, 3] = (1-a)*(1-b)
    A = (M.T - sp.eye(4))
    A = A[:3, :]
    A = A.col_join(sp.ones(1, 4))
    rhs = sp.Matrix([0, 0, 0, 1])
    w = A.LUsolve(rhs)
    return [sp.cancel(sp.together(x)) for x in w]


def margin_limit(p, q):
    """lim_{eps->0} (w_DC - w_CD) for focal p against co-player q."""
    w = pair_stationary(p, q)
    return sp.limit(sp.cancel(w[2] - w[1]), EPS, 0, '+')


def rival_margin_limit(p, sign=1, qs=None):
    """min_q sign*(W_DC - W_CD); 0 iff p is a rival.  Returns (value, argmin q)."""
    qs = VERT if qs is None else qs
    best, arg = None, None
    for q in qs:
        m = sp.nsimplify(sign*margin_limit(p, q))
        if best is None or m < best:
            best, arg = m, q
    return best, arg


# ---------------------------------------------------------------- deficits (limits, convention A)
def efficiency_deficit_limit(p, u, v):
    """pi(p,p) - max(1,(1+u+v)/2) evaluated on the eps->0 limit.  <= 0."""
    X, Y, Z = selfplay_limit(p)
    u, v = sp.nsimplify(u), sp.nsimplify(v)
    pi = X + sp.Rational(1, 2)*Y*(1 + u + v)          # W_CD = W_DC = Y/2
    return sp.simplify(pi - sp.Max(sp.Integer(1), (1 + u + v)/2))


def joint_deficit(p, u, v):
    sign = 1 if (u - v) < 1 else -1
    e = efficiency_deficit_limit(p, u, v)
    r, arg = rival_margin_limit(p, sign)
    return sp.Min(e, r), e, r, arg
