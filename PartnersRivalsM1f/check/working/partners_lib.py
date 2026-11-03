"""Route D machinery: the Nash POLYGON of a memory-one resident, and the field of
partners over the plane of games.

THE ONE STRUCTURAL IDENTITY THIS FILE IS BUILT ON.  Let W(s,p) be the eps -> 0 limit
of the pair chain with deviator s focal against resident p, and V = W(p,p).  Then

    pi(s,p) - pi(p,p) = dCD (u-1) + dDC v - dDD,      d_X := W_X(s,p) - V_X,

because the four W's sum to one, so the constant term W_CC - V_CC = -(dCD+dDC+dDD)
eliminates R = 1 exactly.  Hence for a FIXED p the Nash set in the (u,v) plane is

    NASH(p) = { (u,v) : dCD(s)(u-1) + dDC(s) v <= dDD(s) for each of the 16 s },

an intersection of at most 16 half-planes: a closed convex POLYGON whose facet
coefficients are rational functions of p and do not depend on the game at all.
The 16-deviator quantifier is exact (LaPorte, Hilbe, Glynatsi & Nowak 2026; Levinsky,
Neyman & Zeleny 2020): against a memory-one resident the deviator's problem is a
unichain average-reward MDP over a box, whose optimum is at a vertex of the box.

Everything here is float at eps = 1e-9 unless it says exact; check_exact() re-does any
claim in Fractions.
"""
import numpy as np

from m1lib import SWAP, DET, clamp, stationary_exact, in_eff_cc, in_eff_alt, in_eff_line

EPS = 1e-9
SIG = np.array(DET)                      # (16,4) the deterministic deviators


# ------------------------------------------------------------------ batched chains

def stat_batch(P, Q, eps=EPS):
    """Stationary distributions for focal P[i] against co-player Q[i], batched.

    P, Q : (n,4) arrays of raw strategies (clamping applied here).
    returns (n,4) stationary row vectors.
    """
    P = clamp(np.atleast_2d(P), eps)
    Q = clamp(np.atleast_2d(Q), eps)
    n = P.shape[0]
    M = np.empty((n, 4, 4))
    for s in range(4):
        a, b = P[:, s], Q[:, SWAP[s]]
        M[:, s, 0] = a * b
        M[:, s, 1] = a * (1 - b)
        M[:, s, 2] = (1 - a) * b
        M[:, s, 3] = (1 - a) * (1 - b)
    A = np.transpose(M, (0, 2, 1)) - np.eye(4)[None]
    A[:, 3, :] = 1.0                              # replace the redundant 4th balance
    rhs = np.zeros((n, 4))
    rhs[:, 3] = 1.0
    return np.linalg.solve(A, rhs)


def facets(p, eps=EPS):
    """The 16 half-planes of NASH(p): arrays a, b, c with a(u-1) + b v <= c.

    a = W_CD(s,p) - V_CD,  b = W_DC(s,p) - V_DC,  c = W_DD(s,p) - V_DD.
    """
    p = np.asarray(p, float)
    W = stat_batch(SIG, np.tile(p, (16, 1)), eps)
    V = stat_batch(p[None], p[None], eps)[0]
    return W[:, 1] - V[1], W[:, 2] - V[2], W[:, 3] - V[3]


def nash_gap_uv(p, u, v, eps=EPS):
    """min over the 16 deviators of [ dDD - dCD(u-1) - dDC v ] = pi(p,p) - max pi(s,p).
    >= 0 (up to eps) means (p,p) is a limit Nash equilibrium at the game (u,v)."""
    a, b, c = facets(p, eps)
    return float(np.min(c - a * (u - 1.0) - b * v))


def selfW(p, eps=EPS):
    return stat_batch(np.asarray(p, float)[None], np.asarray(p, float)[None], eps)[0]


# ------------------------------------------------------------- efficiency by region

def is_efficient(p, u, v, tol=0.0):
    """The closed forms of item 2, selected by the side of the switch line."""
    s = u + v - 1.0
    if s < 0:
        return in_eff_cc(p, tol)
    if s > 0:
        return in_eff_alt(p, tol)
    return in_eff_line(p, tol)


# --------------------------------------------------------------------- exact route

def facets_exact(p, eps):
    """The same 16 (a,b,c) in exact Fractions at a given (rational) eps."""
    V = stationary_exact(p, p, eps)
    out = []
    for s in DET_EXACT_LIST:
        W = stationary_exact(s, p, eps)
        out.append((W[1] - V[1], W[2] - V[2], W[3] - V[3]))
    return out


from fractions import Fraction
from itertools import product
DET_EXACT_LIST = [tuple(Fraction(x) for x in b) for b in product([0, 1], repeat=4)]
