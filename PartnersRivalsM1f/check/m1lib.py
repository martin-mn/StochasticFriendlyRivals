"""Shared harness for the FULL (continuous) memory-one strategy space.

Conventions are those of ../../PartnersRivalsM1 and are not redefined anywhere else:

  * base game normalised to (R,S,T,P) = (1, u, 1+v, 0);  T-S = 1+v-u,  2R-(T+S) = 1-u-v
  * a memory-one strategy is p = (p1,p2,p3,p4), the probabilities of playing C after the
    outcomes CC, CD, DC, DD, written FOCAL PLAYER'S ACTION FIRST; states CC=0 CD=1 DC=2 DD=3
  * implementation error eps acts as p -> eps + (1-2 eps) p, so a strategy is a fixed point
    of the closed cube [0,1]^4 and eps is only the error floor ("convention A").  The 16
    binary strategies of ../../PartnersRivalsM1 are the 16 vertices of that cube.
  * w = (w_CC, w_CD, w_DC, w_DD) is the stationary distribution of the four-state pair chain;
    W is its eps -> 0 limit.  pi(p,q) = w.(1,u,1+v,0) and pi(q,p) = w.(1,1+v,u,0), so
    pi(p,q) - pi(q,p) = (T-S) (w_DC - w_CD).

Everything here is available in two arithmetics: float (numpy, fast, for sweeps) and exact
rational (Fraction, slow, for receipts).  The exact route uses no threshold anywhere.
"""

from fractions import Fraction
from itertools import product

import numpy as np

# outcome CC=0 CD=1 DC=2 DD=3 seen from the focal player; SWAP maps a state to the same
# state read from the co-player's side, which is where the co-player's response lives.
SWAP = (0, 2, 1, 3)

# the 16 binary strategies, in ../../PartnersRivalsM1's genome order (LSB first, bit=1 -> C)
DET = [np.array(b, float) for b in product([0, 1], repeat=4)]
DET_EXACT = [tuple(Fraction(x) for x in b) for b in product([0, 1], repeat=4)]


def genome(p):
    """'1010' for TFT: the genome string of a vertex, LSB first."""
    return "".join("1" if x > 0.5 else "0" for x in p)


NAMES = {"0000": "ALLD", "1111": "ALLC", "1010": "TFT", "1001": "WSLS", "1000": "Grim",
         "0101": "ATFT", "0110": "AWSLS"}


# --------------------------------------------------------------------------- float route

def clamp(p, eps):
    return eps + (1 - 2 * eps) * np.asarray(p, float)


def transition(p, q, eps):
    """Row-stochastic 4x4 matrix of the pair chain for focal p against co-player q."""
    p, q = clamp(p, eps), clamp(q, eps)
    M = np.empty((4, 4))
    for s in range(4):
        a, b = p[s], q[SWAP[s]]
        M[s] = (a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b))
    return M


def stationary(p, q, eps):
    """Unique stationary distribution; the chain is ergodic for 0 < eps < 1/2."""
    M = transition(p, q, eps)
    A = np.vstack([M.T - np.eye(4), np.ones(4)])
    return np.linalg.lstsq(A, np.array([0.0, 0.0, 0.0, 0.0, 1.0]), rcond=None)[0]


def payoff(w, u, v):
    """Focal payoff.  The co-player's is payoff_co below."""
    return w[0] + w[1] * u + w[2] * (1 + v)


def payoff_co(w, u, v):
    return w[0] + w[1] * (1 + v) + w[2] * u


def rival_margin(p, sign, eps, coplayers=None):
    """min over co-players of sign*(w_DC - w_CD).  >= 0 (up to eps) means p is a rival.

    sign = +1 for the half-plane T > S (u - v < 1), sign = -1 for T < S.
    The default co-player set is the 16 vertices, which by the MDP argument realises the
    minimum over the whole box [eps,1-eps]^4 -- see ../Text/rivals.md.
    """
    qs = DET if coplayers is None else coplayers
    return min(sign * (w[2] - w[1]) for w in (stationary(p, q, eps) for q in qs))


def self_play(p, eps):
    return stationary(p, p, eps)


def nash_gap(p, u, v, eps, deviators=None):
    """pi(p,p) - max over deviators of pi(sigma,p).  >= 0 (up to eps) means limit-Nash."""
    qs = DET if deviators is None else deviators
    own = payoff(self_play(p, eps), u, v)
    return own - max(payoff(stationary(s, p, eps), u, v) for s in qs)


def efficiency_deficit(p, u, v, eps):
    """pi(p,p) minus the largest average payoff two players can attain: max(R, (T+S)/2)."""
    return payoff(self_play(p, eps), u, v) - max(1.0, (u + 1 + v) / 2)


# --------------------------------------------------------------------------- exact route

def clamp_exact(p, eps):
    eps = Fraction(eps)
    return [eps + (1 - 2 * eps) * Fraction(x) for x in p]


def stationary_exact(p, q, eps):
    """Stationary distribution in exact rationals by Gaussian elimination on Fractions."""
    p, q = clamp_exact(p, eps), clamp_exact(q, eps)
    M = [[None] * 4 for _ in range(4)]
    for s in range(4):
        a, b = p[s], q[SWAP[s]]
        M[s] = [a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)]
    # w is a ROW vector: w M = w.  Equation i reads sum_j w_j (M[j][i] - delta_ij) = 0; three
    # of those plus the normalisation sum_j w_j = 1 determine w.
    A = [[M[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(4)]
         for i in range(3)] + [[Fraction(1)] * 4]
    b = [Fraction(0), Fraction(0), Fraction(0), Fraction(1)]
    n = 4
    for c in range(n):                                   # forward elimination, exact
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        b[c], b[piv] = b[piv], b[c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
                b[r] -= f * b[c]
    w = [Fraction(0)] * n
    for r in reversed(range(n)):                          # back substitution
        w[r] = (b[r] - sum(A[r][k] * w[k] for k in range(r + 1, n))) / A[r][r]
    return w


def akin_residual_exact(p, q, eps):
    """Akin's lemma: w . (p1-1, p2-1, p3, p4) = 0 with the CLAMPED p.  Must be exactly 0."""
    w = stationary_exact(p, q, eps)
    pc = clamp_exact(p, eps)
    return w[0] * (pc[0] - 1) + w[1] * (pc[1] - 1) + w[2] * pc[2] + w[3] * pc[3]


# --------------------------------------------------------------- the conjectured closed form

def in_riv_plus(p, tol=0.0):
    """Riv+ (T > S):  p4 = 0,  p2 + p3 <= 1,  and NOT (p2 = 1 and p3 = 0)."""
    p1, p2, p3, p4 = p
    return p4 <= tol and p2 + p3 <= 1 + tol and not (p2 >= 1 - tol and p3 <= tol)


def in_riv_minus(p, tol=0.0):
    """Riv- (T < S):  p1 = 1,  p2 + p3 >= 1,  and NOT (p2 = 1 and p3 = 0)."""
    p1, p2, p3, p4 = p
    return p1 >= 1 - tol and p2 + p3 >= 1 - tol and not (p2 >= 1 - tol and p3 <= tol)


def in_equaliser(p, tol=0.0):
    """Riv+ cap Riv-: the generalised-tit-for-tat line p = (1, s, 1-s, 0), 0 <= s < 1."""
    return in_riv_plus(p, tol) and in_riv_minus(p, tol)


# ---------------------------------------------------------------- efficiency (self-play)
#
# Emax = max{R, (T+S)/2} = max{1, (1+u+v)/2}, so with W the eps -> 0 self-play limit
#   pi(p,p) - Emax = -W_DD - (W_CD+W_DC)(1-u-v)/2          below the line, u+v < 1
#                  = -W_DD                                  on   the line, u+v = 1
#                  = W_CC(1-E) - W_DD E,  E = (1+u+v)/2 > 1 above the line, u+v > 1
# so p is efficient iff  W_CC = 1  /  W_DD = 0  /  W_CD = W_DC = 1/2  respectively.  The
# three closed forms below are proved in ../Text/efficiency.md; they depend on the game
# only through sign(u+v-1).

def _alt_corner(p2, p3, tol=0.0):
    """(p2,p3) is one of the two alternating corners (0,1) and (1,0): the locus where the
    two-cycle CD <-> DC is closed at eps = 0, its escape rate p2 p3 + (1-p2)(1-p3) = 0."""
    return ((p2 <= tol and p3 >= 1 - tol) or (p2 >= 1 - tol and p3 <= tol))


def in_eff_cc(p, tol=0.0):
    """Efficient BELOW the switch line (u+v < 1):  W_CC = 1, permanent mutual cooperation.

    Eff_CC = { p1 = 1,  (p2,p3) not in {(0,1),(1,0)},  and ( p4 > 0 or max(p2,p3) = 1 ) }.

    p1 = 1 makes CC leak only at rate 2eps-eps^2; the two exclusions are the two traps that
    hold the excursion for Theta(1/eps) steps and so survive the limit -- the alternating
    trap (p2,p3) in {(0,1),(1,0)} and, when p4 = 0, the DD trap, which the excursion avoids
    only if it never reaches DD, i.e. (1-p2)(1-p3) = 0.
    """
    p1, p2, p3, p4 = p
    if p1 < 1 - tol or _alt_corner(p2, p3, tol):
        return False
    return p4 > tol or max(p2, p3) >= 1 - tol


def in_eff_alt(p, tol=0.0):
    """Efficient ABOVE the switch line (u+v > 1):  W_CD = W_DC = 1/2, perfect alternation.

    Eff_alt = { (p2,p3) in {(0,1),(1,0)},  p1 < 1,  p4 > 0,  not (p1 = 0 and p4 = 1) }.

    Non-empty and TWO-dimensional -- it lies in the two 2-faces (p2,p3) = (0,1) and (1,0),
    with only p1 and p4 free -- so the "no alternator" theorem of the 16 binary
    strategies is a vertex artefact: p1 = 1 is excluded (CC is then sticky), p4 = 0 is
    excluded (DD is then sticky), and (p1,p4) = (0,1) is excluded because CC -> DD -> CC is
    then a second closed two-cycle that splits the limit mass.
    """
    p1, p2, p3, p4 = p
    if not _alt_corner(p2, p3, tol):
        return False
    if p1 >= 1 - tol or p4 <= tol:
        return False
    return not (p1 <= tol and p4 >= 1 - tol)


def in_eff_line(p, tol=0.0):
    """Efficient ON the switch line (u+v = 1):  W_DD = 0.  Strictly weaker than either of
    the two above, and their disjoint union is a proper subset of it:

    Eff_line = { p1 = 1 and ( p4 > 0 or (max(p2,p3) = 1 and min(p2,p3) > 0) ) }  cup  Eff_alt.
    """
    p1, p2, p3, p4 = p
    if p1 >= 1 - tol:
        return p4 > tol or (max(p2, p3) >= 1 - tol and min(p2, p3) > tol)
    return in_eff_alt(p, tol)
