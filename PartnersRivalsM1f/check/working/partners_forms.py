#!/usr/bin/env python3
"""The reduced closed forms of the Nash polygon, and their verification.

BELOW the line (p1 = 1, A := 1 - p2, so A in [0,1]).  The 16 inequalities collapse to
six distinct ones; six deviators (1001 WSLS, 1010 TFT, 1011, 1101, 1110, 1111 ALLC) give
the identically-zero facet 0 <= 0, i.e. they always tie with the resident.

  I    p4 v <= A                                              0000 0100 1000 1100
  II   A (u-1) + p3 v <= 0                                    0110 0111
  III  p4 [ A (u-1) + v ] <= A (1 - p3)                       0010
  IV   A (u-1) + K v <= A (1 - p3),   K := p3 + p4 - p3 p4    0011
  V    A (1-p4)(u-1) + K v <= A                               0001
  VI   A (1-p4)(u-1) + p3 v <= A p3                           0101

Every left side is <= 0 and every right side >= 0 when u <= 1 and v <= 0, which is the
proof that EVERY efficient p is a Nash equilibrium on that quadrant.

ABOVE the line, corner (p2,p3) = (0,1), free (p1,p4).  Six distinct inequalities; the
four alternating deviators (*,0,1,*) tie, and two more (0110 WSLS-mirror, 1001 WSLS) are
implied by u+v > 1:

  A'   p4 (1 + v - u) <= 1 + u + v                            0000 0100 1000
  B'   p1 (1 - u + v) <= 2 v                                  1101 1110 1111
  C'   (1+p4) u + (2-p1) v + 1 - p1 - p4 >= 0                 1100
  D'   (p1 p4 + p4 + 1) u + (1 + p4 - p1 p4) v >= p1 p4 + p4 - 1     0001
  E'   (p1 + p4 + 1) u + (3 - p1 - p4) v >= p1 + p4 - 1              0101
  F'   (1 + p4 - p1 p4) u + (3 - 2 p1 - p4 + p1 p4) v >= 2 p1 + p4 - p1 p4 - 1   0111
  (always true above the line:  (1-p1)(u+v+1) + p4 (u+v-1) >= 0     0110 1001)

ABOVE the line, corner (p2,p3) = (1,0): the deviator that answers CD with C locks the
resident into permanent defection and scores u; the one that answers DC with D locks it
into permanent cooperation and scores 1+v.  So max deviation = max(u, 1+v) >= (1+u+v)/2
with equality only on u - v = 1.  That sheet is a partner sheet on the ray T = S and
nowhere else.
"""
import numpy as np

import partners_lib as L

RNG = np.random.default_rng(7)


# ------------------------------------------------------------- reduced forms, below

def slack_below(A, p3, p4, u, v):
    """min slack of the six reduced inequalities, written as RHS - LHS."""
    K = p3 + p4 - p3 * p4
    return np.minimum.reduce([
        A - p4 * v,
        -(A * (u - 1) + p3 * v),
        A * (1 - p3) - p4 * (A * (u - 1) + v),
        A * (1 - p3) - (A * (u - 1) + K * v),
        A - (A * (1 - p4) * (u - 1) + K * v),
        A * p3 - (A * (1 - p4) * (u - 1) + p3 * v),
    ])


def slack_above01(p1, p4, u, v):
    """min slack of the six reduced inequalities at the corner (p2,p3) = (0,1)."""
    return np.minimum.reduce([
        (1 + u + v) - p4 * (1 + v - u),
        2 * v - p1 * (1 - u + v),
        (1 + p4) * u + (2 - p1) * v + 1 - p1 - p4,
        (p1 * p4 + p4 + 1) * u + (1 + p4 - p1 * p4) * v - (p1 * p4 + p4 - 1),
        (p1 + p4 + 1) * u + (3 - p1 - p4) * v - (p1 + p4 - 1),
        (1 + p4 - p1 * p4) * u + (3 - 2 * p1 - p4 + p1 * p4) * v
        - (2 * p1 + p4 - p1 * p4 - 1),
        (1 - p1) * (u + v + 1) + p4 * (u + v - 1),
    ])


# ------------------------------------------------------------------------- checking

def check(n=4000, ngame=40, eps=1e-9):
    """Compare the reduced forms with the raw 16-facet test, sign for sign."""
    worst_b = worst_a = 0.0
    nb = na = 0
    P = np.column_stack([np.ones(n), RNG.random((n, 3))])
    ab, bb, cb = _bank(P)                                       # (3,16,n)
    m = n // 2
    Q = np.empty((n, 4))
    Q[:, 0], Q[:, 3] = RNG.random(n), RNG.random(n)
    Q[:m, 1], Q[:m, 2] = 0.0, 1.0
    Q[m:, 1], Q[m:, 2] = 1.0, 0.0
    aa, ba, ca = _bank(Q)

    # Six facets are identically zero below the line and four above it, so the RAW min
    # slack is <= 0 always: the test is "raw >= -tol", not "raw > 0".  The reduced slack
    # is the raw one rescaled by strictly positive facet denominators, so only
    # MEMBERSHIP must agree.  Points within the tolerance band of either test are
    # excluded from the comparison and counted separately.
    tol, band, amb = 1e-9, 1e-9, 0
    for _ in range(ngame):
        u, v = RNG.uniform(-8, 8, 2)
        raw = (cb - ab * (u - 1) - bb * v).min(axis=0)
        red = slack_below(1 - P[:, 1], P[:, 2], P[:, 3], u, v)
        ok = (np.abs(raw) > band) & (np.abs(red) > band)
        amb += int((~ok).sum())
        nb += int(((raw[ok] >= -tol) != (red[ok] >= -tol)).sum())
        worst_b = max(worst_b, float(np.abs(raw[raw > tol]).max(initial=0.0)))

        if u + v > 1:
            rawa = (ca - aa * (u - 1) - ba * v).min(axis=0)
            reda = slack_above01(Q[:, 0], Q[:, 3], u, v)
            sel = (Q[:, 1] == 0)
            oka = sel & (np.abs(rawa) > band) & (np.abs(reda) > band)
            amb += int((sel & ~oka).sum())
            na += int(((rawa[oka] >= -tol) != (reda[oka] >= -tol)).sum())
    print("   raw slack never exceeds %.2e (the tie facets pin it at 0)" % worst_b)
    print("   %d of %d comparisons fell inside the +-%g band and were skipped"
          % (amb, 2 * n * ngame, band))
    return nb, na


def _bank1(P, eps):
    W = np.empty((16, P.shape[0], 4))
    for i in range(16):
        W[i] = L.stat_batch(np.tile(L.SIG[i], (P.shape[0], 1)), P, eps)
    V = L.stat_batch(P, P, eps)
    D = W - V[None]
    return np.array([D[:, :, 1], D[:, :, 2], D[:, :, 3]])


def _bank(P, eps=1e-6):
    """Richardson extrapolation of the facet coefficients to eps = 0.

    W(eps) is a RATIONAL function of eps, so W(eps) = W0 + W1 eps + O(eps^2) and
    (10 W(eps/10) - W(eps))/9 kills the linear term.  At eps = 1e-6 the two solves are
    still perfectly conditioned and the result is exact to ~1e-11, checked against the
    closed forms of sym_alt.py / sym_cc.py."""
    return (10 * _bank1(P, eps / 10) - _bank1(P, eps)) / 9


if __name__ == "__main__":
    nb, na = check()
    print("reduced-vs-raw sign disagreements: below %d, above(corner 01) %d" % (nb, na))

    # the (1,0) sheet: max deviation is exactly max(u, 1+v)
    bad = 0
    for _ in range(200):
        u, v = RNG.uniform(-8, 8, 2)
        if u + v <= 1:
            continue
        p1, p4 = RNG.random(), RNG.random()
        p = (p1, 1.0, 0.0, p4)
        a, b, c = _bank(np.asarray(p, float)[None])[:, :, 0]
        gap = float(np.min(c - a * (u - 1) - b * v))
        pred = (1 + u + v) / 2 - max(u, 1 + v)
        if abs(gap - pred) > 1e-6:
            bad += 1
            print("   (1,0) sheet mismatch at (%.3f,%.3f) p=(%.3f,%.3f): %.3e vs %.3e"
                  % (u, v, p1, p4, gap, pred))
    print("(1,0) sheet: max deviation = max(u,1+v) -- mismatches %d" % bad)
