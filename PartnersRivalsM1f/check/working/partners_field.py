#!/usr/bin/env python3
"""ROUTE D, step 1: the partner field over the square -8 <= u,v <= 8, numerically.

METHOD, stated honestly as a search.
  * A candidate resident p is EFFICIENT at the game (u,v) iff it lies in Eff_CC
    (u+v<1), Eff_alt (u+v>1) or Eff_line (u+v=1) -- the closed forms of item 2.
    Efficiency therefore depends on the game only through the side of the line, so
    ONE library of candidates serves every game on that side.
  * For a fixed p the Nash condition is 16 AFFINE inequalities in (u,v),
        a_s (u-1) + b_s v <= c_s,   (a,b,c) = W(s,p) - W(p,p) in coords (CD,DC,DD),
    whose coefficients are computed once per candidate at eps = 1e-9 (error O(eps),
    verified against the exact symbolic facets in sym_alt.py).  So the whole field is
    one matrix product: 16 x |library| x |games|.
  * LIBRARY.  Below the line: p1 = 1 and (p2,p3,p4) on a Sobol-like uniform sample of
    the open cube -- every such point is in Eff_CC.  Above the line: the two corners
    (p2,p3) in {(0,1),(1,0)} with (p1,p4) uniform in the open square -- every such
    point is in Eff_alt.  On the line: the below-line library (Eff_CC subset Eff_line).
  * RESOLUTION.  The reported field is therefore a LOWER bound on "a partner exists"
    and the reported fraction is a Monte-Carlo estimate of the partner set's relative
    volume inside the efficient set.  A game where the sample finds nothing is then
    re-attacked by a targeted search (partners_refine.py).
"""
import numpy as np

import partners_lib as L

RNG = np.random.default_rng(20260909)
EPS = 1e-9


def library_below(n):
    """p1 = 1, (p2,p3,p4) uniform in (0,1)^3 : all of these lie in Eff_CC."""
    q = RNG.random((n, 3))
    P = np.column_stack([np.ones(n), q])
    return P


def library_above(n):
    """(p2,p3) at a corner, (p1,p4) uniform in (0,1)^2 : all of these lie in Eff_alt."""
    m = n // 2
    q = RNG.random((n, 2))
    P = np.empty((n, 4))
    P[:, 0], P[:, 3] = q[:, 0], q[:, 1]
    P[:m, 1], P[:m, 2] = 0.0, 1.0
    P[m:, 1], P[m:, 2] = 1.0, 0.0
    return P


def facet_bank(P, eps=EPS):
    """(3, 16, n) array of the facet coefficients a, b, c for every p in P."""
    n = P.shape[0]
    W = np.empty((16, n, 4))
    for i in range(16):
        W[i] = L.stat_batch(np.tile(L.SIG[i], (n, 1)), P, eps)
    V = L.stat_batch(P, P, eps)
    D = W - V[None]                       # (16, n, 4)
    return D[:, :, 1], D[:, :, 2], D[:, :, 3]


def scan(P, U, V, chunk=2000):
    """For every game, the number of library points that are Nash there, and the best
    slack attained.  U, V are flat arrays of games."""
    a, b, c = facet_bank(P)
    ng = U.size
    cnt = np.zeros(ng, np.int64)
    best = np.full(ng, -np.inf)
    for i in range(0, P.shape[0], chunk):
        aa, bb, cc = a[:, i:i + chunk], b[:, i:i + chunk], c[:, i:i + chunk]
        # slack[s, j, g] = cc - aa (u-1) - bb v ; min over s
        s = (cc[:, :, None] - aa[:, :, None] * (U - 1.0)[None, None, :]
             - bb[:, :, None] * V[None, None, :]).min(axis=0)     # (chunk, ng)
        cnt += (s >= -1e-7).sum(axis=0)
        best = np.maximum(best, s.max(axis=0))
    return cnt, best


if __name__ == "__main__":
    import sys
    NG = int(sys.argv[1]) if len(sys.argv) > 1 else 161      # games per axis
    NP = int(sys.argv[2]) if len(sys.argv) > 2 else 8000     # library size per side

    g = np.linspace(-8.0, 8.0, NG)
    UU, VV = np.meshgrid(g, g, indexing="ij")
    U, V = UU.ravel(), VV.ravel()
    below = U + V < 1.0 - 1e-12
    above = U + V > 1.0 + 1e-12

    Pb, Pa = library_below(NP), library_above(NP)
    cb, sb = scan(Pb, U[below], V[below])
    ca, sa = scan(Pa, U[above], V[above])

    FRAC = np.zeros(U.size)
    BEST = np.full(U.size, -np.inf)
    FRAC[below], BEST[below] = cb / NP, sb
    FRAC[above], BEST[above] = ca / NP, sa

    print("games %d x %d on [-8,8]^2 ; library %d per side ; eps = %g" % (NG, NG, NP, EPS))
    print("BELOW the line: %d games; partner found at %d; min fraction %.5f at %s"
          % (below.sum(), (cb > 0).sum(), (cb / NP).min(),
             (U[below][np.argmin(cb)], V[below][np.argmin(cb)])))
    print("ABOVE the line: %d games; partner found at %d; min fraction %.5f at %s"
          % (above.sum(), (ca > 0).sum(), (ca / NP).min(),
             (U[above][np.argmin(ca)], V[above][np.argmin(ca)])))
    bad = np.where(FRAC == 0)[0]
    print("games with NO library partner: %d" % bad.size)
    for k in bad[:40]:
        print("   (u,v) = (%7.3f,%7.3f)  side %s  best slack %.3e"
              % (U[k], V[k], "above" if above[k] else "below", BEST[k]))
    np.savez_compressed("partners_field.npz", u=UU, v=VV,
                        frac=FRAC.reshape(NG, NG), best=BEST.reshape(NG, NG))
    print("wrote partners_field.npz")
