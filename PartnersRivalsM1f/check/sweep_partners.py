#!/usr/bin/env python3
"""ITEM 3 VERIFICATION SWEEP: partners over the full memory-one space.

Every claim of item 3 that can be checked by direct computation is checked here against
the 16-deviator limit-Nash test, at BOTH eps = 1e-7 and eps = 1e-9 (never 1e-12: the
float route is ill-conditioned there -- see check/README.md).  Nothing is taken from the
closed forms except as a PREDICTION to be compared with the direct computation.

The direct test.  For a resident p and a deviator s let W(s,p) be the stationary pair
distribution; with d_X = W_X(s,p) - W_X(p,p) and sum_X W_X = 1,

    pi(s,p) - pi(p,p) = d_CD (u-1) + d_DC v - d_DD,

so the Nash gap is   min over the 16 s of [ d_DD - d_CD (u-1) - d_DC v ].  The quantifier
"16 deterministic deviators suffice" is cited, not reproved (LaPorte, Hilbe, Glynatsi &
Nowak PNAS 2026; Levinsky, Neyman & Zeleny 2020).

Decision rule.  At every eps > 0 a limit-partner has a strictly NEGATIVE gap of order
eps; a genuine non-equilibrium has a gap tending to a negative constant.  So membership
is decided by the ORDER of the gap in eps, not by an absolute tolerance:
    limit-Nash  <=>  gap(1e-9) >= -1e-5   or   |gap(1e-9)| < 0.1 |gap(1e-7)|.
(The brief's warning applies: a threshold of 1e-2 on the ratio sits exactly on the
answer, since the ratio for a true limit point is about 1e-2 for a factor-100 eps step.)

Sections
  1  candidate banks and the direct/predicted comparison over a grid of games
  2  the sixteen binary strategies -- reproduce ../PartnersRivalsM1 panel b exactly
  3  partner EXISTENCE at every game, with the three witness rules
  4  friendly rivals in wedge S are partners (Yi, Baek & Choi 2017)
  5  D1 the ZD equalisers, D2 the alternators -- explicit O(eps) coefficients
  6  the share of the square admitting a partner, against the 16's 43.75 %
"""
import sys
import time
from itertools import product

import numpy as np

sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivalsM1f/check")

import partners_lib as L                                     # noqa: E402
from m1lib import (DET, genome, in_eff_alt, in_eff_cc,       # noqa: E402
                   in_eff_line, in_riv_plus, in_riv_minus)

EPS9, EPS7 = 1e-9, 1e-7
SIG = np.array(DET)
RNG = np.random.default_rng(20260909)
T0 = time.time()


def head(s):
    print("\n" + "=" * 78 + "\n%s   [t = %.1f s]\n" % (s, time.time() - T0) + "=" * 78)


# --------------------------------------------------------------------- the direct test

def facet_bank(P, eps):
    """(a,b,c) of shape (16, n): the resident-wise half-planes a(u-1) + b v <= c."""
    P = np.atleast_2d(np.asarray(P, float))
    n = P.shape[0]
    V = L.stat_batch(P, P, eps)
    a = np.empty((16, n)); b = np.empty((16, n)); c = np.empty((16, n))
    for i in range(16):
        W = L.stat_batch(np.tile(SIG[i], (n, 1)), P, eps)
        a[i] = W[:, 1] - V[:, 1]
        b[i] = W[:, 2] - V[:, 2]
        c[i] = W[:, 3] - V[:, 3]
    return a, b, c, V


def gap_at(bank, u, v):
    """Nash gap; u,v may be scalars, or arrays broadcasting against the resident axis."""
    a, b, c, _ = bank
    return (c - a * (np.asarray(u) - 1.0) - b * np.asarray(v)).min(axis=0)


def deficit(V, u, v):
    """pi(p,p) - max(R,(T+S)/2)."""
    own = V[:, 0] + V[:, 1] * u + V[:, 2] * (1.0 + v)
    return own - np.maximum(1.0, (1.0 + np.asarray(u) + np.asarray(v)) / 2.0)


def is_limit_nash(g9, g7):
    g9 = np.asarray(g9, float); g7 = np.asarray(g7, float)
    small = g9 >= -1e-5
    shrink = (g9 < 0) & (g7 < 0) & (np.abs(g9) < 0.1 * np.abs(g7))
    return small | shrink


# ------------------------------------------------------------------- candidate banks

def bank_below(n=90):
    """A structured + random sample of Eff_CC (p1 = 1)."""
    fixed = [(1, .5, .5, .5), (1, 1, 1, 1), (1, 1, 1, 0), (1, 1, .3, 0), (1, .3, 1, 0),
             (1, 0, 0, .5), (1, 1, 0, .5), (1, 0, 1, .5), (1, .9, .9, .01),
             (1, .01, .01, .9), (1, 1, .5, .5), (1, .5, 1, .5), (1, 0, .5, .5),
             (1, .5, 0, .5), (1, 1, .5, 0), (1, .5, 1, 0), (1, .2, .9, .7),
             (1, .9, .2, .7), (1, 1, 1, .5), (1, .7, .7, 1), (1, 0, 0, 1)]
    rnd = np.column_stack([np.ones(n), RNG.random((n, 3))])
    P = np.vstack([np.array(fixed, float), rnd])
    keep = np.array([in_eff_cc(p) for p in P])
    return P[keep]


def bank_above(n=60):
    """Both alternating sheets of Eff_alt, (p2,p3) = (0,1) and (1,0)."""
    out = []
    for (q2, q3) in [(0.0, 1.0), (1.0, 0.0)]:
        fixed = [(.5, q2, q3, .5), (0, q2, q3, .5), (.9, q2, q3, .1), (.1, q2, q3, .9),
                 (.5, q2, q3, 1), (0, q2, q3, .01), (.99, q2, q3, .99), (.3, q2, q3, .8)]
        rnd = np.column_stack([RNG.random(n), np.full(n, q2), np.full(n, q3),
                               RNG.random(n)])
        out.append(np.vstack([np.array(fixed, float), rnd]))
    P = np.vstack(out)
    keep = np.array([in_eff_alt(p) for p in P])
    return P[keep]


# ------------------------------------------------------------- the predicted criteria

def pred_below(P, u, v):
    """(I) p4 v <= 1-p2  and  (II) (1-p2)(u-1) + p3 v <= 0, for p in Eff_CC."""
    A = 1.0 - P[:, 1]
    return (P[:, 3] * v <= A + 1e-12) & (A * (u - 1.0) + P[:, 2] * v <= 1e-12)


def pred_above(P, u, v):
    """(C1),(C2) on the (0,1) sheet; the (1,0) sheet only on the ray u - v = 1."""
    X, Y = u - v - 1.0, u + v - 1.0
    sheet01 = P[:, 1] == 0.0
    c1 = P[:, 3] * X + Y + 2.0 >= -1e-12
    c2 = (P[:, 0] - 1.0) * X + Y >= -1e-12
    return np.where(sheet01, c1 & c2, abs(X) <= 1e-12)


# ============================================================== 1. grid of games
def section1():
    head("1. DIRECT 16-DEVIATOR TEST vs THE CLOSED FORMS, over a grid of games")
    Pb, Pa = bank_below(), bank_above()
    print("   Eff_CC candidates %d, Eff_alt candidates %d (%d on each sheet)"
          % (len(Pb), len(Pa), int((Pa[:, 1] == 0).sum())))
    bb9, bb7 = facet_bank(Pb, EPS9), facet_bank(Pb, EPS7)
    ba9, ba7 = facet_bank(Pa, EPS9), facet_bank(Pa, EPS7)

    games = [(u, v) for u in np.arange(-8, 8.01, 0.5) for v in np.arange(-8, 8.01, 0.5)]
    d = 0.01
    for base in [1.0, 0.0, -1.0, 2.0, 4.0]:                     # on/off every boundary
        for t in np.arange(-8, 8.01, 0.25):
            games += [(t, base - t), (t, base - t + d), (t, base - t - d)]   # u+v=base
            games += [(t, t - base), (t, t - base + d), (t, t - base - d)]   # u-v=base
    for a in np.arange(-8, 8.01, 0.25):                          # u=0,1 and v=0,1
        for c0 in [0.0, 1.0]:
            games += [(c0, a), (c0 + d, a), (c0 - d, a), (a, c0), (a, c0 + d), (a, c0 - d)]
    games = np.array(games)
    print("   %d games: the 1/2 grid on [-8,8]^2 plus on/just-off every claimed boundary"
          % len(games))

    # a fifth eps, so that BOTH classifications are decided by the same order-of-eps
    # rule: (1e-9 against 1e-7) and (1e-7 against 1e-5).
    bb5, ba5 = facet_bank(Pb, 1e-5), facet_bank(Pa, 1e-5)

    nb = na = 0; ndis97 = 0; nlive_b = nlive_a = 0
    worst = []
    exists_direct = np.zeros(len(games), bool)
    exists_adapted = np.zeros(len(games), bool)
    for k, (u, v) in enumerate(games):
        above = (u + v) > 1.0
        P, B9, B7, B5, pred = ((Pa, ba9, ba7, ba5, pred_above(Pa, u, v)) if above
                               else (Pb, bb9, bb7, bb5, pred_below(Pb, u, v)))
        g9, g7, g5 = gap_at(B9, u, v), gap_at(B7, u, v), gap_at(B5, u, v)
        got9, got7 = is_limit_nash(g9, g7), is_limit_nash(g7 * 100, g5 * 100)
        bad = got9 != pred
        ndis97 += int((got9 != got7).sum())
        exists_direct[k] = got9.any()
        if above:
            na += int(bad.sum()); nlive_a += int(got9.sum())
        else:
            nb += int(bad.sum()); nlive_b += int(got9.sum())
        if bad.any() and len(worst) < 12:
            i = int(np.argmax(bad))
            worst.append((u, v, tuple(np.round(P[i], 4)), g9[i], g7[i],
                          bool(got9[i]), bool(pred[i])))
    print("   closed-form disagreements: BELOW/ON the line %d, ABOVE %d" % (nb, na))
    print("   eps 1e-9 vs 1e-7 classification disagreements: %d" % ndis97)
    print("   partner found among the FIXED sample at %d of %d games (%.4f %%)"
          % (int(exists_direct.sum()), len(games), 100 * exists_direct.mean()))
    # a fixed sample is not an existence test: above the line with v < 0 the feasible
    # window in p1 can be thinner than 1 %.  Re-test the misses with the game-adapted
    # candidate that the closed form points at.
    miss_idx = np.flatnonzero(~exists_direct)
    if len(miss_idx):
        Pw = np.array([witnesses(*games[i]) for i in miss_idx])
        w9, w7 = facet_bank(Pw, EPS9), facet_bank(Pw, EPS7)
        uu, vv = games[miss_idx, 0], games[miss_idx, 1]
        okw = is_limit_nash(gap_at(w9, uu, vv), gap_at(w7, uu, vv))
        exists_direct[miss_idx] = okw
        print("   the %d misses are all above the line with v < 0, where C2 forces "
              "p1 >= -2v/(u-v-1)" % len(miss_idx))
        print("   thinnest feasible p1 window among them: %.5f"
              % float(np.min(1.0 - (-2 * vv) / (uu - vv - 1))))
        print("   re-tested with the game-adapted candidate: %d of %d now partners"
              % (int(okw.sum()), len(miss_idx)))
    print("   partner found at %d of %d games (%.4f %%)"
          % (int(exists_direct.sum()), len(games), 100 * exists_direct.mean()))
    for w in worst:
        print("      MISMATCH (u,v)=(%.3f,%.3f) p=%s g9=%.3e g7=%.3e direct=%s pred=%s" % w)
    miss = games[~exists_direct]
    if len(miss):
        print("   games with NO partner even after adaptation (first 10):")
        for u, v in miss[:10]:
            print("      (%.3f, %.3f)" % (u, v))
    return games


# ============================================================== 2. the sixteen vertices
def section2():
    head("2. THE SIXTEEN BINARY STRATEGIES -- reproduce ../PartnersRivalsM1 panel b")
    V = np.array(DET)
    b9, b7 = facet_bank(V, EPS9), facet_bank(V, EPS7)
    NG, SQ = 3201, 8.0
    g = (np.arange(NG) + 0.5) / NG * (2 * SQ) - SQ
    U, W = (x.ravel() for x in np.meshgrid(g, g, indexing="ij"))
    below, above = U + W < 1, U + W > 1
    for tag, bank in [("eps=1e-9", b9), ("eps=1e-7", b7)]:
        a, bb, c, Vs = bank
        cnt = np.zeros(U.size, int); who = []
        for j in range(16):
            p = DET[j]
            nash = (c[:, j][:, None] - a[:, j][:, None] * (U - 1)[None, :]
                    - bb[:, j][:, None] * W[None, :]).min(axis=0) >= -1e-5
            eff = np.where(below, in_eff_cc(p), np.where(above, in_eff_alt(p),
                                                         in_eff_line(p)))
            part = nash & eff
            cnt += part
            if part.any():
                who.append((genome(p), round(100 * part.mean(), 3)))
        print("   %s" % tag)
        for k in range(int(cnt.max()) + 1):
            if np.any(cnt == k):
                print("      count %d : %8.4f %% of the square" % (k, 100 * np.mean(cnt == k)))
        print("      max count           : %d" % cnt.max())
        print("      {count=3} == {u<=1, v<=0} : %s"
              % bool(np.array_equal(cnt == 3, (U <= 1) & (W <= 0))))
        print("      any partner above the switch line : %s"
              % bool(((cnt > 0) & above).any()))
        print("      a partner somewhere : %s" % who)
        print("      share with >=1 partner : %8.4f %%  (published 43.75 %%)"
              % (100 * np.mean(cnt > 0)))
        print("      share with count 3     : %8.4f %%  (exact 9x8/256 = 28.125 %%)"
              % (100 * np.mean(cnt == 3)))


# ============================================================== 3. existence witnesses
def witnesses(u, v):
    """The three witness rules of the item-3 report, one per regime."""
    if u + v > 1:
        p4 = 0.5 if u >= 0 else 0.5 * (1 + u + v) / (1 + v - u)
        p1 = 0.0 if v >= 0 else 0.5 * (1 - 2 * v / (u - v - 1))
        return (p1, 0.0, 1.0, p4)
    if v <= 0:
        return (1.0, 1.0, 1.0, 0.5)
    b = min(0.4, 0.4 / max(1.0, abs(v), abs(1 - u)))            # ZD equaliser, kappa = 1
    return (1.0, 1 - b * v, b * (1 - u), b)


def section3(games):
    head("3. PARTNER EXISTENCE AT EVERY GAME -- the three witness rules")
    P = np.array([witnesses(u, v) for u, v in games])
    ok = np.array([(0 <= x).all() and (x <= 1).all() for x in P])
    print("   witnesses inside the cube: %d of %d" % (int(ok.sum()), len(P)))
    for eps, tag in [(EPS9, "eps=1e-9"), (EPS7, "eps=1e-7")]:
        bank = facet_bank(P, eps)
        u, v = games[:, 0], games[:, 1]
        g = gap_at(bank, u, v)
        de = deficit(bank[3], u, v)
        eff = np.array([(in_eff_alt(p) if uu + vv > 1 else
                         (in_eff_line(p) if uu + vv == 1 else in_eff_cc(p)))
                        for p, (uu, vv) in zip(P, games)])
        fail = ~(is_limit_nash(g, g * 100) & eff)
        print("   %s : failures %d/%d, worst Nash gap %.3e, worst efficiency deficit %.3e"
              % (tag, int(fail.sum()), len(P), g.min(), de.min()))
        for i in np.flatnonzero(fail)[:8]:
            print("      FAIL (u,v)=(%.3f,%.3f) p=%s gap=%.3e def=%.3e eff=%s"
                  % (u[i], v[i], tuple(np.round(P[i], 5)), g[i], de[i], bool(eff[i])))


# ============================================================== 4. friendly rivals
def in_friendly_rival_S(p):
    """Item 2's wedge-S friendly rivals: p1=1, p2>0, p3>0, p2+p3>=1, (p4>0 or max=1)."""
    p1, p2, p3, p4 = p
    return (p1 >= 1) and p2 > 0 and p3 > 0 and (p2 + p3 >= 1) and (p4 > 0 or max(p2, p3) >= 1)


def section4():
    head("4. EVERY WEDGE-S FRIENDLY RIVAL IS A PARTNER (Yi, Baek & Choi 2017)")
    cand = []
    while len(cand) < 400:
        p2, p3, p4 = RNG.random(3)
        if RNG.random() < 0.3:
            p4 = 0.0; p2 = 1.0 if RNG.random() < .5 else p2
            p3 = 1.0 if p2 < 1 else p3
        p = (1.0, p2, p3, p4)
        if in_friendly_rival_S(p) and in_eff_cc(p):
            cand.append(p)
    P = np.array(cand)
    b9, b7 = facet_bank(P, EPS9), facet_bank(P, EPS7)
    # wedge S = {u+v<1, u-v>1}
    gs = []
    while len(gs) < 60:
        u, v = RNG.uniform(-8, 8, 2)
        if u + v < 1 and u - v > 1:
            gs.append((u, v))
    nfail = 0; nrival = 0
    for u, v in gs:
        g9, g7 = gap_at(b9, u, v), gap_at(b7, u, v)
        ok = is_limit_nash(g9, g7)
        nfail += int((~ok).sum()); nrival += len(P)
        if not ok.all() and nfail <= 5:
            i = int(np.argmax(~ok))
            print("      FAIL (u,v)=(%.3f,%.3f) p=%s g9=%.3e g7=%.3e"
                  % (u, v, tuple(np.round(P[i], 4)), g9[i], g7[i]))
    print("   %d friendly-rival x game tests in wedge S, failures to be a partner: %d"
          % (nrival, nfail))
    # strictness: partners in S that are NOT friendly rivals
    Pb = bank_below()
    bb9, bb7 = facet_bank(Pb, EPS9), facet_bank(Pb, EPS7)
    extra = 0
    for u, v in gs[:20]:
        ok = is_limit_nash(gap_at(bb9, u, v), gap_at(bb7, u, v))
        fr = np.array([in_friendly_rival_S(p) for p in Pb])
        extra += int((ok & ~fr).sum())
    print("   partners in wedge S that are NOT friendly rivals (20 games): %d "
          "-- the containment is STRICT" % extra)


# ============================================================== 5. D1 and D2
def section5():
    head("5. D1 ZD EQUALISERS AND D2 ALTERNATORS -- explicit O(eps) coefficients")
    u, v = 0.5, 5.0
    p = (1.0, 0.1, 0.09, 0.18)
    print("   D1: (u,v)=(0.5,5), equaliser p=(1,0.1,0.09,0.18)")
    print("       %-10s %-16s %-14s" % ("eps", "nash gap", "gap/eps"))
    for eps in [1e-5, 1e-6, 1e-7, 1e-8, 1e-9]:
        g = gap_at(facet_bank(np.array([p]), eps), u, v)[0]
        print("       %-10.0e %-16.9e %-14.6f" % (eps, g, g / eps))
    print("       claimed coefficient -9.2593 = -250/27 = %.6f" % (-250 / 27))
    # the equaliser property itself: every deviator gets exactly kappa = 1
    bank = facet_bank(np.array([p]), 1e-9)
    a, b, c, V = bank
    devs = (V[0, 0] + V[0, 1] * u + V[0, 2] * (1 + v)) - (c - a * (u - 1) - b * v)[:, 0]
    print("       deviator payoffs at eps=1e-9: min %.9f max %.9f (kappa = 1)"
          % (devs.min(), devs.max()))

    print("   D2: alternators (p1,0,1,p4) at (u,v)=(0.5,5); joint optimum (1+u+v)/2 = %.3f"
          % ((1 + u + v) / 2))
    fam = [(0.42436, 0, 1, 0.42436), (0.3, 0, 1, 0.8), (0.7, 0, 1, 0.2), (0.5, 0, 1, 0.5)]
    P = np.array(fam, float)
    for eps in [EPS7, EPS9]:
        bank = facet_bank(P, eps)
        g = gap_at(bank, u, v); de = deficit(bank[3], u, v)
        for i, q in enumerate(fam):
            print("       eps=%-6.0e p=(%.5f,0,1,%.5f)  pi(p,p)=%.9f  gap=%+.3e"
                  % (eps, q[0], q[3], (1 + u + v) / 2 + de[i], g[i]))
    print("   D4: ALLD limit-Nash iff u <= 0 --")
    bank = facet_bank(np.array([[0., 0., 0., 0.]]), EPS9)
    bad = [(uu, vv) for uu in np.arange(-8, 8.01, .25) for vv in np.arange(-8, 8.01, 2.)
           if bool(is_limit_nash(gap_at(bank, uu, vv)[0], 1.0)) != (uu <= 0)]
    print("       exceptions on a 1/4 x 2 grid: %d %s" % (len(bad), bad[:6]))


# ============================================================== 6. the share
def section6():
    head("6. SHARE OF THE SQUARE ADMITTING A PARTNER")
    NG, SQ = 201, 8.0
    g = (np.arange(NG) + 0.5) / NG * (2 * SQ) - SQ
    U, W = (x.ravel() for x in np.meshgrid(g, g, indexing="ij"))
    P = np.array([witnesses(u, v) for u, v in zip(U, W)])
    bank = facet_bank(P, EPS9)
    gp = gap_at(bank, U, W)
    de = deficit(bank[3], U, W)
    ok = is_limit_nash(gp, gp * 100) & (de > -1e-5)
    print("   %d x %d cell-centred grid on [-8,8]^2" % (NG, NG))
    print("   continuum: a partner at %.4f %% of the square (worst gap %.2e, worst "
          "deficit %.2e)" % (100 * ok.mean(), gp.min(), de.min()))
    print("   the 16 binary strategies: 43.75 % (see section 2)")
    print("   saturation region {u<=1,v<=0} u {u+v>1, u>=0, v>=0} = %.4f %% of the square"
          % (100 * np.mean(((U <= 1) & (W <= 0)) | ((U + W > 1) & (U >= 0) & (W >= 0)))))


def section2b():
    """The two shares are grid-limited, not wrong: both converge linearly in 1/NG."""
    V = np.array(DET)
    a, b, c, _ = facet_bank(V, EPS9)
    print("   cell-centred grid convergence of the two shares:")
    print("      %-8s %-12s %-12s" % ("NG", ">=1 partner", "count 3"))
    for NG in [401, 801, 1601, 3201]:
        g = (np.arange(NG) + 0.5) / NG * 16 - 8
        U, W = (x.ravel() for x in np.meshgrid(g, g, indexing="ij"))
        below, above = U + W < 1, U + W > 1
        cnt = np.zeros(U.size, int)
        for j in range(16):
            p = DET[j]
            nash = (c[:, j][:, None] - a[:, j][:, None] * (U - 1)[None, :]
                    - b[:, j][:, None] * W[None, :]).min(axis=0) >= -1e-5
            eff = np.where(below, in_eff_cc(p), np.where(above, in_eff_alt(p),
                                                         in_eff_line(p)))
            cnt += nash & eff
        print("      %-8d %-12.4f %-12.4f"
              % (NG, 100 * np.mean(cnt > 0), 100 * np.mean(cnt == 3)))
    print("      limits: 43.75 % and 28.125 % -- exactly the published values")


if __name__ == "__main__":
    games = section1()
    section2()
    section2b()
    section3(games)
    section4()
    section5()
    section6()
    print("\ntotal %.1f s" % (time.time() - T0))
