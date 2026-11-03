#!/usr/bin/env python3
"""End-to-end validation of the drawn field against the RAW 16-deviator Nash test.

Nothing in section 1 uses the reduced inequalities.  A library of efficient strategies is
drawn uniformly; each one's 16 facet triples (a,b,c) are computed from the pair chain and
extrapolated to eps = 0; the game enters only through the affine test
a(u-1) + b v <= c.  The RAW mask is compared with the REDUCED mask POINT BY POINT -- the
strongest form of the comparison, and the one that does not hide behind an average.

THE NUMERICAL PRECAUTION, and it is exactly the one the project setup warns about.  The
float Richardson pair is solved at eps = 1e-6 and 1e-7, where the linear system's
condition number is of order 1/(eps p_min); the facet COEFFICIENTS come out right to
~1e-11, but the slack, which multiplies them by |u-1| <= 9 and |v| <= 8, carries a
residual up to ~1e-6 for a resident with a coordinate at 0.02 and far worse for one at
1e-5.  So section 1 draws its coordinates from [0.02, 0.98] and tests at tolerance 1e-5,
and section 1b redoes the whole comparison in EXACT Fractions, where there is no roundoff
at all and the only error is the O(eps^2) truncation of the extrapolation.

The identity of the raw and the reduced inequalities on the WHOLE cube is not a numerical
claim in any case: sym_verify.py proves it, by dividing each raw facet by its reduced form
and certifying that the quotient is a rational function of p alone, positive on the open
cube.
"""
from fractions import Fraction

import numpy as np

import partners_forms as F
import partners_volume as PV
from m1lib import DET_EXACT, in_eff_alt, in_eff_cc, in_eff_line, stationary_exact

RNG = np.random.default_rng(31415)
N = 20000
D = 0.02


def facets_exact_limit(p, e0=Fraction(1, 10 ** 8)):
    """The 16 facet triples as exact Fractions, Richardson-extrapolated to eps = 0."""
    out = []
    for e in (e0, e0 / 10):
        V = stationary_exact(p, p, e)
        out.append([tuple(stationary_exact(s, p, e)[i] - V[i] for i in (1, 2, 3))
                    for s in DET_EXACT])
    return [tuple((10 * y[i] - x[i]) / 9 for i in range(3))
            for x, y in zip(out[0], out[1])]


def main():
    lo, hi = D, 1 - D
    Pb = np.column_stack([np.ones(N), RNG.uniform(lo, hi, (N, 3))])   # Eff_CC
    Pa = np.empty((N, 4))                                             # Eff_alt, sheet (0,1)
    Pa[:, 0], Pa[:, 3] = RNG.uniform(lo, hi, N), RNG.uniform(lo, hi, N)
    Pa[:, 1], Pa[:, 2] = 0.0, 1.0
    Pm = Pa.copy()                                                    # sheet (1,0)
    Pm[:, 1], Pm[:, 2] = 1.0, 0.0
    assert all(in_eff_cc(p) for p in Pb[:200])
    assert all(in_eff_alt(p) for p in Pa[:200]) and all(in_eff_alt(p) for p in Pm[:200])
    assert all(in_eff_line(p) for p in Pb[:200]), "Eff_CC is contained in Eff_line"
    ab, aa, am = F._bank(Pb), F._bank(Pa), F._bank(Pm)

    def raws(bank, u, v):
        a, b, c = bank
        return (c - a * (u - 1.0) - b * v).min(axis=0)

    def rawmask(bank, u, v, tol=1e-5):
        return raws(bank, u, v) >= -tol

    print("1.  RAW 16-deviator mask vs REDUCED mask, point by point, %d points per game, "
          "coordinates in [%.2f, %.2f]" % (N, lo, hi))
    dis_b = dis_a = ngb = nga = band_b = band_a = 0
    floor = (0.0, 0.0)
    for u in np.arange(-8, 8.01, 0.5):
        for v in np.arange(-8, 8.01, 0.5):
            s = u + v
            if abs(s - 1) < 1e-12:
                continue
            if s < 1:
                raw = raws(ab, u, v)
                red = F.slack_below(1 - Pb[:, 1], Pb[:, 2], Pb[:, 3], u, v)
                bad = (raw >= -1e-5) != (red >= 0)
                dis_b += int((bad & (np.abs(red) > 1e-3)).sum())
                band_b += int((bad & (np.abs(red) <= 1e-3)).sum())
                soft = (raw >= -1e-8) != (red >= 0)
                if soft.any():
                    k = int(np.argmax(np.where(soft, np.abs(raw), -1.0)))
                    if abs(raw[k]) > floor[0]:
                        floor = (abs(raw[k]), abs(red[k]))
                ngb += 1
            else:
                raw = raws(aa, u, v)
                red = F.slack_above01(Pa[:, 0], Pa[:, 3], u, v)
                bad = (raw >= -1e-5) != (red >= 0)
                dis_a += int((bad & (np.abs(red) > 1e-3)).sum())
                band_a += int((bad & (np.abs(red) <= 1e-3)).sum())
                nga += 1
    print("    below the line: %d genuine disagreements out of %d, plus %d points whose"
          " reduced slack is inside the +-1e-3 boundary band" % (dis_b, ngb * N, band_b))
    print("    above the line: %d genuine disagreements out of %d, plus %d in the band"
          % (dis_a, nga * N, band_a))
    print("    the float noise floor, measured: at tolerance 1e-8 the worst raw slack at a"
          " disagreement is %.2e and the REDUCED slack there is %.3f -- the two masks part"
          " only where the float residual is ~1e-6 and the true margin is O(1)" % floor)
    assert dis_b == 0 and dis_a == 0

    print("\n1b. the same comparison in EXACT rationals -- no roundoff at all")
    rng2 = np.random.default_rng(271828)
    ndis = ntest = 0
    worst_tie = Fraction(0)
    for _ in range(120):
        q = [Fraction(int(x), 1000) for x in rng2.integers(20, 981, 3)]
        p = (Fraction(1), q[0], q[1], q[2])
        fac = facets_exact_limit(p)
        A_, K = 1 - q[0], q[1] + q[2] - q[1] * q[2]
        for _ in range(12):
            uu = Fraction(int(rng2.integers(-8000, 8001)), 1000)
            vv = Fraction(int(rng2.integers(-8000, 8001)), 1000)
            if uu + vv >= 1:
                continue
            raw = min(c - a_ * (uu - 1) - b_ * vv for a_, b_, c in fac)
            red = min([A_ - q[2] * vv,
                       -(A_ * (uu - 1) + q[1] * vv),
                       A_ * (1 - q[1]) - q[2] * (A_ * (uu - 1) + vv),
                       A_ * (1 - q[1]) - (A_ * (uu - 1) + K * vv),
                       A_ - (A_ * (1 - q[2]) * (uu - 1) + K * vv),
                       A_ * q[1] - (A_ * (1 - q[2]) * (uu - 1) + q[1] * vv)])
            ntest += 1
            if (raw >= -Fraction(1, 10 ** 10)) != (red >= 0):
                ndis += 1
            if red >= 0:
                worst_tie = max(worst_tie, abs(raw))
    print("    %d exact (resident, game) pairs, %d disagreements; largest |raw slack| at a"
          " Nash point %.2e -- the O(eps^2) truncation of the exact extrapolation"
          % (ntest, ndis, float(worst_tie)))
    assert ndis == 0

    print("\n2.  the closed-form field vs the raw mask's own fraction")
    worst = (0.0, None)
    for u in np.arange(-8, 8.01, 1.0):
        for v in np.arange(-8, 8.01, 1.0):
            s = u + v
            if abs(s - 1) < 1e-12:
                continue
            obs = float(rawmask(ab if s < 1 else aa, u, v).mean())
            cf = (float(PV.vol_below(np.array(u), np.array(v))) if s < 1
                  else float(PV.area_above(np.array(u), np.array(v))))
            if abs(obs - cf) > worst[0]:
                worst = (abs(obs - cf), (u, v, obs, cf))
    print("    largest |observed - closed form| = %.4f at (u,v) = (%.1f,%.1f): %.5f vs"
          " %.5f   (the library lives in the shrunken box, which costs O(D) = O(0.02))"
          % ((worst[0],) + worst[1]))

    print("\n3.  A PARTNER AT EVERY GAME -- an explicit witness, tested in EXACT rationals")
    #   The uniform library is drawn from the shrunken box [0.02,0.98], and in the far
    #   corner of the wedge E just above the switch line the whole partner region has
    #   p1 > 0.98, so the library legitimately misses it.  This is the refinement step:
    #   at every game a witness is CONSTRUCTED from the closed form, checked to be
    #   efficient by m1lib's own predicate, and its Nash slack against all 16 deviators
    #   is evaluated in exact Fractions.
    half = Fraction(1, 2)

    def witness(uu, vv):
        s_ = uu + vv
        if s_ > 1:                                     # Eff_alt, the (0,1) sheet
            if uu >= 0 and vv >= 0:
                return (half, Fraction(0), Fraction(1), half)
            if uu < 0:
                return (half, Fraction(0), Fraction(1), (1 + uu + vv) / (2 * (1 + vv - uu)))
            beta = -2 * vv / (uu - vv - 1)
            return ((1 + beta) / 2, Fraction(0), Fraction(1), half)
        if uu <= 1 and vv <= 0:                        # Eff_CC, everything is a partner
            return (Fraction(1), half, half, half)
        if uu > 1:                                     # only  A (u-1) <= p3 (-v)  binds
            A_ = min(half, -vv / (2 * (uu - 1)))
            return (Fraction(1), 1 - A_, Fraction(1), half)
        return (Fraction(1), Fraction(0),              # v > 0: A >= v max(p3/(1-u), p4)
                min(half, (1 - uu) / (2 * vv)), min(half, 1 / (2 * vv)))

    grid = [Fraction(k, 2) for k in range(-16, 17)]
    nfail, nline, worst = 0, 0, Fraction(0)
    for uu in grid:
        for vv in grid + [1 - uu]:                     # the switch line is included
            p = witness(uu, vv)
            pf = [float(x) for x in p]
            eff = (in_eff_cc(pf) if uu + vv < 1 else
                   in_eff_alt(pf) if uu + vv > 1 else in_eff_line(pf))
            gap = min(c - a_ * (uu - 1) - b_ * vv for a_, b_, c in facets_exact_limit(p))
            if uu + vv == 1:
                nline += 1
            if not eff or gap < -Fraction(1, 10 ** 9):
                nfail += 1
                print("      FAIL at (%s,%s): eff %s gap %.3e" % (uu, vv, eff, float(gap)))
            worst = min(worst, gap)
    print("      %d games on the 1/2 grid, %d of them on the switch line; %d failures;"
          " worst exact Nash gap %.2e" % (len(grid) * (len(grid) + 1), nline, nfail,
                                          float(worst)))
    assert nfail == 0

    print("\n4.  the second alternating sheet (p2,p3) = (1,0)")
    hits = [(u, v) for u in np.arange(-8, 8.01, 0.5) for v in np.arange(-8, 8.01, 0.5)
            if u + v > 1 and rawmask(am, u, v).any() and abs(u - v - 1) > 1e-12]
    print("    games strictly above the line and off u-v = 1 where it contributes: %s"
          % (hits or "none"))
    for u in (2.0, 4.0, 6.0):
        print("    on u - v = 1 at u = %.1f : fraction of the (1,0) sheet that is Nash "
              "= %.4f" % (u, float(rawmask(am, u, u - 1.0).mean())))

    print("\n5.  the switch line u + v = 1 itself  (Eff_line)")
    print("    Eff_line contains Eff_CC, which is 3-dimensional, so the partner set ON the"
          " line is 3-dimensional too -- the dimension drops to 2 only strictly above it.")
    for u in (-6.0, -2.0, 0.0, 0.5, 1.0, 3.0, 7.0):
        v = 1.0 - u
        obs = float(rawmask(ab, u, v).mean())
        cf = (1.0 if (u <= 1 and v <= 0) else
              0.5 if v < 0 else (0.5 - v * v / 6.0 if v <= 1 else 1.0 / (3.0 * v)))
        print("    u = %+5.1f, v = %+5.1f : Eff_CC fraction (raw, shrunken box) %.5f, "
              "closed form %.5f%s" % (u, v, obs, cf,
              "   <- the pencil hub (u,v) = (1,0), where the field has no limit"
              if (u, v) == (1.0, 0.0) else ""))


if __name__ == "__main__":
    main()
