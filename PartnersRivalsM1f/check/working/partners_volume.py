#!/usr/bin/env python3
"""ROUTE D: the measure of the partner set at every game, in closed form, and the
refinement of every boundary of the field by bisection.

TWO INDEPENDENT ROUTES TO THE SAME NUMBER.

  * QUADRATURE.  Each of the six raw inequalities below the line is AFFINE IN A = 1 - p2
    with zero intercept, so for a fixed (p3,p4) the feasible A-set is an exactly computable
    interval; above the line each of the seven is affine in p4, so for a fixed p1 the
    feasible p4-set is an interval.  Integrating the interval length on a grid is
    deterministic -- no sampling -- and Richardson-refines in the grid step.
  * CLOSED FORM.  The six inequalities collapse, in each of three regimes, to one or two:

    BELOW  u+v < 1,  Eff_CC = {p1 = 1},  A = 1 - p2, coordinates (p2,p3,p4) in [0,1]^3
      u <= 1, v <= 0 :  every efficient strategy is a partner        vol = 1
      u > 1          :  A (u-1) <= p3 (-v)                           vol = 1 + (u-1)/(2v)
      v > 0          :  A >= v max( p3/(1-u), p4 )                   vol = J / (m v)
                        with m = v/(1-u), M = min(m,1), V = min(v,1) and
                        J = M V - M V^2/2 - M^3/6   (M <= V)
                          = M V - V M^2/2 - V^3/6   (M >= V)

    ABOVE  u+v > 1,  Eff_alt = two sheets (p2,p3) in {(0,1),(1,0)}, coordinates (p1,p4)
      the (1,0) sheet is never a partner sheet off the ray u - v = 1
      u >= 0, v >= 0 :  the whole (0,1) sheet                        area = 1
      u < 0          :  p4 <= (1+u+v)/(1+v-u)                        area = (1+u+v)/(1+v-u)
      v < 0          :  p1 >= -2v/(u-v-1)                            area = (u+v-1)/(u-v-1)

  THE PROOF that the six collapse to two is one sign.  Below the line with v > 0 every one
  of the six is a LOWER bound on A; writing tau = 1-u the six bounds are
     L1 = p4 v,  L2 = p3 v/tau,  L3 = p4 v/(1-p3+p4 tau),  L4 = K v/(1-p3+tau),
     L5 = K v/(1+(1-p4) tau),  L6 = p3 v/(p3+(1-p4) tau),   K = p3+p4-p3 p4,
  and each of L3..L6 <= max(L1,L2) reduces, after one cancellation, to the single
  comparison p3 vs p4 tau -- which is also the comparison that decides which of L1, L2 is
  the larger.  With v < 0 (so u > 1) the bound II implies the other five outright.
  Above the line the same happens: the worst corner of the (p1,p4) square makes C', D',
  E', F' hold with slack exactly (u+v-1) or (u+v-1)(u+v+1)/(...) > 0.
"""
import numpy as np

# --------------------------------------------------------------------- closed forms

def vol_below(u, v):
    """Relative 3-volume of the partner set inside Eff_CC (u+v <= 1)."""
    u, v = np.asarray(u, float), np.asarray(v, float)
    tau = 1.0 - u
    out = np.ones(np.broadcast(u, v).shape)
    out = np.where(u > 1.0, 1.0 + (u - 1.0) / (2.0 * np.where(v == 0, 1.0, v)), out)
    pos = v > 0.0
    m = np.where(pos & (tau > 0), v / np.where(tau > 0, tau, 1.0), 0.0)
    M, V = np.minimum(m, 1.0), np.minimum(v, 1.0)
    J = np.where(M <= V, M * V - M * V ** 2 / 2 - M ** 3 / 6,
                 M * V - V * M ** 2 / 2 - V ** 3 / 6)
    out = np.where(pos, J / np.where(m * v > 0, m * v, 1.0), out)
    return out


def area_above(u, v):
    """Relative 2-area of the partner set inside ONE alternating sheet (u+v >= 1)."""
    u, v = np.asarray(u, float), np.asarray(v, float)
    out = np.ones(np.broadcast(u, v).shape)
    d1 = np.where(np.abs(1 + v - u) > 0, 1 + v - u, 1.0)
    d2 = np.where(np.abs(u - v - 1) > 0, u - v - 1, 1.0)
    out = np.where(u < 0, (1 + u + v) / d1, out)
    out = np.where(v < 0, (u + v - 1) / d2, out)
    return out


def frac_of_efficient(u, v):
    """measure(partner) / measure(efficient).  Eff_CC has 3-volume 1; Eff_alt is TWO unit
    sheets, so the fraction above the line is half the area of the one live sheet."""
    s = u + v
    return np.where(s > 1, area_above(u, v) / 2.0, vol_below(u, v))


def level(u, v):
    """The discrete field: 0 none (never), 1 dim 2 partial, 2 dim 2 saturated,
    3 dim 3 partial, 4 dim 3 saturated."""
    u, v = np.asarray(u, float), np.asarray(v, float)
    above = u + v > 1
    sat_a = (u >= 0) & (v >= 0)
    sat_b = (u <= 1) & (v <= 0)
    return np.where(above, np.where(sat_a, 2, 1), np.where(sat_b, 4, 3))


# ---------------------------------------------------------- deterministic quadrature

def _interval(alpha, beta, lo, hi):
    """Intersect {A : alpha A >= beta} with [lo,hi], elementwise."""
    pos, neg = alpha > 0, alpha < 0
    zer = ~(pos | neg)
    lo = np.where(pos, np.maximum(lo, beta / np.where(pos, alpha, 1.0)), lo)
    hi = np.where(neg, np.minimum(hi, beta / np.where(neg, alpha, 1.0)), hi)
    hi = np.where(zer & (beta > 0), lo - 1.0, hi)
    return lo, hi


def vol_below_quad(u, v, n=2000):
    """Deterministic midpoint quadrature of the partner volume from the SIX raw
    inequalities, each of which is affine in A with zero intercept."""
    g = (np.arange(n) + 0.5) / n
    p3, p4 = np.meshgrid(g, g, indexing="ij")
    K = p3 + p4 - p3 * p4
    lo = np.zeros_like(p3)
    hi = np.ones_like(p3)
    for alpha, beta in [
            (np.ones_like(p3), p4 * v),                                   # I
            (-(u - 1) * np.ones_like(p3), p3 * v),                        # II
            ((1 - p3) - p4 * (u - 1), p4 * v),                            # III
            ((1 - p3) - (u - 1), K * v),                                  # IV
            (1 - (1 - p4) * (u - 1), K * v),                              # V
            (p3 - (1 - p4) * (u - 1), p3 * v)]:                           # VI
        lo, hi = _interval(alpha, beta, lo, hi)
    return float(np.maximum(hi - lo, 0.0).mean())


def area_above_quad(u, v, n=2000):
    """The same for the (0,1) sheet: seven inequalities, each affine in p4."""
    p1 = (np.arange(n) + 0.5) / n
    one = np.ones_like(p1)
    rows = [(-(1 + v - u) * one, -(1 + u + v) * one),                        # A'
            (0 * one, -(2 * v - p1 * (1 - u + v))),                          # B'
            ((u - 1) * one, -(u + (2 - p1) * v + 1 - p1)),                   # C'
            (p1 * (u - v - 1) + (u + v - 1), -(u + v + 1) * one),            # D'
            ((u - v - 1) * one, -(p1 * (u - v - 1) + (u + 3 * v + 1))),      # E'
            ((1 - p1) * (u - v - 1), -(-2 * p1 * (v + 1) + (u + 3 * v + 1))),  # F'
            ((u + v - 1) * one, -((1 - p1) * (u + v + 1)))]                  # G'
    lo, hi = np.zeros_like(p1), np.ones_like(p1)
    for alpha, beta in rows:
        lo, hi = _interval(alpha, beta, lo, hi)
    return float(np.maximum(hi - lo, 0.0).mean())


# ------------------------------------------------------------------------ bisection

def bisect(pred, a, b, iters=80):
    """Last point where pred is True, along the segment a -> b (pred(a) True)."""
    a, b = np.array(a, float), np.array(b, float)
    for _ in range(iters):
        m = 0.5 * (a + b)
        if pred(*m):
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def contour_point(target, u, vlo, vhi, iters=200):
    """v with frac_of_efficient(u,v) = target, by bisection in v (field decreasing in v)."""
    for _ in range(iters):
        vm = 0.5 * (vlo + vhi)
        if float(frac_of_efficient(np.array(u), np.array(vm))) > target:
            vlo = vm
        else:
            vhi = vm
    return 0.5 * (vlo + vhi)


if __name__ == "__main__":
    print("1.  closed form vs deterministic quadrature")
    worst_b = worst_a = (0.0, None)
    for u in np.arange(-8, 8.01, 0.5):
        for v in np.arange(-8, 8.01, 0.5):
            if u + v < 1:
                d = abs(vol_below_quad(u, v, 1200) - float(vol_below(u, v)))
                if d > worst_b[0]:
                    worst_b = (d, (u, v))
            elif u + v > 1:
                d = abs(area_above_quad(u, v, 400000) - float(area_above(u, v)))
                if d > worst_a[0]:
                    worst_a = (d, (u, v))
    print("    below: max |quad - closed| = %.3e at %s   (midpoint rule, n = 1200^2)"
          % worst_b)
    print("    above: max |quad - closed| = %.3e at %s   (midpoint rule, n = 4e5)"
          % worst_a)

    print("\n2.  the boundaries of the discrete field, by bisection")
    #   saturation below the line: v = 0 for u < 1, and u = 1 for v < 0
    f = lambda a, b: float(frac_of_efficient(np.array(a), np.array(b)))
    for u0 in (-6.0, -2.0, 0.0, 0.5, 0.9):
        vb = bisect(lambda a, b: f(a, b) >= 1.0, (u0, -3.0), (u0, 0.9), 200)
        print("    saturation edge at u = %+6.2f : v = %+.16e   (exact v = 0)" % (u0, vb[1]))
    for v0 in (-0.5, -2.0, -6.0):
        ub = bisect(lambda a, b: f(a, b) >= 1.0, (-3.0, v0), (7.0, v0), 200)
        print("    saturation edge at v = %+6.2f : u = %+.16e   (exact u = 1)" % (v0, ub[0]))
    for v0 in (2.0, 5.0, 8.0):
        ub = bisect(lambda a, b: f(a, b) >= 0.5, (3.0, v0), (-6.0, v0), 200)
        print("    saturation edge at v = %+6.2f : u = %+.16e   (exact u = 0, above)"
              % (v0, ub[0]))

    print("\n3.  are the level sets of the graded field straight?")
    #   Put five points on one level set by bisection in v, fit the best straight line
    #   through them, report the largest residual.  1e-15 says the level set IS a line;
    #   1e-2 says it is not.
    F = lambda a_, b_: float(frac_of_efficient(np.array(a_), np.array(b_)))

    def contour(u0, target, vhot, vcold, it=200):
        """vhot has F > target, vcold has F < target."""
        assert F(u0, vhot) > target > F(u0, vcold), \
            "bad bracket at u = %g: F(%g) = %g, F(%g) = %g" % (
                u0, vhot, F(u0, vhot), vcold, F(u0, vcold))
        for _ in range(it):
            vm = 0.5 * (vhot + vcold)
            vhot, vcold = (vm, vcold) if F(u0, vm) > target else (vhot, vm)
        return 0.5 * (vhot + vcold)

    E = 1e-9
    tests = [
        ("BELOW  u > 1  (v < 0)", 0.90, [1.1, 1.4, 1.8, 2.1, 2.5],
         lambda u0: (-8.0, 1.0 - u0 - E)),
        ("BELOW  0 < u < 1     ", 0.80, [0.1, 0.3, 0.5, 0.7, 0.9],
         lambda u0: (E, 1.0 - u0 - E)),
        ("BELOW  u < 0  (v > 0)", 0.30, [-1.0, -2.0, -3.5, -5.0, -7.0],
         lambda u0: (E, 1.0 - u0 - E)),
        ("ABOVE  v < 0  (u > 1)", 0.10, [2.5, 4.0, 5.5, 7.0, 8.0],
         lambda u0: (-8.0, 1.0 - u0 + E)),
        ("ABOVE  u < 0  (v > 1)", 0.25, [-1.2, -1.6, -2.0, -2.3, -2.6],
         lambda u0: (8.0, 1.0 - u0 + E)),
    ]
    for tag, target, us, brk in tests:
        P = np.array([[u0, contour(u0, target, *brk(u0))] for u0 in us])
        A = np.column_stack([P[:, 0], np.ones(len(P))])
        coef = np.linalg.lstsq(A, P[:, 1], rcond=None)[0]
        r = np.abs(A @ coef - P[:, 1]).max()
        print("    %-22s level %.2f   best line v = %+.9f u %+.9f   max residual %.2e  %s"
              % (tag, target, coef[0], coef[1], r, "LINE" if r < 1e-12 else "NOT A LINE"))
        print("        %s" % np.array2string(P, precision=9).replace("\n", "\n        "))

    print("\n    the two pencils")
    for tag, target, us, brk, hub in [
            ("BELOW u>1 ", 0.90, [1.1, 1.4, 1.8, 2.1, 2.5],
             lambda u0: (-8.0, 1.0 - u0 - E), (1.0, 0.0)),
            ("ABOVE v<0 ", 0.10, [2.5, 4.0, 5.5, 7.0, 8.0],
             lambda u0: (-8.0, 1.0 - u0 + E), (1.0, 0.0)),
            ("ABOVE u<0 ", 0.25, [-1.2, -1.6, -2.0, -2.3, -2.6],
             lambda u0: (8.0, 1.0 - u0 + E), (0.0, -1.0))]:
        P = np.array([[u0, contour(u0, target, *brk(u0))] for u0 in us])
        sl = (P[:, 1] - hub[1]) / (P[:, 0] - hub[0])
        print("      %s level %.2f : slope seen from %s = %s   spread %.2e"
              % (tag, target, hub, np.array2string(sl, precision=12), sl.max() - sl.min()))

    print("\n    the curved branch (below the line, v > 0), identified exactly:")
    print("      0 < u < 1      :  v (3 + t^2) = 6 t (1-c)                t = 1-u")
    print("      u <= 0, v <= 1 :  v (1 + 3 t^2) = 6 t^2 (1-c)")
    print("      u <= 0, v > 1  :  v^2 + 6 c t^2 v - 3 t^2 = 0")
    for u0, cs in ((0.9, (0.90, 0.70, 0.55)), (0.5, (0.90, 0.70, 0.55)),
                   (0.1, (0.90, 0.70, 0.55)), (-0.5, (0.80, 0.50, 0.30)),
                   (-2.0, (0.80, 0.50, 0.30)), (-5.0, (0.80, 0.50, 0.30)),
                   (-7.0, (0.80, 0.50, 0.30))):
        for c in cs:
            vv = contour(u0, c, E, 1.0 - u0 - E)
            t = 1.0 - u0
            if u0 > 0:
                pred = 6 * t * (1 - c) / (3 + t * t)
            elif vv <= 1.0:
                pred = 6 * t * t * (1 - c) / (1 + 3 * t * t)
            else:
                pred = 3 * t * t * (-c + np.sqrt(c * c + 1.0 / (3 * t * t)))
            print("      u = %+5.2f  c = %.2f :  v = %.15f   curve %.15f   |diff| %.2e"
                  % (u0, c, vv, pred, abs(vv - pred)))

    print("\n4.  shares of the square -8 <= u,v <= 8")
    N = 6401
    g = (np.arange(N) + 0.5) / N * 16.0 - 8.0
    U, V = np.meshgrid(g, g, indexing="ij")
    lv = level(U, V)
    for k, nm in ((0, "no partner"), (1, "dim 2, partial"), (2, "dim 2, whole sheet"),
                  (3, "dim 3, partial"), (4, "dim 3, all of Eff_CC")):
        print("    %-22s %8.4f %%" % (nm, 100 * np.mean(lv == k)))
    print("    a partner exists       %8.4f %%" % (100 * np.mean(lv > 0)))
    print("    saturated (every efficient strategy is a partner) %8.4f %%   (= 135.5/256)"
          % (100 * np.mean((lv == 2) | (lv == 4))))
    print("    mean fraction of the efficient set that is a partner: %.6f"
          % float(frac_of_efficient(U, V).mean()))
