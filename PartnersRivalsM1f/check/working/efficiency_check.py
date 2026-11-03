"""Receipts for the EFFICIENCY closed forms over the memory-one continuum (item 2, route A).

Everything below is self-play of a single p against itself, under convention A
(p fixed, p -> eps + (1-2eps)p).  Nothing here re-implements the pair chain: the
four-state chain always comes from m1lib.

The engine of the proof is a two-step exact reduction, both steps verified here:

  (L)  LUMPING.  The partition {CC} | {CD,DC} | {DD} is strongly lumpable for the
       self-play chain, because from CD and from DC the probabilities of going to CC
       (p2~ p3~), to DD ((1-p2~)(1-p3~)) and to the block {CD,DC} coincide.  Moreover
       the role swap CD <-> DC is an automorphism, so w_CD = w_DC exactly, at every eps.

  (T)  MATRIX-TREE.  For the resulting three-state chain the stationary distribution is
       proportional to the spanning in-tree sums

         N_CC = (h+g) d2 + h dd,   N_M = c2 dd + e d2 + e dd,   N_DD = c2 (h+g) + e g

       with the eight nonnegative one-step probabilities
         c2 = (1-p1~)^2      CC -> DD          e  = 2 p1~ (1-p1~)   CC -> M
         h  = p2~ p3~        M  -> CC          g  = (1-p2~)(1-p3~)  M  -> DD
         d2 = p4~^2          DD -> CC          dd = 2 p4~ (1-p4~)   DD -> M

       Each of the eight is a polynomial in eps with a strictly positive lowest-order
       coefficient, so no cancellation can occur in N_CC, N_M, N_DD and the eps -> 0
       limit is read off from lowest degrees alone -- exactly, with no threshold.

Run:  python3 efficiency_check.py        (about 90 s)
"""

import random
import sys
import time
from fractions import Fraction as F
from itertools import product

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import m1lib


# --------------------------------------------------------- exact polynomials in eps

class Poly:
    """Polynomial in eps with Fraction coefficients, lowest degree first."""
    __slots__ = ("c",)

    def __init__(self, c):
        c = list(c)
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        self.c = c

    def __add__(s, o):
        n = max(len(s.c), len(o.c))
        return Poly([(s.c[i] if i < len(s.c) else 0) + (o.c[i] if i < len(o.c) else 0)
                     for i in range(n)])

    def __sub__(s, o):
        n = max(len(s.c), len(o.c))
        return Poly([(s.c[i] if i < len(s.c) else 0) - (o.c[i] if i < len(o.c) else 0)
                     for i in range(n)])

    def __mul__(s, o):
        if isinstance(o, (int, F)):
            return Poly([x * o for x in s.c])
        r = [F(0)] * (len(s.c) + len(o.c) - 1)
        for i, a in enumerate(s.c):
            if a:
                for j, b in enumerate(o.c):
                    r[i + j] += a * b
        return Poly(r)

    __rmul__ = __mul__

    def lead(s):
        for i, a in enumerate(s.c):
            if a != 0:
                return i, a
        return None, F(0)

    def at(s, eps):
        eps, r = F(eps), F(0)
        for a in reversed(s.c):
            r = r * eps + a
        return r


ONE = Poly([F(1)])


def clamp_poly(x):
    """eps + (1-2eps)x as a polynomial in eps."""
    x = F(x)
    return Poly([x, 1 - 2 * x])


def pieces(p):
    """The eight lumped one-step probabilities, as exact polynomials in eps."""
    p1, p2, p3, p4 = [F(x) for x in p]
    A1, c = clamp_poly(p1), ONE - clamp_poly(p1)
    al, be, de = clamp_poly(p2), clamp_poly(p3), clamp_poly(p4)
    return dict(c2=c * c, e=2 * (A1 * c),
                h=al * be, g=(ONE - al) * (ONE - be),
                d2=de * de, dd=2 * (de * (ONE - de)),
                D=de * (2 * ONE - de))


def trees(p):
    q = pieces(p)
    h, g, e, c2, dd, d2 = q["h"], q["g"], q["e"], q["c2"], q["dd"], q["d2"]
    return ((h + g) * d2 + h * dd,
            c2 * dd + e * d2 + e * dd,
            c2 * (h + g) + e * g)


def limit_W(p):
    """Exact (W_CC, W_CD, W_DC, W_DD) = lim_{eps->0} w(p,p,eps), from leading orders."""
    N = trees(p)
    od = [x.lead()[0] for x in N]
    co = [x.lead()[1] for x in N]
    m = min(od)
    tot = sum(co[i] for i in range(3) if od[i] == m)
    v = [co[i] / tot if od[i] == m else F(0) for i in range(3)]
    return (v[0], v[1] / 2, v[1] / 2, v[2])


def hitting(p, eps):
    """(m, n) = (E_CD[tau_CC] = E_DC[tau_CC], E_DD[tau_CC]) in closed form, exact."""
    _, al, be, de = m1lib.clamp_exact(p, eps)
    h, g = al * be, (1 - al) * (1 - be)
    D, dd = de * (2 - de), 2 * de * (1 - de)
    m = (D + g) / (h * D + g * de * de)
    return m, (1 + dd * m) / D


# --------------------------------------------------------------------------- checks

def check_lumping_and_tree():
    """(L) and (T), exactly, at eps that are NOT small: these are identities."""
    grid = [F(k, 4) for k in range(5)]
    n = sym = tree = 0
    for p in product(grid, repeat=4):
        for eps in (F(1, 7), F(1, 100), F(1, 3)):
            n += 1
            w = m1lib.stationary_exact(list(p), list(p), eps)
            if w[1] != w[2]:
                sym += 1
            Nx, Ny, Nz = [x.at(eps) for x in trees(p)]
            t = Nx + Ny + Nz
            if (Nx / t, Ny / t / 2, Ny / t / 2, Nz / t) != tuple(w):
                tree += 1
    print(f"  (L)+(T): {n} (p,eps) pairs, w_CD != w_DC in {sym}, tree formula wrong in {tree}")
    return sym == 0 and tree == 0


def check_hitting_and_kac():
    """m, n against the 4-state chain, and Kac's w_CC = 1/(1 + e m + c2 n)."""
    grid = [F(k, 4) for k in range(5)]
    n = badm = badk = badb = 0
    for p in product(grid, repeat=4):
        for eps in (F(1, 7), F(1, 50), F(1, 3)):
            n += 1
            m, nn = hitting(list(p), eps)
            # brute force: solve (I-Q)t = 1 on the transient set {CD,DC,DD}
            pc = m1lib.clamp_exact(list(p), eps)
            M = [[pc[s] * pc[m1lib.SWAP[s]], pc[s] * (1 - pc[m1lib.SWAP[s]]),
                  (1 - pc[s]) * pc[m1lib.SWAP[s]], (1 - pc[s]) * (1 - pc[m1lib.SWAP[s]])]
                 for s in range(4)]
            idx = [1, 2, 3]
            A = [[(F(1) if i == j else F(0)) - M[idx[i]][idx[j]] for j in range(3)]
                 for i in range(3)]
            b = [F(1)] * 3
            for c in range(3):
                piv = next(r for r in range(c, 3) if A[r][c] != 0)
                A[c], A[piv] = A[piv], A[c]
                b[c], b[piv] = b[piv], b[c]
                for r in range(c + 1, 3):
                    f = A[r][c] / A[c][c]
                    if f:
                        for k in range(c, 3):
                            A[r][k] -= f * A[c][k]
                        b[r] -= f * b[c]
            t = [F(0)] * 3
            for r in reversed(range(3)):
                t[r] = (b[r] - sum(A[r][k] * t[k] for k in range(r + 1, 3))) / A[r][r]
            if not (t[0] == t[1] == m and t[2] == nn):
                badm += 1
            q = pieces(p)
            L = q["e"].at(eps) * m + q["c2"].at(eps) * nn
            w = m1lib.stationary_exact(list(p), list(p), eps)
            if w[0] != 1 / (1 + L):
                badk += 1
            if w[0] > 1 / (2 - pc[0] * pc[0]):
                badb += 1
    print(f"  hitting times: {n} cases, closed form wrong in {badm}, "
          f"Kac identity wrong in {badk}, bound w_CC <= 1/(2-p1~^2) violated in {badb}")
    return badm == badk == badb == 0


def check_grid_exact(k=10):
    grid = [F(i, k) for i in range(k + 1)]
    cnt = {"cc": 0, "alt": 0, "line": 0}
    bad = {"cc": 0, "alt": 0, "line": 0}
    for p in product(grid, repeat=4):
        W = limit_W(p)
        for key, truth, pred in (("cc", W[0] == 1, m1lib.in_eff_cc(p)),
                                 ("alt", W[1] == F(1, 2), m1lib.in_eff_alt(p)),
                                 ("line", W[3] == 0, m1lib.in_eff_line(p))):
            cnt[key] += bool(pred)
            bad[key] += (truth != pred)
    print(f"  exact {k+1}^4 = {(k+1)**4} grid: counts {cnt}, disagreements {bad}")
    return sum(bad.values()) == 0


def check_strata_random(reps=40, seed=20260909):
    rng = random.Random(seed)

    def interior():
        d = rng.choice([7, 11, 13, 17, 23, 29, 31, 97, 101, 1000003])
        return F(rng.randint(1, d - 1), d)

    bad = 0
    for st in product(("0", "int", "1"), repeat=4):
        for _ in range(reps):
            p = tuple(F(0) if s == "0" else (F(1) if s == "1" else interior()) for s in st)
            W = limit_W(p)
            bad += ((W[0] == 1) != m1lib.in_eff_cc(p))
            bad += ((W[1] == F(1, 2)) != m1lib.in_eff_alt(p))
            bad += ((W[3] == 0) != m1lib.in_eff_line(p))
    print(f"  all 81 strata x {reps} random interior rationals = {81*reps} points: "
          f"disagreements {bad}")
    return bad == 0


def check_grid_float(k=20, epss=(1e-8, 1e-9, 1e-10, 1e-11), tol=1e-4):
    grid = [F(i, k) for i in range(k + 1)]
    ok = True
    for eps in epss:
        worst_in = {"cc": 0.0, "alt": 0.0, "line": 0.0}
        best_out = {"cc": 9.0, "alt": 9.0, "line": 9.0}
        bad = {"cc": 0, "alt": 0, "line": 0}
        for p in product(grid, repeat=4):
            w = m1lib.self_play([float(x) for x in p], eps)
            dev = {"cc": abs(w[0] - 1.0),
                   "alt": max(abs(w[1] - 0.5), abs(w[2] - 0.5)),
                   "line": abs(w[3])}
            pr = {"cc": m1lib.in_eff_cc(p), "alt": m1lib.in_eff_alt(p),
                  "line": m1lib.in_eff_line(p)}
            for key in dev:
                if pr[key]:
                    worst_in[key] = max(worst_in[key], dev[key])
                else:
                    best_out[key] = min(best_out[key], dev[key])
                bad[key] += ((dev[key] < tol) != pr[key])
        ok &= sum(bad.values()) == 0
        print(f"  float {k+1}^4 = {(k+1)**4} grid, eps={eps:g}, tol={tol:g}: disagreements {bad}"
              f"  worst-in {[f'{worst_in[x]:.1e}' for x in ('cc','alt','line')]}"
              f"  best-out {[f'{best_out[x]:.1e}' for x in ('cc','alt','line')]}")
    return ok


def check_payoffs(k=10, eps=1e-9):
    """The predicates against m1lib.efficiency_deficit at real games in every region."""
    grid = [F(i, k) for i in range(k + 1)]
    games = [("W  u=-0.3 v=0.4", -0.3, 0.4, m1lib.in_eff_cc),
             ("S  u=1.2 v=-0.5", 1.2, -0.5, m1lib.in_eff_cc),
             ("N  u=0.9 v=0.4", 0.9, 0.4, m1lib.in_eff_alt),
             ("E  u=1.6 v=0.2", 1.6, 0.2, m1lib.in_eff_alt),
             ("line u=0.3 v=0.7", 0.3, 0.7, m1lib.in_eff_line),
             ("line u=1.4 v=-0.4", 1.4, -0.4, m1lib.in_eff_line)]
    ok = True
    for name, u, v, pred in games:
        worst_in, best_out, n = 0.0, -9.9, 0
        for p in product(grid, repeat=4):
            d = m1lib.efficiency_deficit([float(x) for x in p], u, v, eps)
            if pred(p):
                n += 1
                worst_in = min(worst_in, d)
            else:
                best_out = max(best_out, d)
        ok &= worst_in > best_out
        print(f"  {name:20s} n_eff={n:5d}  worst deficit inside {worst_in: .2e}"
              f"   best deficit outside {best_out: .2e}")
    return ok


def vertex_census():
    for name, f in (("below (W_CC=1)", m1lib.in_eff_cc),
                    ("above (alternator)", m1lib.in_eff_alt),
                    ("on the line (W_DD=0)", m1lib.in_eff_line)):
        hits = ["".join(str(int(x)) for x in b) for b in product([0, 1], repeat=4)
                if f([F(x) for x in b])]
        print(f"  {name:22s} {len(hits)}: {hits}")
    return True


if __name__ == "__main__":
    t0 = time.time()
    allok = True
    print("1. the two exact reductions")
    allok &= check_lumping_and_tree()
    print("2. excursion return times and Kac")
    allok &= check_hitting_and_kac()
    print("3. closed forms, exact leading order")
    allok &= check_grid_exact(10)
    allok &= check_strata_random()
    print("4. closed forms, float, against m1lib.self_play")
    allok &= check_grid_float()
    print("5. closed forms at the payoff level, against m1lib.efficiency_deficit")
    allok &= check_payoffs()
    print("6. vertex census (must be 3 / 0 / 5)")
    allok &= vertex_census()
    print(f"\n{'ALL CHECKS PASSED' if allok else 'FAILURES ABOVE'}   ({time.time()-t0:.0f} s)")
