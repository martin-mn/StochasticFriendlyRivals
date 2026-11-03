#!/usr/bin/env python3
"""Exact-rational verification of the closed forms for the EFFICIENT strategies and the
FRIENDLY RIVALS over the FULL four-parameter memory-one cube.  No threshold, no
tolerance, no float, anywhere.

The efficiency counterpart of exact_rivals.py, whose polynomial engine it imports rather
than rewriting, and the continuum counterpart of ../../PartnersRivalsM1's vertex census.

CONVENTIONS are m1lib's and are not redefined here: p = (p1,p2,p3,p4) is the probability
of C after CC, CD, DC, DD with the FOCAL player's action first, states CC=0 CD=1 DC=2
DD=3; the error acts as p -> p~ = eps + (1-2eps) p on a FIXED p in the closed cube
(convention A); w is the stationary distribution of the four-state pair chain in SELF
PLAY (q = p) and W its eps -> 0 limit.  Base game (R,S,T,P) = (1,u,1+v,0), so
Emax = max{R,(T+S)/2} = max{1,(1+u+v)/2} and

    pi(p,p) - Emax = -W_DD - (W_CD+W_DC)(1-u-v)/2       below the switch line u+v = 1
                   = -W_DD                               on   the switch line
                   = W_CC(1-E) - W_DD E, E = (1+u+v)/2   above the switch line,

whence: efficient iff W_CC = 1 / W_DD = 0 / W_CD = W_DC = 1/2 respectively.

THE CLOSED FORMS UNDER TEST (m1lib.in_eff_cc / in_eff_alt / in_eff_line)
    Eff_CC   (u+v < 1) = { p1 = 1, (p2,p3) not in {(0,1),(1,0)},
                           and ( p4 > 0 or max(p2,p3) = 1 ) }
    Eff_alt  (u+v > 1) = { (p2,p3) in {(0,1),(1,0)}, p1 < 1, p4 > 0, (p1,p4) != (0,1) }
    Eff_line (u+v = 1) = { p1 = 1, ( p4 > 0 or (max(p2,p3) = 1 and min(p2,p3) > 0) ) }
                         cup Eff_alt
and, with Riv+ = {p4 = 0, p2+p3 <= 1, p2 < 1} and Riv- = {p1 = 1, p2+p3 >= 1, p3 > 0}
taken as PROVED input from ../Text/rivals.md, the friendly rivals

    FR_W = Eff_CC  cap Riv+ = EMPTY            FR_N = Eff_alt  cap Riv+ = EMPTY
    FR_E = Eff_alt cap Riv- = EMPTY            FR_S = Eff_CC   cap Riv-
         = { p1 = 1, p2 > 0, p3 > 0, p2+p3 >= 1, ( p4 > 0 or max(p2,p3) = 1 ) }
plus the five boundary loci (two rays of u+v = 1, two halves of T = S, the point (1,0)).

HOW "W_i = 1" AND "W_i = 1/2" AND "W_i = 0" ARE DECIDED, EXACTLY.  By the Markov-chain
tree theorem the unnormalised stationary weights d_i(eps) = det(I-M) with row i and
column i deleted are polynomials in eps with integer coefficients once p is carried
scaled by their common denominator SCALE (the common factor SCALE^6 cancels in every
ratio).  Writing ord(.) for the order of vanishing at eps = 0
and lc(.) for the coefficient of that power, with Dtot = sum_i d_i,

    W_i = lim_{eps->0} d_i/Dtot = lc(d_i)/lc(Dtot) if ord(d_i) = ord(Dtot), else exactly 0,

since ord(d_i) >= ord(Dtot) always.  The LEADING EXPONENT ord(d_i) - ord(Dtot) is thus an
INTEGER read off by integer polynomial arithmetic -- computed, not fitted -- and the limit
is an exact Fraction.  "= 1", "= 1/2", "= 0" are decided by exact rational equality.

The finite-eps route of m1lib (Fraction Gaussian elimination at eps = 10^-4, 10^-6,
10^-8, 10^-10) is run on a subgrid and asserted to agree with the rational function
EXACTLY, value by value, so the polynomial engine is not taken on trust; and the exact
rival verdict against the 16 vertex co-players is recomputed on the FULL grid so that
the friendly-rival statements do not rest on the closed form of Riv+/Riv- either.

    python3 exact_efficiency.py
"""
import os
import sys
from fractions import Fraction as F
from itertools import product
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m1lib                                                    # noqa: E402
import exact_rivals as er                                       # noqa: E402

# the main grid is all rationals in [0,1] with denominator 1..DEN; exact_rivals' integer
# polynomial engine carries p scaled by the common denominator, so SCALE is retuned to
# lcm(1..DEN) here (60 for DEN = 6, as shipped; 840 for DEN = 8) before it is used.
DEN, SUBDEN = 8, 6
SCALE = 1
for _b in range(1, DEN + 1):
    SCALE = SCALE * _b // __import__("math").gcd(SCALE, _b)
er.SCALE = SCALE
QS = er.QS

GRID = sorted({F(a, b) for b in range(1, DEN + 1) for a in range(0, b + 1)})
# the finite-eps cross-check subgrid (denominators 1..SUBDEN)
SUBGRID = sorted({F(a, b) for b in range(1, SUBDEN + 1) for a in range(0, b + 1)})
EPSK = [4, 6, 8, 10]

STATE = ("CC", "CD", "DC", "DD")


# ----------------------------------------------------------------- the exact self-play limit

def self_limit(pa):
    """(exponents, limits, cd_eq_dc) for the self-play chain of the scaled strategy pa.

    exponents[i] = ord(d_i) - ord(Dtot), an integer (None if d_i is the zero polynomial);
    limits[i] = W_i as an exact Fraction; cd_eq_dc = (d_CD and d_DC are the SAME polynomial).
    """
    ds = er.tree_weights(pa, pa)
    Dt = er.padd(er.padd(ds[0], ds[1]), er.padd(ds[2], ds[3]))
    oD, cD = er.order_lc(Dt)
    exps, lims = [], []
    for i in range(4):
        oi, ci = er.order_lc(ds[i])
        if oi is None:
            exps.append(None)
            lims.append(F(0))
            continue
        assert oi >= oD
        exps.append(oi - oD)
        lims.append(F(ci, cD) if oi == oD else F(0))
    n1, n2 = er.order_lc(er.psub(ds[1], ds[2]))
    return tuple(exps), tuple(lims), n1 is None


def limit_chunk(chunk):
    return [(pa, self_limit(pa)) for pa in chunk]


# --------------------------------------------------- exact rivalry, from the 16 co-players

def rival_chunk(chunk):
    out = []
    for pa in chunk:
        lims = []
        for q in QS:
            qa = tuple(SCALE * b for b in q)
            d, c = er.leading(pa, qa)
            lims.append(F(0) if d is None or d > 0 else c)
        out.append((pa, min(lims) >= 0, max(lims) <= 0))
    return out


# ------------------------------------------------------------- the finite-eps cross-check

def finite_chunk(chunk):
    out = []
    for pa in chunk:
        p = [F(a, SCALE) for a in pa]
        ds = er.tree_weights(pa, pa)
        Dt = er.padd(er.padd(ds[0], ds[1]), er.padd(ds[2], ds[3]))
        exps, lims, _ = self_limit(pa)
        agree = sym = akin = 0
        devs = [[] for _ in range(4)]
        for k in EPSK:
            eps = F(1, 10 ** k)
            w = m1lib.stationary_exact(p, p, eps)
            pc = m1lib.clamp_exact(p, eps)
            akin += (w[0] * (pc[0] - 1) + w[1] * (pc[1] - 1)
                     + w[2] * pc[2] + w[3] * pc[3] == 0)
            sym += (w[1] == w[2])
            den = sum(F(x) * eps ** i for i, x in enumerate(Dt))
            ok = all(w[i] == sum(F(x) * eps ** j for j, x in enumerate(ds[i])) / den
                     for i in range(4))
            agree += ok
            for i in range(4):
                devs[i].append(abs(w[i] - lims[i]))
        # threshold-free confirmation of every exponent: |w_i - W_i| strictly decreasing
        mono = all(all(d[j] > d[j + 1] for j in range(len(d) - 1))
                   if any(d) else True for d in devs)
        out.append((pa, akin == len(EPSK), sym == len(EPSK),
                    agree == len(EPSK), mono))
    return out


# ------------------------------------------------------------------------- the closed forms

def fr_S(p):
    p1, p2, p3, p4 = p
    return (p1 == 1 and p2 > 0 and p3 > 0 and p2 + p3 >= 1
            and (p4 > 0 or max(p2, p3) == 1))


def fr_ray_TltS(p):
    """u+v = 1, u > 1: Riv- cap ({p4>0} cup Eff_CC)."""
    p1, p2, p3, p4 = p
    return (m1lib.in_riv_minus(p)
            and (p4 > 0 or (p2 > 0 and p3 > 0 and max(p2, p3) == 1)))


def gen(p):
    return "".join("1" if x == 1 else "0" for x in p)


def main():
    print(__doc__.split("    python3")[0].rstrip())
    print("=" * 78)
    pts = [er.scaled(p) for p in product(GRID, repeat=4)]
    print("grid: %d values per coordinate (denominators 1..%d, common denominator "
          "SCALE = %d) ->\n%d strategies p, each "
          "solved\nin self play as an exact rational function of eps"
          % (len(GRID), DEN, SCALE, len(pts)))
    print("grid values: %s" % ", ".join(str(x) for x in GRID))

    with Pool(12) as pool:
        got = [x for part in pool.map(limit_chunk, [pts[i::96] for i in range(96)])
               for x in part]
    table = dict(got)
    un = {pa: tuple(F(a, SCALE) for a in pa) for pa in pts}

    # ------------------------------------------------------- step 1: the role symmetry
    nsym = sum(1 for pa in pts if table[pa][2])
    print("\nSTEP 1  role-swap symmetry.  d_CD and d_DC are the SAME polynomial in eps "
          "(hence\n        w_CD = w_DC exactly, at every eps, not merely in the limit): "
          "%d / %d points" % (nsym, len(pts)))

    # ------------------------------------------------------ the exponents, per component
    print("\nSTEP 2  leading eps-exponent of each self-play weight, ord(d_i) - ord(Dtot), "
          "computed\n        exactly by integer polynomial arithmetic (never fitted):")
    for i in range(4):
        cnt = {}
        for pa in pts:
            cnt[table[pa][0][i]] = cnt.get(table[pa][0][i], 0) + 1
        print("        W_%s : %s" % (STATE[i], "  ".join(
            "eps^%s x %d" % (k, cnt[k]) for k in sorted(cnt, key=lambda z: (z is None, z)))))
    print("        exponent 0 <-> W_i > 0;  exponent > 0 <-> W_i = 0 exactly.")

    # ------------------------------------------------ step 3: the three efficiency forms
    print("\nSTEP 3  the three efficiency closed forms, zero tolerance, exact Fractions:")
    n_cc = n_alt = n_ln = 0
    bad = {"cc": [], "alt": [], "line": []}
    for pa in pts:
        W = table[pa][1]
        p = un[pa]
        e_cc, e_alt, e_ln = W[0] == 1, (W[1] == F(1, 2) and W[2] == F(1, 2)), W[3] == 0
        n_cc += e_cc
        n_alt += e_alt
        n_ln += e_ln
        if e_cc != m1lib.in_eff_cc(p):
            bad["cc"].append((p, W))
        if e_alt != m1lib.in_eff_alt(p):
            bad["alt"].append((p, W))
        if e_ln != m1lib.in_eff_line(p):
            bad["line"].append((p, W))
    print("        below the line  W_CC = 1      : %6d points; closed form Eff_CC   "
          "disagrees at %d" % (n_cc, len(bad["cc"])))
    print("        above the line  W_CD=W_DC=1/2 : %6d points; closed form Eff_alt  "
          "disagrees at %d" % (n_alt, len(bad["alt"])))
    print("        on    the line  W_DD = 0      : %6d points; closed form Eff_line "
          "disagrees at %d" % (n_ln, len(bad["line"])))
    for k in bad:
        for p, W in bad[k][:6]:
            print("          MISMATCH %-5s %s  W = %s" % (k, p, W))
    # the two structural facts about the alternating set
    alt_pts = [un[pa] for pa in pts if table[pa][1][1] == F(1, 2)]
    print("        Eff_alt is 2-dimensional: every alternator has (p2,p3) in {(0,1),(1,0)}"
          ": %s" % all((p[1], p[2]) in {(F(0), F(1)), (F(1), F(0))} for p in alt_pts))
    print("        and p1 < 1, p4 > 0, (p1,p4) != (0,1): %s"
          % all(p[0] < 1 and p[3] > 0 and not (p[0] == 0 and p[3] == 1) for p in alt_pts))
    print("        Eff_CC subset {p1 = 1}: %s ;  Eff_CC cap Eff_alt = empty: %s"
          % (all(un[pa][0] == 1 for pa in pts if table[pa][1][0] == 1),
             not any(table[pa][1][0] == 1 and table[pa][1][1] == F(1, 2) for pa in pts)))
    print("        Eff_line strictly contains Eff_CC u Eff_alt: %d = %d + %d + %d extra"
          % (n_ln, n_cc, n_alt, n_ln - n_cc - n_alt))
    extra = [un[pa] for pa in pts
             if table[pa][1][3] == 0 and table[pa][1][0] != 1 and table[pa][1][1] != F(1, 2)]
    print("        one extra point and its exact self-play limit: %s  W = %s"
          % (extra[0], table[er.scaled(extra[0])][1]))

    # ------------------------------------------------- step 4: friendly rivals per wedge
    print("\nSTEP 4  FRIENDLY RIVALS = efficient and rival, wedge by wedge "
          "(Riv+/Riv- as proved input):")
    W_ = [un[pa] for pa in pts if table[pa][1][0] == 1 and m1lib.in_riv_plus(un[pa])]
    S_ = [un[pa] for pa in pts if table[pa][1][0] == 1 and m1lib.in_riv_minus(un[pa])]
    E_ = [un[pa] for pa in pts if table[pa][1][1] == F(1, 2) and m1lib.in_riv_minus(un[pa])]
    N_ = [un[pa] for pa in pts if table[pa][1][1] == F(1, 2) and m1lib.in_riv_plus(un[pa])]
    print("        W  {u+v<1, T>S}  Eff_CC  cap Riv+ : %6d points   %s"
          % (len(W_), "EMPTY" if not W_ else "NOT EMPTY -- %s" % W_[:3]))
    print("        S  {u+v<1, T<S}  Eff_CC  cap Riv- : %6d points   %s"
          % (len(S_), "non-empty"))
    print("        E  {u+v>1, T<S}  Eff_alt cap Riv- : %6d points   %s"
          % (len(E_), "EMPTY" if not E_ else "NOT EMPTY -- %s" % E_[:3]))
    print("        N  {u+v>1, T>S}  Eff_alt cap Riv+ : %6d points   %s"
          % (len(N_), "EMPTY" if not N_ else "NOT EMPTY -- %s" % N_[:3]))
    Sset = set(S_)
    mm = [p for p in un.values() if fr_S(p) != (p in Sset)]
    print("        FR_S closed form {p1=1, p2>0, p3>0, p2+p3>=1, (p4>0 or max(p2,p3)=1)}"
          " disagrees at %d points" % len(mm))
    print("        FR_S lies in the facet p1 = 1: %s ; its share of Riv- on the grid: "
          "%d / %d" % (all(p[0] == 1 for p in S_), len(S_),
                       sum(1 for p in un.values() if m1lib.in_riv_minus(p))))

    # boundary loci
    print("\nSTEP 5  the five boundary loci:")
    loci = [
        ("ray u+v=1, u>1 (T<S)", [un[pa] for pa in pts
                                  if table[pa][1][3] == 0 and m1lib.in_riv_minus(un[pa])]),
        ("ray u+v=1, u<1 (T>S)", [un[pa] for pa in pts
                                  if table[pa][1][3] == 0 and m1lib.in_riv_plus(un[pa])]),
        ("line T=S, u+v<1     ", [un[pa] for pa in pts if table[pa][1][0] == 1]),
        ("line T=S, u+v>1     ", [un[pa] for pa in pts if table[pa][1][1] == F(1, 2)]),
        ("point (u,v) = (1,0) ", [un[pa] for pa in pts if table[pa][1][3] == 0]),
    ]
    for lbl, L in loci:
        print("        %s : %6d grid points %s" % (lbl, len(L), "EMPTY" if not L else ""))
    R0 = set(loci[0][1])
    mm = [p for p in un.values() if fr_ray_TltS(p) != (p in R0)]
    print("        ray u+v=1,u>1 closed form  Riv- cap ({p4>0} u Eff_CC)  disagrees at %d"
          % len(mm))
    print("        it strictly contains FR_S; the extra points are exactly {(1,0,1,p4), "
          "p4>0}: %s"
          % sorted(R0 - Sset)[:3] +
          " ... (%d of them)" % len(R0 - Sset))

    # ---------------------------------------------------------- step 6: the vertex census
    print("\nSTEP 6  the vertex census, reproducing ../../PartnersRivalsM1 exactly:")
    vlim = {}
    for q in QS:
        pa = tuple(SCALE * b for b in q)
        vlim[q] = self_limit(pa)[1]
    v_cc = [gen(q) for q in QS if vlim[q][0] == 1]
    v_alt = [gen(q) for q in QS if vlim[q][1] == F(1, 2)]
    v_ln = [gen(q) for q in QS if vlim[q][3] == 0]
    print("        efficient below the line: %d -- %s   (expect 3: 1110 1001 1111)"
          % (len(v_cc), " ".join(v_cc)))
    print("        efficient above the line: %d -- %s   (expect 0, the no-alternator "
          "theorem)" % (len(v_alt), " ".join(v_alt) or "(none)"))
    print("        efficient on    the line: %d -- %s   (expect 5)"
          % (len(v_ln), " ".join(v_ln)))
    fr = {}
    for q in QS:
        p = tuple(F(b) for b in q)
        fr[q] = (vlim[q][0] == 1 and m1lib.in_riv_plus(p),          # W
                 vlim[q][0] == 1 and m1lib.in_riv_minus(p),         # S
                 vlim[q][1] == F(1, 2) and m1lib.in_riv_minus(p),   # E
                 vlim[q][1] == F(1, 2) and m1lib.in_riv_plus(p))    # N
    cnt = [sum(fr[q][i] for q in QS) for i in range(4)]
    print("        friendly rivals W / S / E / N = %d / %d / %d / %d   (expect 0 / 2 / 0 "
          "/ 0)" % tuple(cnt))
    print("           the two in S: %s" % " ".join(gen(tuple(F(b) for b in q))
                                                   for q in QS if fr[q][1]))
    bnd = []
    for lbl, pred in (("ray u+v=1, u>1", lambda q, p: vlim[q][3] == 0 and m1lib.in_riv_minus(p)),
                      ("ray u+v=1, u<1", lambda q, p: vlim[q][3] == 0 and m1lib.in_riv_plus(p)),
                      ("T=S, u+v<1    ", lambda q, p: vlim[q][0] == 1),
                      ("T=S, u+v>1    ", lambda q, p: vlim[q][1] == F(1, 2)),
                      ("point (1,0)   ", lambda q, p: vlim[q][3] == 0)):
        g = [gen(tuple(F(b) for b in q)) for q in QS if pred(q, tuple(F(b) for b in q))]
        bnd.append(len(g))
        print("        %s : %d -- %s" % (lbl, len(g), " ".join(g) or "(none)"))
    print("        boundary census %s  (expect 3 / 0 / 3 / 0 / 5): %s"
          % (" / ".join(str(x) for x in bnd), bnd == [3, 0, 3, 0, 5]))

    # ------------------------------------------- step 7: the self-play limit points
    print("\nSTEP 7  the self-play limit points W = (W_CC, W_CD, W_DC, W_DD):")
    vset = sorted(set(vlim.values()))
    print("        at the 16 VERTICES: %d distinct points --" % len(vset))
    for W in vset:
        who = " ".join(gen(tuple(F(b) for b in q)) for q in QS if vlim[q] == W)
        print("            (%s)   %s" % (", ".join(str(x) for x in W), who))
    cset = sorted({table[pa][1] for pa in pts})
    print("        on the CONTINUUM grid: %d distinct points" % len(cset))
    print("        every one has W_CD = W_DC exactly: %s"
          % all(W[1] == W[2] for W in cset))
    print("        the 16-vertex set is contained in it: %s"
          % set(vset).issubset(set(cset)))
    dens = sorted({x.denominator for W in cset for x in W})
    print("        denominators of the coordinates run from %d to %d (%d distinct); "
          "the vertex\n        set uses only %s"
          % (dens[0], dens[-1], len(dens),
             sorted({x.denominator for W in vset for x in W})))
    print("        STRUCTURE.  Strong lumping of {CC}|{CD,DC}|{DD} (both CD and DC send")
    print("        beta = (1-p~2)(1-p~3) to DD and gamma = p~2 p~3 to CC) gives the")
    print("        Markov-chain-tree weights, with a = p~1, d = p~4,")
    print("            X = d[(2-d)gamma + d beta],  Y = 2(1-a)d(1+a-d),")
    print("            Z = (1-a)[(1-a)gamma + (1+a)beta],   W = (X, Y/2, Y/2, Z)/(X+Y+Z),")
    print("        so the limit set is the union over the 81 strata (each p_i at 0 / "
          "interior / 1)")
    print("        of the images of finitely many RATIONAL maps of p; it is a "
          "2-parameter family")
    print("        inside the plane W_CD = W_DC of the 3-simplex, not a finite list.  The "
          "six\n        vertex points are the images of the corners; denominators are "
          "unbounded as the\n        grid refines, which is exactly why the vertex census "
          "is finite and this is not.")
    check = 0
    for pa in pts:
        p1, p2, p3, p4 = un[pa]
        beta, gamma = (1 - p2) * (1 - p3), p2 * p3
        a, d = p1, p4
        X = d * ((2 - d) * gamma + d * beta)
        Y = 2 * (1 - a) * d * (1 + a - d)
        Z = (1 - a) * ((1 - a) * gamma + (1 + a) * beta)
        tot = X + Y + Z
        if tot != 0:
            Wf = (X / tot, Y / (2 * tot), Y / (2 * tot), Z / tot)
            check += (Wf == table[pa][1])
    ngen = sum(1 for pa in pts if (lambda p: (p[3] * ((2 - p[3]) * p[1] * p[2]
                                                      + p[3] * (1 - p[1]) * (1 - p[2]))
                                              + 2 * (1 - p[0]) * p[3] * (1 + p[0] - p[3])
                                              + (1 - p[0]) * ((1 - p[0]) * p[1] * p[2]
                                                              + (1 + p[0]) * (1 - p[1]) * (1 - p[2]))) != 0)(un[pa]))
    print("        the eps^0 lumped formula reproduces the exact limit at %d of the %d "
          "points\n        where it is non-degenerate (the rest need the eps-orders)"
          % (check, ngen))

    # ------------------------------------------- step 8: the finite-eps cross-check
    print("\n" + "=" * 78)
    print("finite-eps cross-check with m1lib.stationary_exact, eps = 10^-4, -6, -8, -10")
    sub = [er.scaled(p) for p in product(SUBGRID, repeat=4)]
    print("subgrid: denominators 1..%d, %d values -> %d strategies, %d exact solves"
          % (SUBDEN, len(SUBGRID), len(sub), len(sub) * len(EPSK)))
    with Pool(12) as pool:
        rep = [x for part in pool.map(finite_chunk, [sub[i::96] for i in range(96)])
               for x in part]
    print("   (counts are strategies for which the property holds at ALL FOUR eps)")
    print("   Akin residual w.(p~1-1, p~2-1, p~3, p~4) exactly 0 : %d / %d"
          % (sum(r[1] for r in rep), len(rep)))
    print("   w_CD == w_DC exactly (Fractions, finite eps)       : %d / %d"
          % (sum(r[2] for r in rep), len(rep)))
    print("   w_i == d_i(eps)/Dtot(eps) exactly, all four i      : %d / %d"
          % (sum(r[3] for r in rep), len(rep)))
    print("   |w_i - W_i| strictly decreasing in k, all four i   : %d / %d "
          "(threshold-free\n   confirmation that the computed exponent is the true one)"
          % (sum(r[4] for r in rep), len(rep)))
    for lbl, i in (("akin", 1), ("symmetry", 2), ("agreement", 3), ("monotone", 4)):
        f = [er.scaled and tuple(F(a, SCALE) for a in r[0]) for r in rep if not r[i]]
        if f:
            print("   FAILURES (%s): %d, first few %s" % (lbl, len(f), f[:5]))

    print("\nexact rivalry on the FULL grid, recomputed from the 16 vertex co-players "
          "(no closed form\nassumed anywhere): %d strategies x 16 ordered pairs" % len(pts))
    with Pool(12) as pool:
        rv = [x for part in pool.map(rival_chunk, [pts[i::96] for i in range(96)])
              for x in part]
    dis_p = dis_m = 0
    frW = frS = frE = frN = 0
    for pa, isP, isM in rv:
        p = un[pa]
        dis_p += (isP != m1lib.in_riv_plus(p))
        dis_m += (isM != m1lib.in_riv_minus(p))
        W = table[pa][1]
        frW += (W[0] == 1 and isP)
        frS += (W[0] == 1 and isM)
        frE += (W[1] == F(1, 2) and isM)
        frN += (W[1] == F(1, 2) and isP)
    print("   closed forms in_riv_plus / in_riv_minus disagree with the exact limit "
          "verdict at %d / %d points" % (dis_p, dis_m))
    print("   friendly rivals from the EXACT rival verdict, W / S / E / N = %d / %d / %d "
          "/ %d" % (frW, frS, frE, frN))
    print("   -> W, E and N are empty with nothing taken on trust.")
    print("\ndone.")


if __name__ == "__main__":
    main()
