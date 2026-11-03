#!/usr/bin/env python3
"""Exact-rational verification of the closed form for the RIVAL sets over the FULL
four-parameter memory-one cube.  No threshold, no tolerance, no float, anywhere.

The continuum counterpart of ../../PartnersRivalsM1/check/m1pairs.py, which did the
same job for the 16 vertices of the cube.

CONVENTIONS are m1lib's and are not redefined here: p = (p1,p2,p3,p4) is the
probability of C after CC, CD, DC, DD with the FOCAL player's action first, states
CC=0 CD=1 DC=2 DD=3; the implementation error acts as p -> p~ = eps + (1-2eps) p on a
FIXED p in the closed cube (convention A); w is the stationary distribution of the
four-state pair chain and W its eps -> 0 limit; pi(p,q) - pi(q,p) = (T-S)(w_DC - w_CD).

THE CLOSED FORM UNDER TEST
    Riv+ (T > S) = { p : p4 = 0, p2 + p3 <= 1, p2 < 1 }
    Riv- (T < S) = { p : p1 = 1, p2 + p3 >= 1, p3 > 0 }
both being the stated prism minus the closed edge E = {(p1,1,0,0)}.

HOW "VANISHES IN THE LIMIT" IS DECIDED, EXACTLY.  Every entry of the transition matrix
is a polynomial in eps with rational coefficients, so by the Markov-chain tree theorem
the unnormalised stationary weights d_i(eps) = det (I - M) with row i and column i
deleted are polynomials in eps with rational coefficients, and

    w_DC(eps) - w_CD(eps) = N(eps) / Dtot(eps),   N = d_DC - d_CD,  Dtot = sum_i d_i,

exactly, as a rational function.  Writing ord(.) for the order of vanishing at eps = 0
and lc(.) for the coefficient of that leading power,

    margin(eps) = (lc N / lc Dtot) * eps^(ord N - ord Dtot) * (1 + O(eps)),

so the LEADING EXPONENT d := ord N - ord Dtot and the LEADING COEFFICIENT
c := lc N / lc Dtot are read off exactly by integer polynomial arithmetic.  The limit is
c when d = 0 and exactly 0 when d > 0 (d < 0 never occurs: |margin| <= 1).  Nothing is
fitted and nothing is compared against a tolerance.  The polynomials are carried in
SCALED INTEGER arithmetic -- every p_i on the grid has denominator dividing 60, and the
common factor 60^6 cancels in N / Dtot -- which is what makes the full grid affordable.

The finite-eps route of m1lib (Fraction Gaussian elimination at eps = 10^-4, 10^-6,
10^-8, 10^-10) is run on a subgrid and asserted to agree with the rational function
EXACTLY, value by value, so the polynomial engine is not taken on trust.

    python3 exact_rivals.py
"""
import os
import sys
from fractions import Fraction as F
from itertools import product
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m1lib                                                    # noqa: E402

SWAP = (0, 2, 1, 3)
SCALE = 60                          # common denominator of every grid value

# the grid: all rationals in [0,1] with denominator 1..6
GRID = sorted({F(a, b) for b in range(1, 7) for a in range(0, b + 1)})
# the finite-eps cross-check subgrid (denominators 1..4), 7^4 = 2401 points
SUBGRID = sorted({F(a, b) for b in range(1, 5) for a in range(0, b + 1)})
EPSK = [4, 6, 8, 10]

QS = [tuple(int(x) for x in q) for q in product([0, 1], repeat=4)]


# ------------------------------------------------------- integer polynomials in eps

def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return r


def padd(a, b):
    if len(a) < len(b):
        a, b = b, a
    r = list(a)
    for i, x in enumerate(b):
        r[i] += x
    return r


def psub(a, b):
    r = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        r[i] += x
    for i, x in enumerate(b):
        r[i] -= x
    return r


def cl(a):
    """SCALE * (eps + (1 - 2 eps) p) for p = a / SCALE."""
    return (a, SCALE - 2 * a)


def ncl(a):
    """SCALE * (1 - clamped p)."""
    return (SCALE - a, -(SCALE - 2 * a))


def det3(M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    return padd(psub(pmul(a, psub(pmul(e, i), pmul(f, h))),
                     pmul(b, psub(pmul(d, i), pmul(f, g)))),
                pmul(c, psub(pmul(d, h), pmul(e, g))))


def tree_weights(pa, qa):
    """Unnormalised stationary weights d_CC, d_CD, d_DC, d_DD as integer polynomials.

    pa, qa are the strategies scaled by SCALE.  The overall factor SCALE^6 common to the
    four is irrelevant: only ratios are ever used.
    """
    SS = SCALE * SCALE
    L = [None] * 4
    for s in range(4):
        a, na = cl(pa[s]), ncl(pa[s])
        b, nb = cl(qa[SWAP[s]]), ncl(qa[SWAP[s]])
        row = (pmul(a, b), pmul(a, nb), pmul(na, b), pmul(na, nb))
        L[s] = [psub([SS] if s == j else [0], row[j]) for j in range(4)]
    out = []
    for k in range(4):
        idx = [i for i in range(4) if i != k]
        out.append(det3([[L[i][j] for j in idx] for i in idx]))
    return out


def order_lc(P):
    """(order of vanishing at eps = 0, coefficient of that power); (None, 0) if P == 0."""
    for i, c in enumerate(P):
        if c:
            return i, c
    return None, 0


def akin_poly(pa, ds):
    """Akin's residual as a polynomial: sum_i d_i (p~_i - 1{i in {CC,CD}}), must be 0."""
    acc = [0]
    for i in range(4):
        t = cl(pa[i])
        if i < 2:
            t = (t[0] - SCALE, t[1])
        acc = padd(acc, pmul(ds[i], t))
    return acc


def margin_poly(pa, qa):
    """((ord, lc) of N, (ord, lc) of Dtot, akin residual polynomial)."""
    ds = tree_weights(pa, qa)
    N = psub(ds[2], ds[1])
    Dt = padd(padd(ds[0], ds[1]), padd(ds[2], ds[3]))
    return order_lc(N), order_lc(Dt), akin_poly(pa, ds)


def leading(pa, qa):
    """(exponent d, coefficient c) of w_DC - w_CD; the limit is c if d == 0 else 0.

    d is None exactly when the margin is the zero rational function.
    """
    (oN, cN), (oD, cD), ak = margin_poly(pa, qa)
    assert not any(ak), "Akin residual is not identically zero"
    if oN is None:
        return None, F(0)
    d = oN - oD
    assert d >= 0, "margin unbounded at eps -> 0"
    return d, F(cN, cD)


# ------------------------------------------------------------------- the grid sweep

def sweep_chunk(chunk):
    """For each scaled p in the chunk: the 16 (d, c) pairs against the vertex q's."""
    res = []
    for pa in chunk:
        res.append((pa, [leading(pa, q) for q in
                         [tuple(SCALE * b for b in q) for q in QS]]))
    return res


def scaled(p):
    return tuple(int(x * SCALE) for x in p)


def limit_of(dc):
    d, c = dc
    return F(0) if d is None or d > 0 else c


def main():
    print(__doc__.split("python3")[0].rstrip())
    print("=" * 78)
    print("grid: %d values per coordinate (denominators 1..6) -> %d strategies p, "
          "each against the 16 deterministic vertices = %d ordered pairs"
          % (len(GRID), len(GRID) ** 4, len(GRID) ** 4 * 16))
    print("grid values: %s" % ", ".join(str(x) for x in GRID))

    pts = [scaled(p) for p in product(GRID, repeat=4)]
    chunks = [pts[i::48] for i in range(48)]
    with Pool(12) as pool:
        got = [x for part in pool.map(sweep_chunk, chunks) for x in part]
    table = dict(got)
    print("\nAkin's residual w . (p~1-1, p~2-1, p~3, p~4) is exactly the zero polynomial "
          "for all %d pairs" % (len(pts) * 16))

    # ---------------------------------------------------------- exponent structure
    exps = {}
    for pa in pts:
        for dc in table[pa]:
            exps[dc[0]] = exps.get(dc[0], 0) + 1
    print("\nleading exponent of eps in w_DC - w_CD, over all pairs:")
    for k in sorted(exps, key=lambda z: (z is None, z)):
        lbl = "identically 0" if k is None else "eps^%d" % k
        print("   %-14s %8d pairs%s" % (lbl, exps[k],
                                        "   (limit non-zero)" if k == 0 else
                                        "   (vanishes in the limit)" if k else ""))

    # -------------------------------------------------------------- rival verdicts
    unscaled = {pa: tuple(F(a, SCALE) for a in pa) for pa in pts}
    bad_p, bad_m = [], []
    riv_p, riv_m = [], []
    for pa in pts:
        lims = [limit_of(dc) for dc in table[pa]]
        isP = min(lims) >= 0
        isM = max(lims) <= 0
        p = unscaled[pa]
        if isP:
            riv_p.append(p)
        if isM:
            riv_m.append(p)
        if isP != m1lib.in_riv_plus(p):
            bad_p.append((p, min(lims)))
        if isM != m1lib.in_riv_minus(p):
            bad_m.append((p, max(lims)))
    print("\nrivals on the grid:  |Riv+| = %d,  |Riv-| = %d,  |Riv+ cap Riv-| = %d "
          "out of %d grid strategies" % (len(riv_p), len(riv_m),
                                         len(set(riv_p) & set(riv_m)), len(pts)))
    print("closed form m1lib.in_riv_plus  disagrees with the exact limit at %d points%s"
          % (len(bad_p), "" if bad_p else "   <-- closed form exact on the grid"))
    for p, m in bad_p[:10]:
        print("     %s  min margin %s" % (str(p), m))
    print("closed form m1lib.in_riv_minus disagrees with the exact limit at %d points%s"
          % (len(bad_m), "" if bad_m else "   <-- closed form exact on the grid"))
    for p, m in bad_m[:10]:
        print("     %s  max margin %s" % (str(p), m))

    inter = sorted(set(riv_p) & set(riv_m))
    ok_line = all(p[0] == 1 and p[3] == 0 and p[1] + p[2] == 1 and p[1] < 1
                  for p in inter)
    print("Riv+ cap Riv- is exactly the generalised-TFT segment (1,s,1-s,0), s < 1: %s "
          "(%d points, s = %s)" % (ok_line, len(inter),
                                   ", ".join(str(p[1]) for p in inter)))

    # ------------------------------------------------------------- the two witnesses
    i_alld, i_0001, i_allc, i_0111 = QS.index((0, 0, 0, 0)), QS.index((0, 0, 0, 1)), \
        QS.index((1, 1, 1, 1)), QS.index((0, 1, 1, 1))
    two_p = [pa for pa in pts if limit_of(table[pa][i_alld]) >= 0
             and limit_of(table[pa][i_0001]) >= 0]
    two_m = [pa for pa in pts if limit_of(table[pa][i_allc]) <= 0
             and limit_of(table[pa][i_0111]) <= 0]
    one_p = [pa for pa in pts if limit_of(table[pa][i_alld]) >= 0]
    print("\nquantifier: {ALLD 0000, 0001} alone reproduce Riv+ exactly: %s   "
          "(ALLD alone gives %d points, %d too many)"
          % (sorted(unscaled[pa] for pa in two_p) == sorted(riv_p),
             len(one_p), len(one_p) - len(riv_p)))
    print("            {ALLC 1111, 0111} alone reproduce Riv- exactly: %s"
          % (sorted(unscaled[pa] for pa in two_m) == sorted(riv_m)))

    # ------------------------------------------------------- the deterministic census
    names = {"0000": "ALLD", "1111": "ALLC", "1010": "TFT", "1000": "Grim"}
    vP, vM = [], []
    for q in QS:
        pa = tuple(SCALE * b for b in q)
        lims = [limit_of(dc) for dc in table[pa]]
        g = "".join(str(b) for b in q)
        if min(lims) >= 0:
            vP.append(g)
        if max(lims) <= 0:
            vM.append(g)
    print("\nvertex census (genome LSB first, 1 = C):")
    print("   rivals T > S: %s" % " ".join(vP))
    print("   rivals T < S: %s" % " ".join(vM))
    print("   both:         %s   %s" % (" ".join(sorted(set(vP) & set(vM))),
                                        ", ".join(names.get(g, g) for g in sorted(set(vP) & set(vM)))))
    print("   reproduces ../../PartnersRivalsM1 exactly: %s"
          % (sorted(vP) == sorted(["0000", "1000", "0010", "1010"])
             and sorted(vM) == sorted(["1010", "1110", "1011", "1111"])))

    # ------------------------------------------------------------- the C/D mirror
    mirror_q = {i: QS.index((1 - q[3], 1 - q[2], 1 - q[1], 1 - q[0]))
                for i, q in enumerate(QS)}
    mirror_ok = all(table[tuple(SCALE - a for a in (pa[3], pa[2], pa[1], pa[0]))][mirror_q[i]]
                    == (table[pa][i][0], -table[pa][i][1])
                    for pa in pts for i in range(16))
    print("\nC/D relabelling m(p) = (1-p4,1-p3,1-p2,1-p1) negates the margin, exactly, "
          "leading term and all: %s" % mirror_ok)

    # --------------------------------- the exact witness values of the theorem, clause 8
    def one(x):
        return F(x, SCALE)

    w_alld = w_alld_E = w_0001 = 0
    n_alld = n_alld_E = n_0001 = 0
    for pa in pts:
        p1, p2, p3, p4 = (one(a) for a in pa)
        if p4 > 0:                                   # ALLD witness, generic branch
            n_alld += 1
            w_alld += (limit_of(table[pa][i_alld]) == -p4 / (1 - p2 + p4))
        if p2 == 1 and p4 == 0:                      # ALLD witness, degenerate corner
            n_alld_E += 1
            w_alld_E += (limit_of(table[pa][i_alld]) == -(1 + p3) / (3 + p3 - p1))
        if p4 == 0 and p2 < 1 and p3 > 0:            # the 0001 witness
            n_0001 += 1
            w_0001 += (limit_of(table[pa][i_0001])
                       == (1 - p2 - p3) / (p3 + 2 * (1 - p2)))
    print("\nexact witness values (theorem clause 8), over the grid:")
    print("   vs ALLD, p4 > 0:            -p4/(1-p2+p4)            %d / %d exact"
          % (w_alld, n_alld))
    print("   vs ALLD, p2 = 1, p4 = 0:    -(1+p3)/(3+p3-p1)        %d / %d exact"
          % (w_alld_E, n_alld_E))
    print("   vs 0001, p4=0, p2<1, p3>0:  (1-p2-p3)/(p3+2(1-p2))   %d / %d exact"
          % (w_0001, n_0001))

    # ------------------------------------------------- the distinct non-zero limits
    vals = sorted({abs(limit_of(dc)) for pa in pts for dc in table[pa]} - {F(0)})
    vert = sorted({abs(limit_of(table[tuple(SCALE * b for b in q)][i]))
                   for q in QS for i in range(16)} - {F(0)})
    seven = [F(1, 5), F(1, 4), F(1, 3), F(2, 5), F(1, 2), F(2, 3), F(1)]
    print("\ndistinct non-zero limit values of |W_DC - W_CD|")
    print("   over the 256 VERTEX pairs: %d values -- %s"
          % (len(vert), ", ".join(str(v) for v in vert)))
    print("   which is the seven of ../../PartnersRivalsM1: %s" % (vert == seven))
    print("   over the whole grid (%d pairs): %d values, smallest %s, largest %s"
          % (len(pts) * 16, len(vals), vals[0], vals[-1]))
    dens = sorted({v.denominator for v in vals})
    print("   %d distinct denominators, from %d to %d; %d values have denominator <= 12"
          % (len(dens), dens[0], dens[-1],
             sum(1 for v in vals if v.denominator <= 12)))
    print("   with denominator <= 6: %s" % ", ".join(str(v)
                                                     for v in vals if v.denominator <= 6))
    print("   the seven vertex values are all present: %s"
          % set(seven).issubset(set(vals)))
    print("   STRUCTURE: the set is not a fixed finite list -- it is the union of the "
          "value sets of\n   finitely many RATIONAL FUNCTIONS of p, one per co-player, "
          "evaluated on the grid.  The\n   two witness families above already generate "
          "%d and %d of the %d values; the rest come\n   from the remaining 14 vertex "
          "co-players.  Denominators are unbounded as the grid\n   refines (they are "
          "denominators of ratios of polynomials in p), which is why the\n   vertex "
          "census is a finite seven-element list and the continuum census is not."
          % (len({abs(v) for v in (-one(pa[3]) / (1 - one(pa[1]) + one(pa[3]))
                                   for pa in pts if pa[3] > 0)} - {F(0)}),
             len({abs((1 - one(pa[1]) - one(pa[2])) / (one(pa[2]) + 2 * (1 - one(pa[1]))))
                  for pa in pts if pa[3] == 0 and one(pa[1]) < 1 and pa[2] > 0} - {F(0)}),
             len(vals)))

    # --------------------------------------- the finite-eps route, exact agreement
    print("\n" + "=" * 78)
    print("finite-eps cross-check with m1lib.stationary_exact, eps = 10^-4, -6, -8, -10")
    print("subgrid: %s  ->  %d strategies, %d pairs, %d exact solves"
          % (", ".join(str(x) for x in SUBGRID), len(SUBGRID) ** 4,
             len(SUBGRID) ** 4 * 16, len(SUBGRID) ** 4 * 16 * len(EPSK)))
    sub = [tuple(SUBGRID[i] for i in idx)
           for idx in product(range(len(SUBGRID)), repeat=4)]
    with Pool(12) as pool:
        rep = pool.map(finite_chunk, [sub[i::48] for i in range(48)])
    n_pairs = n_akin = n_agree = n_sign = n_mono = n_alld = 0
    fail_akin, fail_agree, fail_sign, fail_mono = [], [], [], []
    fixed = {}                       # p -> per-eps list of margins over the 16 q
    for part in rep:
        for (p, q, ak_ok, agree_ok, sign_ok, mono_ok, alld_ok, mm) in part:
            n_pairs += 1
            n_akin += ak_ok
            n_agree += agree_ok
            n_sign += sign_ok
            n_mono += mono_ok
            if q == (0, 0, 0, 0):
                n_alld += alld_ok
            fixed.setdefault(p, [[] for _ in EPSK])
            for i, m in enumerate(mm):
                fixed[p][i].append(m)
            if not ak_ok:
                fail_akin.append((p, q))
            if not agree_ok:
                fail_agree.append((p, q))
            if not sign_ok:
                fail_sign.append((p, q))
            if not mono_ok:
                fail_mono.append((p, q))
    print("   (counts are pairs for which the property holds at ALL FOUR eps)")
    print("   Akin residual exactly 0 (Fractions, clamped p):        %d / %d pairs"
          % (n_akin, n_pairs))
    print("   Gaussian-elimination w_DC - w_CD == N(eps)/Dtot(eps):  %d / %d pairs "
          "(exact equality of Fractions)" % (n_agree, n_pairs))
    print("   sign(margin at eps_k) == sign(leading coefficient):    %d / %d pairs"
          % (n_sign, n_pairs))
    print("   |margin/(c eps^d) - 1| strictly decreasing in k:       %d / %d pairs "
          "(threshold-free confirmation of the exponent)" % (n_mono, n_pairs))
    for lbl, f in (("akin", fail_akin), ("agreement", fail_agree),
                   ("sign", fail_sign), ("monotone", fail_mono)):
        if f:
            print("   FAILURES (%s), first few: %s" % (lbl, f[:5]))
    print("   vs ALLD, margin == -(1-2eps) sum_s w_s p_s exactly:     %d / %d strategies"
          % (n_alld, len(fixed)))
    for i, k in enumerate(EPSK):
        rp = [p for p in fixed if min(fixed[p][i]) >= 0]
        rm = [p for p in fixed if max(fixed[p][i]) <= 0]
        print("   at the FIXED eps = 10^-%-2d :  Riv+(eps) = %s,  Riv-(eps) = %s"
              % (k, [tuple(str(x) for x in p) for p in rp],
                 [tuple(str(x) for x in p) for p in rm]))
    print("   -- rivalry is strictly an eps -> 0 notion (theorem clause 9a)")
    print("\ndone.")


def finite_chunk(chunk):
    out = []
    for p in chunk:
        pa = scaled(p)
        for q in QS:
            qa = tuple(SCALE * b for b in q)
            ds = tree_weights(pa, qa)
            N = psub(ds[2], ds[1])
            Dt = padd(padd(ds[0], ds[1]), padd(ds[2], ds[3]))
            (oN, cN), (oD, cD) = order_lc(N), order_lc(Dt)
            d = None if oN is None else oN - oD
            c = F(0) if oN is None else F(cN, cD)
            ak_ok = agree_ok = sign_ok = alld_ok = 0
            devs = []
            mm = []
            for k in EPSK:
                eps = F(1, 10 ** k)
                w = m1lib.stationary_exact(list(p), [F(b) for b in q], eps)
                pc = m1lib.clamp_exact(list(p), eps)
                ak = w[0] * (pc[0] - 1) + w[1] * (pc[1] - 1) + w[2] * pc[2] + w[3] * pc[3]
                ak_ok += (ak == 0)
                m = w[2] - w[1]
                mm.append(m)
                if q == (0, 0, 0, 0):
                    # theorem clause 9(a): against clamped ALLD the pair distribution
                    # factorises and the margin is exactly -(1-2eps) sum_s w_s p_s
                    alld_ok += (m == -(1 - 2 * eps) * sum(w[s] * F(p[s])
                                                          for s in range(4)))
                num = sum(F(x) * eps ** i for i, x in enumerate(N))
                den = sum(F(x) * eps ** i for i, x in enumerate(Dt))
                agree_ok += (m == num / den)
                if d is None:
                    sign_ok += (m == 0)
                else:
                    sign_ok += ((m > 0) == (c > 0) and (m < 0) == (c < 0))
                    devs.append(abs(m / (c * eps ** d) - 1))
            mono_ok = 1 if d is None else int(
                all(devs[i] > devs[i + 1] for i in range(len(devs) - 1)))
            out.append((p, q, ak_ok == len(EPSK), agree_ok == len(EPSK),
                        sign_ok == len(EPSK), mono_ok, alld_ok == len(EPSK), mm))
    return out


if __name__ == "__main__":
    main()
