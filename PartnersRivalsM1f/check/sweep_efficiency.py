#!/usr/bin/env python3
"""Exhaustive FLOAT sweep of the closed forms for the EFFICIENT set and the FRIENDLY-RIVAL
set over the full four-parameter memory-one cube.  Companion to sweep_rivals.py, which does
the same job for Riv+ / Riv-; the rival closed forms are consumed here as a proved input.

CLOSED FORMS UNDER TEST (conventions in m1lib.py; efficiency depends on the game only
through sign(u+v-1), rivalry only through sign(T-S) = sign(1+v-u)).

    Eff_CC   (u+v < 1) = { p1 = 1, (p2,p3) not in {(0,1),(1,0)}, (p4 > 0 or max(p2,p3) = 1) }
    Eff_alt  (u+v > 1) = { (p2,p3) in {(0,1),(1,0)}, p1 < 1, p4 > 0, (p1,p4) != (0,1) }
    Eff_line (u+v = 1) = { p1 = 1, (p4 > 0 or (max(p2,p3) = 1 and min(p2,p3) > 0)) } u Eff_alt

    Riv+ (T > S) = { p4 = 0, p2+p3 <= 1, p2 < 1 }     Riv- (T < S) = { p1 = 1, p2+p3 >= 1, p3 > 0 }

    FRIENDLY RIVAL = efficient and rival.  Predicted: EMPTY in W, N, E; in S equal to
    Eff_CC n Riv- = { p1=1, p2>0, p3>0, p2+p3>=1, (p4>0 or max(p2,p3)=1) }; on the five
    boundary loci as printed in section [3].

WHAT IS MEASURED.  The self-play limit W(p,p) is computed at eps = 1e-10 and 1e-12 and each
p is classified as all-C (W_CC = 1), perfect alternation (W_CD = W_DC = 1/2), W_DD = 0, or
none of these; that verdict is compared with the closed forms evaluated EXACTLY (on integer
numerators for grid points, on exact 0/1/interior floats for stratified random points).
Friendly rivalhood is then compared against the direct payoff computation: efficiency
deficit pi(p,p) - max(1,(1+u+v)/2) and rival margin min_q sign(T-S)(w_DC - w_CD).

NUMERICS.  The stationary vector is NOT obtained from a linear solve -- at eps = 1e-12 the
pair chain is nearly reducible and an LU/lstsq solve loses ~12 digits.  We reuse the
Markov-chain-tree route of sweep_rivals.py: every unnormalised weight is a sum of products of
non-negative transition entries, so it carries full relative accuracy and no cancellation.
Section [0] is the receipt against exact Fractions.

Runtime ~5 min on 12 cores.
"""

import os
import sys
import time
from itertools import product
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m1lib import (DET, SWAP, clamp, genome, NAMES, in_eff_alt, in_eff_cc,   # noqa: E402
                   in_eff_line, in_riv_minus, in_riv_plus, stationary_exact)
from sweep_rivals import stationary_tree                                      # noqa: E402

EPSILONS = (1e-10, 1e-12)
TOL = 1e-6                 # verdict threshold on W_CC / W_DD / |W_CD - 1/2| and on margins
GRID_N = 20                # (GRID_N+1)^4 = 194481 strategies
N_RANDOM = 200_000         # stratified random points, spread over the 81 strata
CHUNK = 4096

QMAT = np.array([np.asarray(q, float) for q in DET])          # (16,4) the co-players


# ------------------------------------------------------------------ vectorised pair chain

def chain_self(P, eps):
    """(N,4,4) transition matrices of the pair chain in SELF-PLAY, focal P[n] vs itself."""
    a = clamp(P, eps)
    b = a[:, SWAP]
    M = np.empty(P.shape + (4,))
    M[..., 0] = a * b
    M[..., 1] = a * (1 - b)
    M[..., 2] = (1 - a) * b
    M[..., 3] = (1 - a) * (1 - b)
    return M


def chain_vs_vertices(P, eps):
    """(N,16,4,4): focal P[n] against the k-th deterministic co-player."""
    a = clamp(P, eps)[:, None, :]
    b = clamp(QMAT, eps)[:, SWAP][None, :, :]
    M = np.empty(np.broadcast_shapes(a.shape, b.shape) + (4,))
    M[..., 0] = a * b
    M[..., 1] = a * (1 - b)
    M[..., 2] = (1 - a) * b
    M[..., 3] = (1 - a) * (1 - b)
    return M


def self_limit(P, eps):
    """(N,4) self-play stationary distribution."""
    return stationary_tree(chain_self(P, eps))


def rival_margins(P, eps):
    """(mplus, mminus): min over the 16 vertices of +-(w_DC - w_CD)."""
    d = stationary_tree(chain_vs_vertices(P, eps))
    d = d[..., 2] - d[..., 1]
    return d.min(axis=1), -d.max(axis=1)


def payoff_self(w, u, v):
    return w[..., 0] + w[..., 1] * u + w[..., 2] * (1 + v)


def eff_deficit(w, u, v):
    return payoff_self(w, u, v) - max(1.0, (1.0 + u + v) / 2.0)


# ------------------------------------------------------- measured classification of W(p,p)

def classify(w, tol=TOL):
    """(is_allC, is_alt, is_ddzero) from the measured self-play limit."""
    allc = w[..., 0] > 1 - tol
    alt = (np.abs(w[..., 1] - 0.5) < tol) & (np.abs(w[..., 2] - 0.5) < tol)
    ddz = w[..., 3] < tol
    return allc, alt, ddz


# ------------------------------------------- the closed forms, evaluated without rounding

def forms_int(NUM, den):
    """Eff_CC / Eff_alt / Eff_line / Riv+ / Riv- on integer numerators over `den`."""
    i1, i2, i3, i4 = NUM[:, 0], NUM[:, 1], NUM[:, 2], NUM[:, 3]
    corner = ((i2 == 0) & (i3 == den)) | ((i2 == den) & (i3 == 0))
    imax, imin = np.maximum(i2, i3), np.minimum(i2, i3)
    ecc = (i1 == den) & (~corner) & ((i4 > 0) | (imax == den))
    ealt = corner & (i1 < den) & (i4 > 0) & ~((i1 == 0) & (i4 == den))
    eline = ((i1 == den) & ((i4 > 0) | ((imax == den) & (imin > 0)))) | ealt
    rplus = (i4 == 0) & (i2 + i3 <= den) & (i2 < den)
    rminus = (i1 == den) & (i2 + i3 >= den) & (i3 > 0)
    return ecc, ealt, eline, rplus, rminus


def forms_float(P):
    """The same, for points whose coordinates are exact 0/1 or strictly interior floats."""
    p1, p2, p3, p4 = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    corner = ((p2 == 0.0) & (p3 == 1.0)) | ((p2 == 1.0) & (p3 == 0.0))
    pmax, pmin = np.maximum(p2, p3), np.minimum(p2, p3)
    ecc = (p1 == 1.0) & (~corner) & ((p4 > 0.0) | (pmax == 1.0))
    ealt = corner & (p1 < 1.0) & (p4 > 0.0) & ~((p1 == 0.0) & (p4 == 1.0))
    eline = ((p1 == 1.0) & ((p4 > 0.0) | ((pmax == 1.0) & (pmin > 0.0)))) | ealt
    rplus = (p4 == 0.0) & (p2 + p3 <= 1.0) & (p2 < 1.0)
    rminus = (p1 == 1.0) & (p2 + p3 >= 1.0) & (p3 > 0.0)
    return ecc, ealt, eline, rplus, rminus


# ------------------------------------------------------------------------ the nine games
#   (label, u, v, sign of T-S, which efficiency form applies, name of the FR closed form)

GAMES = [
    ("W  wedge  (u+v<1, T>S)          ", -1.0, 1.0, +1, "cc", "empty"),
    ("S  wedge  (u+v<1, T<S)          ", 0.8, -0.5, -1, "cc", "Eff_CC n Riv-"),
    ("E  wedge  (u+v>1, T<S)          ", 1.5, 0.2, -1, "alt", "empty"),
    ("N  wedge  (u+v>1, T>S)          ", 0.2, 1.5, +1, "alt", "empty"),
    ("ray u+v=1, u<1  (W|N, T>S)      ", 0.5, 0.5, +1, "line", "empty"),
    ("ray u+v=1, u>1  (S|E, T<S)      ", 1.5, -0.5, -1, "line", "Eff_line n Riv-"),
    ("line T=S, u+v<1 (W|S, vacuous)  ", 0.2, -0.8, 0, "cc", "Eff_CC"),
    ("line T=S, u+v>1 (N|E, vacuous)  ", 1.5, 0.5, 0, "alt", "Eff_alt"),
    ("point (u,v)=(1,0)  (vacuous)    ", 1.0, 0.0, 0, "line", "Eff_line"),
]


def fr_predicted(ecc, ealt, eline, rplus, rminus, kind, sgn):
    """Closed-form friendly-rival predicate for a game.  sgn = sign(T-S); sgn = 0 is the
    line T = S, where pi(p,q) = pi(q,p) identically and EVERY p is a rival."""
    eff = {"cc": ecc, "alt": ealt, "line": eline}[kind]
    if sgn == 0:
        return eff
    return eff & (rplus if sgn > 0 else rminus)


def measured_margin(mplus, mminus, sgn):
    """The rival margin actually in force at a game.  On T = S the payoff difference
    (T-S)(w_DC - w_CD) vanishes identically, so the margin is 0 whatever w does."""
    if sgn > 0:
        return mplus
    if sgn < 0:
        return mminus
    return np.zeros_like(mplus)


# ------------------------------------------------------------------------ the order test
#
# A disagreement at fixed eps is not yet a refutation: the closed forms describe the eps -> 0
# LIMIT, and the approach to that limit is not uniform in p.  On p = (1,0,0,p4), for instance,
# the exact weights give W_DD ~ 2 eps / p4^2, so any tolerance is breached once p4 <~ sqrt(eps)
# even though the limit is 0.  The test below decides each disagreement by measuring the
# deviation at eps and again at eps * 1e-4:
#
#   closed form says IN  -> the deviation must COLLAPSE (limit 0), factor <= 0.1
#   closed form says OUT -> the deviation must PERSIST  (limit > 0), factor >= 0.5
#
# Anything else is a genuine disagreement and is printed in full.

def deviation(w, cat):
    """Distance from the defining equality of a category; the closed form says it -> 0."""
    if cat == "allC":
        return 1.0 - w[..., 0]
    if cat == "alt":
        return np.maximum(np.abs(w[..., 1] - 0.5), np.abs(w[..., 2] - 0.5))
    return w[..., 3]


def order_test(P, cat, inset, eps):
    """(dev1, dev2, verdict) for a block of disagreeing points."""
    d1 = deviation(self_limit(P, eps), cat)
    d2 = deviation(self_limit(P, eps * 1e-4), cat)
    ok = np.where(inset, d2 <= 0.1 * d1 + 1e-300, d2 >= 0.5 * d1)
    return d1, d2, ok


def report_disagreements(rows, eps, fh=sys.stdout):
    """Print the per-category summary and adjudicate every disagreement by the order test."""
    for nm, want, got, ex, nb in rows:
        print(f"        {nm:>7}: closed form {want:6d}   measured {got:6d}   "
              f"disagreements {nb}", file=fh)
        if not ex:
            continue
        Pb = np.array([e[0] for e in ex])
        inset = np.array([e[2] for e in ex])
        d1, d2, ok = order_test(Pb, nm if nm != "W_DD=0" else "ddz", inset, eps)
        n_art = int(ok.sum())
        print(f"            order test at eps and eps*1e-4: {n_art}/{len(ex)} are finite-eps "
              f"resolution artefacts (the deviation moves the way the closed form says)",
              file=fh)
        for i in range(min(3, len(ex))):
            print(f"              example p={ex[i][0]}  in-set={inset[i]}  "
                  f"dev(eps)={d1[i]:.3e}  dev(eps*1e-4)={d2[i]:.3e}", file=fh)
        genuine = np.nonzero(~ok)[0]
        if len(genuine):
            print(f"            *** {len(genuine)} GENUINE DISAGREEMENTS:", file=fh)
            for i in genuine[:20]:
                print(f"              p={ex[i][0]}  W={ex[i][1]}  form={ex[i][2]} "
                      f"measured={ex[i][3]}  dev {d1[i]:.3e} -> {d2[i]:.3e}", file=fh)
        if nb > len(ex):
            print(f"            (adjudicated the first {len(ex)} of {nb})", file=fh)


# --------------------------------------------------------------------- one block of work

def block_selfplay(args):
    """Classify a block against the closed forms; returns counts and disagreements."""
    P, ecc, ealt, eline, label = args
    out = {}
    for eps in EPSILONS:
        w = self_limit(P, eps)
        allc, alt, ddz = classify(w)
        rows = []
        for name, got, want in (("allC", allc, ecc), ("alt", alt, ealt), ("W_DD=0", ddz, eline)):
            bad = np.nonzero(got != want)[0]
            rows.append((name, int(want.sum()), int(got.sum()),
                         [(tuple(P[i]), tuple(np.round(w[i], 14)), bool(want[i]), bool(got[i]))
                          for i in bad[:400]], len(bad)))
        out[eps] = dict(n=P.shape[0], rows=rows,
                        sym=float(np.abs(w[:, 1] - w[:, 2]).max()))
    return out


def merge_selfplay(parts):
    tot = {}
    for pr in parts:
        for eps, d in pr.items():
            t = tot.setdefault(eps, dict(n=0, rows=None, sym=0.0))
            t["n"] += d["n"]
            t["sym"] = max(t["sym"], d["sym"])
            if t["rows"] is None:
                t["rows"] = [[nm, w, g, list(ex), nb] for nm, w, g, ex, nb in d["rows"]]
            else:
                for r, (nm, w, g, ex, nb) in zip(t["rows"], d["rows"]):
                    r[1] += w
                    r[2] += g
                    r[4] += nb
                    r[3] += ex[:max(0, 400 - len(r[3]))]
    return tot


def block_friendly(args):
    """Friendly-rival predicate versus direct computation, for all nine games."""
    P, forms, eps = args
    ecc, ealt, eline, rplus, rminus = forms
    w = self_limit(P, eps)
    mplus, mminus = rival_margins(P, eps)
    res = []
    for (label, u, v, sgn, kind, _fr) in GAMES:
        dfc = eff_deficit(w, u, v)
        m = measured_margin(mplus, mminus, sgn)
        got = (dfc > -TOL) & (m > -TOL)
        want = fr_predicted(ecc, ealt, eline, rplus, rminus, kind, sgn)
        bad = np.nonzero(got != want)[0]
        res.append((int(want.sum()), int(got.sum()), len(bad),
                    [(tuple(P[i]), float(dfc[i]), float(m[i]), bool(want[i]), bool(got[i]))
                     for i in bad[:400]]))
    # rival-predicate cross-check (the known error-floor edge shows up here, not above)
    rbad_p = int(((mplus > -TOL) != rplus).sum())
    rbad_m = int(((mminus > -TOL) != rminus).sum())
    return res, rbad_p, rbad_m, P.shape[0]


def block_deficit(args):
    """sup over a block of the joint deficit min(efficiency deficit, rival margin)."""
    P, eps, games = args
    w = self_limit(P, eps)
    mplus, mminus = rival_margins(P, eps)
    out = []
    for (label, u, v, sgn) in games:
        d = np.minimum(eff_deficit(w, u, v), mplus if sgn > 0 else mminus)
        i = int(np.argmax(d))
        out.append((float(d[i]), tuple(P[i])))
    return out


# ---------------------------------------------------------------------------- point sets

def grid_points(n):
    idx = np.arange(n + 1, dtype=np.int64)
    NUM = np.stack(np.meshgrid(idx, idx, idx, idx, indexing="ij"), axis=-1).reshape(-1, 4)
    return NUM, NUM.astype(float) / n


def stratified(total, seed=7):
    """Random points from all 81 strata: each coordinate exactly 0, exactly 1, or uniform."""
    rng = np.random.default_rng(seed)
    per = total // 81
    blocks = []
    for pat in product(range(3), repeat=4):
        B = np.empty((per, 4))
        for j, c in enumerate(pat):
            B[:, j] = 0.0 if c == 0 else (1.0 if c == 1 else rng.uniform(1e-12, 1 - 1e-12, per))
        blocks.append(B)
    return np.vstack(blocks)


# -------------------------------------------------------------------------------- driver

def main():
    fh = sys.stdout
    t0 = time.time()
    print("=" * 90)
    print("SWEEP OF THE EFFICIENT SET AND THE FRIENDLY-RIVAL SET OVER THE FULL CUBE")
    print(f"eps = {EPSILONS},  verdict tolerance = {TOL:g},  co-players = the 16 vertices")
    print("=" * 90)

    # ---------------------------------------------------------------- [0] accuracy receipt
    print("\n[0] accuracy receipt: the tree route against exact Fractions (self-play)")
    rng = np.random.default_rng(3)
    tests = [rng.random(4) for _ in range(15)] + \
            [np.array(x, float) for x in ((1, .3, 1, 0), (0, 0, 1, .5), (1, 1, 1, 0),
                                          (.5, 0, 1, .5), (1, .7, .4, .6))]
    for eps in EPSILONS:
        err = 0.0
        for p in tests:
            we = stationary_exact([float(x) for x in p], [float(x) for x in p], eps)
            wf = self_limit(p[None, :], eps)[0]
            err = max(err, max(abs(float(a) - b) for a, b in zip(we, wf)))
        print(f"    eps={eps:.0e}   max |w_tree - w_exact| over 20 self-play points: {err:.3e}")

    pool = Pool(12)

    # ------------------------------------------------------------- [1] the full 21^4 grid
    print(f"\n[1] full grid at resolution 1/{GRID_N}: {(GRID_N+1)**4} strategies")
    NUM, P = grid_points(GRID_N)
    ecc, ealt, eline, rplus, rminus = forms_int(NUM, GRID_N)
    args = [(P[i:i + CHUNK], ecc[i:i + CHUNK], ealt[i:i + CHUNK], eline[i:i + CHUNK], "grid")
            for i in range(0, P.shape[0], CHUNK)]
    tot = merge_selfplay(pool.map(block_selfplay, args))
    for eps in EPSILONS:
        d = tot[eps]
        print(f"    eps={eps:.0e}   n={d['n']}   max|w_CD - w_DC| = {d['sym']:.3e}"
              "   (Step 1: exactly 0 in exact arithmetic)")
        report_disagreements(d["rows"], eps)

    # --------------------------------------------------------- [2] 200000 stratified points
    print(f"\n[2] {81 * (N_RANDOM // 81)} random points, {N_RANDOM // 81} from each of the "
          "81 strata (each coordinate 0, 1, or uniform)")
    R = stratified(N_RANDOM)
    fcc, falt, fline, frp, frm = forms_float(R)
    args = [(R[i:i + CHUNK], fcc[i:i + CHUNK], falt[i:i + CHUNK], fline[i:i + CHUNK], "rand")
            for i in range(0, R.shape[0], CHUNK)]
    tot = merge_selfplay(pool.map(block_selfplay, args))
    for eps in EPSILONS:
        d = tot[eps]
        print(f"    eps={eps:.0e}   n={d['n']}   max|w_CD - w_DC| = {d['sym']:.3e}")
        report_disagreements(d["rows"], eps)

    # ------------------------------------------- [3] friendly rivals, nine games, directly
    eps3 = 1e-10
    NUM2, P2 = grid_points(10)
    T = np.vstack([P2, R[::7]])
    f_int = forms_int(NUM2, 10)
    f_flt = forms_float(R[::7])
    forms = tuple(np.concatenate([a, b]) for a, b in zip(f_int, f_flt))
    print(f"\n[3] friendly-rival closed form vs direct computation at eps={eps3:.0e}, "
          f"on {T.shape[0]} points (11^4 grid + stratified sample)")
    args = [(T[i:i + CHUNK], tuple(f[i:i + CHUNK] for f in forms), eps3)
            for i in range(0, T.shape[0], CHUNK)]
    parts = pool.map(block_friendly, args)
    ng = len(GAMES)
    agg = [[0, 0, 0, []] for _ in range(ng)]
    rbp = rbm = 0
    for res, a, b, _n in parts:
        rbp += a
        rbm += b
        for k, (w_, g_, nb, ex) in enumerate(res):
            agg[k][0] += w_
            agg[k][1] += g_
            agg[k][2] += nb
            agg[k][3] += ex[:max(0, 400 - len(agg[k][3]))]
    for k, (label, u, v, sgn, kind, frname) in enumerate(GAMES):
        w_, g_, nb, ex = agg[k]
        print(f"    {label} u={u:+.2f} v={v:+.2f}  FR = {frname}")
        print(f"        closed form {w_:5d}   direct {g_:5d}   disagreements {nb}")
        if ex:
            # adjudicate exactly as in [1]/[2]: the joint deficit min(effdef, margin) must
            # collapse when the closed form says friendly rival, and persist when it does not
            Pb = np.array([e[0] for e in ex])
            inset = np.array([e[3] for e in ex])
            jd = []
            for e2 in (eps3, eps3 * 1e-4):
                wb = self_limit(Pb, e2)
                mp, mm = rival_margins(Pb, e2)
                jd.append(np.minimum(eff_deficit(wb, u, v),
                                     measured_margin(mp, mm, sgn)))
            ok = np.where(inset, jd[1] >= 0.1 * jd[0], jd[1] <= 0.5 * jd[0])
            print(f"            order test at eps and eps*1e-4: {int(ok.sum())}/{len(ex)} are "
                  f"finite-eps resolution artefacts")
            for i in range(min(3, len(ex))):
                print(f"              example p={ex[i][0]}  in-set={inset[i]}  "
                      f"joint deficit {jd[0][i]:.3e} -> {jd[1][i]:.3e}")
            gen = np.nonzero(~ok)[0]
            if len(gen):
                print(f"            *** {len(gen)} GENUINE DISAGREEMENTS:")
                for i in gen[:20]:
                    print(f"              p={ex[i][0]} effdef={ex[i][1]:.3e} "
                          f"margin={ex[i][2]:.3e} form={ex[i][3]} direct={ex[i][4]}")
    print(f"    rival-predicate cross-check on the same points: Riv+ mismatches {rbp}, "
          f"Riv- mismatches {rbm}")

    # ----------------------------------------------------------------- [4] vertex census
    print("\n[4] vertex census (the 16 deterministic strategies), measured, eps=1e-10")
    V = np.array([np.asarray(q, float) for q in DET])
    w = self_limit(V, 1e-10)
    allc, alt, ddz = classify(w)
    mplus, mminus = rival_margins(V, 1e-10)
    lab = [genome(v) + (" " + NAMES[genome(v)] if genome(v) in NAMES else "") for v in V]
    print("    efficient below the switch line :",
          sorted(lab[i] for i in range(16) if allc[i]))
    print("    efficient above the switch line :",
          sorted(lab[i] for i in range(16) if alt[i]))
    print("    efficient on   the switch line  :",
          sorted(lab[i] for i in range(16) if ddz[i]))
    fcc, falt, fline, frp, frm = forms_float(V)
    print("    closed forms agree:",
          bool((allc == fcc).all() and (alt == falt).all() and (ddz == fline).all()))
    for (label, u, v, sgn, kind, frname) in GAMES:
        dfc = eff_deficit(w, u, v)
        m = measured_margin(mplus, mminus, sgn)
        got = (dfc > -TOL) & (m > -TOL)
        want = fr_predicted(fcc, falt, fline, frp, frm, kind, sgn)
        print(f"    {label} friendly rivals: {sorted(lab[i] for i in range(16) if got[i])}"
              f"   (closed form {sorted(lab[i] for i in range(16) if want[i])})")

    # ------------------------------------------------ [5] sup of the joint deficit by mesh
    print("\n[5] sup over the mesh (hZ)^4 of the joint deficit d = min(efficiency deficit, "
          "rival margin)")
    dgames = [("W", -1.0, 1.0, +1), ("S", 0.8, -0.5, -1),
              ("N", 0.2, 1.5, +1), ("E", 1.5, 0.2, -1)]
    print("    predicted mesh optima:  W -h/2 at (1,h,1,0);  S 0 (attained);  N and E")
    print("    -h/(1+h) at (p1,0,1,h) resp. (1,0,1,h).  A deviation is flagged, not hidden.")
    for n in (5, 10, 20, 40):
        h = 1.0 / n
        pred = {"W": -h / 2, "S": 0.0, "N": -h / (1 + h), "E": -h / (1 + h)}
        _NUM, PG = grid_points(n)
        args = [(PG[i:i + 20000], 1e-10, dgames) for i in range(0, PG.shape[0], 20000)]
        parts = pool.map(block_deficit, args)
        print(f"    h = 1/{n}   ({PG.shape[0]} mesh points)")
        for k, (nm, u, v, sgn) in enumerate(dgames):
            v_, pt = max((pr[k] for pr in parts), key=lambda z: z[0])
            flag = "matches" if abs(v_ - pred[nm]) < 1e-9 else \
                   f"DEVIATES from the predicted {pred[nm]:+.6f}"
            print(f"        {nm}:  sup d = {v_:+.8f}  at p = "
                  f"{tuple(round(x, 4) for x in pt)}   {flag}")

    pool.close()
    print(f"\nTOTAL {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
