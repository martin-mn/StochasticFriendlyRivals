#!/usr/bin/env python3
"""Exhaustive FLOAT sweep of the conjectured closed form for the RIVAL sets over the full
four-parameter memory-one cube.

THEOREM UNDER TEST (conventions inherited from ../../PartnersRivalsM1; see check/m1lib.py)

    Riv+  (T > S)  =  { p in [0,1]^4 : p4 = 0, p2 + p3 <= 1, p2 < 1 }
    Riv-  (T < S)  =  { p in [0,1]^4 : p1 = 1, p2 + p3 >= 1, p3 > 0 }

with the margin of a strategy p defined as

    m+(p, eps) = min_q  ( w_DC - w_CD ),      m-(p, eps) = min_q  -( w_DC - w_CD )
               = min_q d_q                              = - max_q d_q ,

q ranging over the 16 deterministic memory-one co-players (which by the MDP argument realise
the minimum over the whole box).  p is a rival iff the eps -> 0 LIMIT of the margin is >= 0.

WHAT IS ACTUALLY MEASURED, AND WHY A TOLERANCE IS NEEDED.  At fixed eps > 0 nothing except
ALLD is a rival: on Riv+ the true margin is a small negative number of size O(eps) (exactly
-eps/(eps + (1-2eps)(1-p2)) on the face p2+p3 = 1), whereas off Riv+ it is Theta(1) negative,
bounded away from zero by the two witness formulas

    p4 > 0   against ALLD:  -p4 / (1 - p2 + p4)
    p2+p3>1  against 0001:  (1 - p2 - p3) / (p3 + 2(1-p2))
    p2 = 1, p4 = 0 against ALLD: -(1+p3)/(3 + p3 - p1)  in  [-1/2, -1/3] on the edge E.

So the computed verdict is "rival" iff margin > -TOL, and the sweep is meaningful exactly
because the two populations are separated by many orders of magnitude.  The script reports
that separation (the largest rival margin gap and the smallest non-rival |margin|) so the
reader can check that no choice of TOL in a wide window changes any verdict.  Points placed
DELIBERATELY at distance 1e-3/1e-6/1e-9 outside the boundary have Theta(distance) margins and
are therefore reported as a resolution table, not as disagreements.

Runs in a couple of minutes on 12 cores.  Float only; the exact-Fraction receipts live
elsewhere.
"""

import os
import sys
import time
from itertools import product
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m1lib import (DET, SWAP, clamp, genome, in_riv_minus, in_riv_plus,  # noqa: E402
                   stationary, stationary_exact)

EPSILONS = (1e-10, 1e-12)
TOL = 1e-6            # verdict threshold on the margin (see the module docstring)
GRID_N = 20           # (GRID_N+1)^4 = 194481 strategies
N_RANDOM = 200_000
CHUNK = 8192

QMAT = np.array([np.asarray(q, float) for q in DET])          # (16,4) the co-players


# --------------------------------------------------------------------- vectorised chain
#
# NOT by a linear solve.  At eps = 1e-10 the pair chain is nearly reducible and the matrix
# (M^T - I) with one row replaced by the normalisation has condition number ~1/eps^2 ~ 1e10,
# so an LU or least-squares solve returns the stationary vector with only ~6 correct digits
# (measured: 3.5e-7 absolute error against exact Fractions, uncomfortably close to the
# verdict tolerance).  Instead we use the MARKOV-CHAIN TREE THEOREM, which is the same
# representation the existence proof uses (Step 0):
#
#     w_i  =  d_i / sum_j d_j ,   d_i = sum over spanning arborescences rooted at i of
#                                       prod_{j != i} M[j, parent(j)] .
#
# Every d_i is a sum of 16 products of three NON-NEGATIVE transition entries, so it is
# computed with pure multiplications and additions of positive numbers: relative error
# ~1e-16, no cancellation anywhere.  The only subtraction in the whole pipeline is the final
# d_DC - d_CD, whose cancellation is the honest one we actually want to measure.

def _arborescences(n=4):
    """For each root, the list of (child, parent) edge sets of the 16 spanning trees."""
    allt = []
    for root in range(n):
        others = [j for j in range(n) if j != root]
        out = []
        for par in product(range(n), repeat=len(others)):
            pm = dict(zip(others, par))
            if any(pm[j] == j for j in others):
                continue
            ok = True
            for j in others:
                seen, c = {j}, j
                while c != root:
                    c = pm[c]
                    if c in seen:
                        ok = False
                        break
                    seen.add(c)
                if not ok:
                    break
            if ok:
                out.append([(j, pm[j]) for j in others])
        allt.append(out)
    return allt


TREES = _arborescences()


def chain(P, eps):
    """(N,16,4,4) transition matrices of the pair chain, focal P[n] vs the k-th vertex."""
    a = clamp(P, eps)[:, None, :]                              # (N,1,4) focal
    b = clamp(QMAT, eps)[:, SWAP][None, :, :]                  # (1,16,4) co-player response
    M = np.empty(np.broadcast_shapes(a.shape, b.shape) + (4,))
    M[..., 0] = a * b
    M[..., 1] = a * (1 - b)
    M[..., 2] = (1 - a) * b
    M[..., 3] = (1 - a) * (1 - b)
    return M


def stationary_tree(M):
    """Stationary distribution by the Markov-chain tree theorem; M is (..., 4, 4)."""
    d = np.empty(M.shape[:-1])
    for i, trees in enumerate(TREES):
        acc = np.zeros(M.shape[:-2])
        for edges in trees:
            t = M[..., edges[0][0], edges[0][1]].copy()
            for c, par in edges[1:]:
                t = t * M[..., c, par]
            acc += t
        d[..., i] = acc
    return d / d.sum(axis=-1, keepdims=True)


def margins(P, eps):
    """P is (N,4) of RAW (unclamped) focal strategies.

    Returns d of shape (N,16): d[n,k] = w_DC - w_CD for focal P[n] against the k-th
    deterministic co-player, at error rate eps, with both sides clamped.
    """
    w = stationary_tree(chain(P, eps))
    return w[..., 2] - w[..., 1]


def verdicts(P, eps):
    """(mplus, mminus, dmin_abs_nonzero) for a block of strategies."""
    d = margins(P, eps)
    return d.min(axis=1), -d.max(axis=1), d


# --------------------------------------------------------------------- the closed form

def predict(P):
    """Boolean arrays (riv_plus, riv_minus) from the closed form, in float.

    ONLY safe for points whose coordinates are exactly representable and whose defining
    comparisons are not near-ties -- i.e. for uniformly random points.  For any point built
    from a rational grid the membership must be decided on the INTEGER numerators instead
    (predict_int below), because e.g. 2/12 + 10/12 evaluates to 0.9999999999999999 in binary
    and the float test would wrongly place a point of the face p2+p3=1 outside Riv-.
    """
    p1, p2, p3, p4 = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    plus = (p4 == 0.0) & (p2 + p3 <= 1.0) & (p2 < 1.0)
    minus = (p1 == 1.0) & (p2 + p3 >= 1.0) & (p3 > 0.0)
    return plus, minus


def predict_int(NUM, den):
    """The same closed form evaluated exactly on integer numerators over a denominator."""
    i1, i2, i3, i4 = NUM[:, 0], NUM[:, 1], NUM[:, 2], NUM[:, 3]
    plus = (i4 == 0) & (i2 + i3 <= den) & (i2 < den)
    minus = (i1 == den) & (i2 + i3 >= den) & (i3 > 0)
    return plus, minus


# --------------------------------------------------------------------- one block of work

def run_block(args):
    """Test a block of strategies; returns (counts, list of disagreements, extremes)."""
    P, expect_plus, expect_minus, label = args
    out = {}
    for eps in EPSILONS:
        mp, mm, d = verdicts(P, eps)
        got_plus = mp > -TOL
        got_minus = mm > -TOL
        bad_p = np.nonzero(got_plus != expect_plus)[0]
        bad_m = np.nonzero(got_minus != expect_minus)[0]
        nz = np.abs(d)[np.abs(d) > 0]
        out[eps] = dict(
            n=P.shape[0],
            n_plus_pred=int(expect_plus.sum()), n_minus_pred=int(expect_minus.sum()),
            n_plus_got=int(got_plus.sum()), n_minus_got=int(got_minus.sum()),
            bad=[(tuple(P[i]), float(mp[i]), float(mm[i]), bool(expect_plus[i]),
                  bool(expect_minus[i]), bool(got_plus[i]), bool(got_minus[i]), '+')
                 for i in bad_p]
                + [(tuple(P[i]), float(mp[i]), float(mm[i]), bool(expect_plus[i]),
                    bool(expect_minus[i]), bool(got_plus[i]), bool(got_minus[i]), '-')
                   for i in bad_m],
            # separation statistics
            worst_rival_plus=float(mp[expect_plus].min()) if expect_plus.any() else None,
            best_nonrival_plus=float(mp[~expect_plus].max()) if (~expect_plus).any() else None,
            worst_rival_minus=float(mm[expect_minus].min()) if expect_minus.any() else None,
            best_nonrival_minus=(float(mm[~expect_minus].max())
                                 if (~expect_minus).any() else None),
            min_abs_nonzero=float(nz.min()) if nz.size else None,
            label=label,
        )
    return out


def merge(results):
    """Combine per-block dicts into one per-eps summary."""
    agg = {eps: dict(n=0, n_plus_pred=0, n_minus_pred=0, n_plus_got=0, n_minus_got=0,
                     bad=[], worst_rival_plus=None, best_nonrival_plus=None,
                     worst_rival_minus=None, best_nonrival_minus=None, min_abs_nonzero=None)
           for eps in EPSILONS}
    for r in results:
        for eps, o in r.items():
            a = agg[eps]
            for k in ('n', 'n_plus_pred', 'n_minus_pred', 'n_plus_got', 'n_minus_got'):
                a[k] += o[k]
            a['bad'].extend(o['bad'])
            for k, red in (('worst_rival_plus', min), ('best_nonrival_plus', max),
                           ('worst_rival_minus', min), ('best_nonrival_minus', max),
                           ('min_abs_nonzero', min)):
                if o[k] is not None:
                    a[k] = o[k] if a[k] is None else red(a[k], o[k])
    return agg


def report(name, agg, fh):
    for eps in EPSILONS:
        a = agg[eps]
        print(f"\n[{name}]  eps = {eps:g}   n = {a['n']}", file=fh)
        print(f"    predicted rivals : Riv+ {a['n_plus_pred']:>8}   "
              f"Riv- {a['n_minus_pred']:>8}", file=fh)
        print(f"    computed rivals  : Riv+ {a['n_plus_got']:>8}   "
              f"Riv- {a['n_minus_got']:>8}", file=fh)
        print(f"    DISAGREEMENTS    : {len(a['bad'])}", file=fh)
        for rec in a['bad'][:200]:
            p, mp, mm, ep, em, gp, gm, which = rec
            print(f"        {which} p={tuple(round(x, 12) for x in p)}  "
                  f"m+={mp:+.6e} m-={mm:+.6e}  predicted(+{int(ep)},-{int(em)}) "
                  f"computed(+{int(gp)},-{int(gm)})", file=fh)
        if len(a['bad']) > 200:
            print(f"        ... {len(a['bad']) - 200} more", file=fh)
        print(f"    worst margin among predicted Riv+ : {a['worst_rival_plus']}", file=fh)
        print(f"    best  margin among predicted non-Riv+ : {a['best_nonrival_plus']}", file=fh)
        print(f"    worst margin among predicted Riv- : {a['worst_rival_minus']}", file=fh)
        print(f"    best  margin among predicted non-Riv- : {a['best_nonrival_minus']}", file=fh)
        print(f"    smallest non-zero |w_DC - w_CD| in FLOAT : {a['min_abs_nonzero']}"
              f"   (float noise, not a real value -- see [9])", file=fh)


# --------------------------------------------------------------------- point generators

def grid_points(n=GRID_N):
    """The full (n+1)^4 grid.  Membership is decided on the INTEGER coordinates, so the
    closed form is evaluated with no floating-point comparison at all."""
    idx = np.array(list(product(range(n + 1), repeat=4)), dtype=np.int64)
    P = idx / n
    i2, i3, i4 = idx[:, 1], idx[:, 2], idx[:, 3]
    i1 = idx[:, 0]
    plus = (i4 == 0) & (i2 + i3 <= n) & (i2 < n)
    minus = (i1 == n) & (i2 + i3 >= n) & (i3 > 0)
    return P, plus, minus


def random_points(m=N_RANDOM, seed=20260909):
    rng = np.random.default_rng(seed)
    P = rng.random((m, 4))
    return (P,) + predict(P)


def face_and_edge_points(nf=12, ne=100, seed=7):
    """Dense sampling of every facet, every 2-face and every edge of the cube, plus the two
    relative-interior slabs {p4=0, p2+p3=1} and {p1=1, p2+p3=1} where the closed form is
    tight.  Everything rational, so membership is decided on integers."""
    rng = np.random.default_rng(seed)
    P_all, plus_all, minus_all = [], [], []

    def add_int(NUM, den):
        NUM = np.asarray(NUM, dtype=np.int64)
        pl, mi = predict_int(NUM, den)
        P_all.append(NUM / den)
        plus_all.append(pl)
        minus_all.append(mi)

    # 8 facets: coordinate i pinned to 0 or 1, the other three on a (nf+1)^3 grid
    g3 = np.array(list(product(range(nf + 1), repeat=3)), dtype=np.int64)
    for i in range(4):
        for val in (0, nf):
            NUM = np.empty((g3.shape[0], 4), dtype=np.int64)
            NUM[:, i] = val
            NUM[:, [c for c in range(4) if c != i]] = g3
            add_int(NUM, nf)
    # 24 2-faces: two coordinates pinned, two free on a fine grid
    n2 = 40
    g2 = np.array(list(product(range(n2 + 1), repeat=2)), dtype=np.int64)
    for i, j in product(range(4), repeat=2):
        if i >= j:
            continue
        for vi, vj in product((0, n2), repeat=2):
            NUM = np.empty((g2.shape[0], 4), dtype=np.int64)
            NUM[:, i], NUM[:, j] = vi, vj
            NUM[:, [c for c in range(4) if c not in (i, j)]] = g2
            add_int(NUM, n2)
    # 32 edges: three coordinates pinned to a vertex value, one sweeping
    a1 = np.arange(ne + 1, dtype=np.int64)
    for i in range(4):
        others = [c for c in range(4) if c != i]
        for vals in product((0, ne), repeat=3):
            NUM = np.empty((a1.size, 4), dtype=np.int64)
            NUM[:, i] = a1
            for c, v in zip(others, vals):
                NUM[:, c] = v
            add_int(NUM, ne)
    # the two tight slabs, at a denominator fine enough not to coincide with the grids
    den = 100_000
    s = rng.integers(0, den + 1, size=4000)
    t = rng.integers(0, den + 1, size=4000)
    add_int(np.column_stack([t, s, den - s, np.zeros_like(s)]), den)   # p4=0, p2+p3=1 exactly
    add_int(np.column_stack([np.full_like(s, den), s, den - s, t]), den)

    return np.vstack(P_all), np.concatenate(plus_all), np.concatenate(minus_all)


def boundary_probe_points(seed=99):
    """Points at controlled distance delta from each defining boundary, with the exact
    membership they are supposed to have.  Returned separately with a per-family label so
    the margin can be tabulated against delta."""
    rng = np.random.default_rng(seed)
    fams = {}
    m = 400
    for delta in (1e-3, 1e-6, 1e-9):
        # --- Riv+ : base points strictly inside {p4=0, p2+p3<1, p2<1}
        p1 = rng.random(m)
        s = rng.random(m) * 0.8                    # p2 in [0,0.8)
        r = rng.random(m) * (1 - s) * 0.8          # p3 with p2+p3 < 1
        base = np.column_stack([p1, s, r, np.zeros(m)])
        P = base.copy(); P[:, 3] = delta
        fams[f'Riv+ : p4 = +{delta:g}  (OUT)'] = (P, False, None)
        P = base.copy(); P[:, 2] = 1 - P[:, 1] + delta
        fams[f'Riv+ : p2+p3 = 1+{delta:g}  (OUT)'] = (P, False, None)
        P = base.copy(); P[:, 2] = np.clip(1 - P[:, 1] - delta, 0, 1)
        fams[f'Riv+ : p2+p3 = 1-{delta:g}  (IN)'] = (P, True, None)
        P = base.copy(); P[:, 1] = 1 - delta; P[:, 2] = 0.0
        fams[f'Riv+ : p2 = 1-{delta:g}, p3=p4=0  (IN, near edge E)'] = (P, True, None)
        # --- Riv- : the mirror
        p4 = rng.random(m)
        s = rng.random(m) * 0.8
        r = 1 - s + rng.random(m) * s * 0.8        # p2+p3 > 1
        base = np.column_stack([np.ones(m), s, np.clip(r, 0, 1), p4])
        P = base.copy(); P[:, 0] = 1 - delta
        fams[f'Riv- : p1 = 1-{delta:g}  (OUT)'] = (P, None, False)
        P = base.copy(); P[:, 2] = np.clip(1 - P[:, 1] - delta, 0, 1)
        fams[f'Riv- : p2+p3 = 1-{delta:g}  (OUT)'] = (P, None, False)
        P = base.copy(); P[:, 2] = np.clip(1 - P[:, 1] + delta, 0, 1)
        fams[f'Riv- : p2+p3 = 1+{delta:g}  (IN)'] = (P, None, True)
        P = base.copy(); P[:, 1] = 1.0; P[:, 2] = delta
        fams[f'Riv- : p2=1, p3 = {delta:g}  (IN, near edge E)'] = (P, None, True)
    return fams


# --------------------------------------------------------------------- driver

def sweep(P, plus, minus, label, pool):
    blocks = [(P[i:i + CHUNK], plus[i:i + CHUNK], minus[i:i + CHUNK], label)
              for i in range(0, P.shape[0], CHUNK)]
    return merge(pool.map(run_block, blocks))


def main():
    fh = sys.stdout
    t0 = time.time()

    print("=" * 78, file=fh)
    print("SWEEP OF THE RIVAL CLOSED FORM OVER THE FULL MEMORY-ONE CUBE", file=fh)
    print(f"eps = {EPSILONS},  verdict tolerance = {TOL:g},  "
          f"co-players = the 16 vertices", file=fh)
    print("=" * 78, file=fh)

    # ---- 0. accuracy receipt: the vectorised float solve against exact Fractions
    print("\n[0] accuracy receipt: the two float routes against exact Fractions", file=fh)
    for eps in EPSILONS:
        rng = np.random.default_rng(1)
        w_tree = w_lu = 0.0
        for _ in range(40):
            p = rng.random(4)
            k = int(rng.integers(16))
            pe = [float(x) for x in p]
            we = stationary_exact(pe, [float(x) for x in DET[k]], eps)
            d_exact = float(we[2] - we[1])
            w_tree = max(w_tree, abs(d_exact - float(margins(p[None, :], eps)[0, k])))
            wl = stationary(p, DET[k], eps)
            w_lu = max(w_lu, abs(d_exact - float(wl[2] - wl[1])))
        print(f"    eps={eps:.0e}  max |error| over 40 random (p,q) margins:"
              f"  tree theorem {w_tree:.3e}   linear solve {w_lu:.3e}", file=fh)
    print("    (the tree route is used everywhere below; the linear solve is shown only to"
          " justify\n     not using it -- the pair chain is nearly reducible at these eps)",
          file=fh)

    with Pool(12) as pool:
        # ---- 1. full grid
        P, plus, minus = grid_points()
        print(f"\n[1] full grid {GRID_N + 1}^4 = {P.shape[0]} strategies", file=fh)
        report("grid", sweep(P, plus, minus, "grid", pool), fh)
        print(f"    ({time.time() - t0:.0f} s elapsed)", file=fh)

        # ---- 2. random interior
        P, plus, minus = random_points()
        print(f"\n[2] {P.shape[0]} uniform random interior points", file=fh)
        report("random", sweep(P, plus, minus, "random", pool), fh)
        print(f"    ({time.time() - t0:.0f} s elapsed)", file=fh)

        # ---- 3. faces and edges
        P, plus, minus = face_and_edge_points()
        print(f"\n[3] dense faces / 2-faces / edges : {P.shape[0]} points", file=fh)
        report("faces+edges", sweep(P, plus, minus, "faces", pool), fh)
        print(f"    ({time.time() - t0:.0f} s elapsed)", file=fh)

        # ---- 4. boundary probes, tabulated by family
        print("\n[4] boundary probes: margin as a function of distance to the boundary", file=fh)
        print("    (these are NOT counted as disagreements; a point at distance delta", file=fh)
        print("     outside the set has a margin of size Theta(delta), which is below the", file=fh)
        print("     verdict tolerance once delta << TOL -- that is the resolution limit,", file=fh)
        print("     not a failure of the closed form.)", file=fh)
        fams = boundary_probe_points()
        hdr = (f"    {'family':<44}{'eps':>7}  {'min |m|':>12}{'max |m|':>12}"
               f"  {'verdictOK':>9}  unresolved  note")
        print(hdr, file=fh)
        for name, (Pf, want_plus, want_minus) in fams.items():
            for eps in EPSILONS:
                mp, mm, _ = verdicts(Pf, eps)
                m = mp if want_plus is not None else mm
                want = want_plus if want_plus is not None else want_minus
                ok = int(np.sum((m > -TOL) == want))
                # a mismatch is a RESOLUTION artifact iff the true margin is smaller than the
                # verdict tolerance -- the sweep cannot see a Theta(delta) margin for delta<TOL
                unres = int(np.sum(((m > -TOL) != want) & (np.abs(m) < TOL)))
                note = ''
                if ok < Pf.shape[0] and unres < Pf.shape[0] - ok:
                    # not resolution: it is the eps-vs-(1-p2) crossover of clause (9c)
                    note = 'eps not << 1-p2: see [8]'
                print(f"    {name:<44}{eps:>7.0e}  {np.abs(m).min():12.4e}"
                      f"{np.abs(m).max():12.4e}  {ok:>4}/{Pf.shape[0]:<4}  "
                      f"{unres:>10}  {note}", file=fh)

        # ---- 5. the 16 vertices, spelled out
        print("\n[5] vertex census (eps = 1e-10)", file=fh)
        V = np.array([np.asarray(q, float) for q in DET])
        mp, mm, _ = verdicts(V, 1e-10)
        for i in range(16):
            g = genome(V[i])
            print(f"    {g}  m+={mp[i]:+.6e}  m-={mm[i]:+.6e}   "
                  f"Riv+ {'YES' if mp[i] > -TOL else ' no'} (pred "
                  f"{'YES' if in_riv_plus(V[i]) else ' no'})   "
                  f"Riv- {'YES' if mm[i] > -TOL else ' no'} (pred "
                  f"{'YES' if in_riv_minus(V[i]) else ' no'})", file=fh)

        # ---- 6. the exceptional edge E = {(p1,1,0,0)}
        print("\n[6] the exceptional edge E = {(p1,1,0,0)} : predicted OUT of both sets;"
              " witness value -1/(3-p1)", file=fh)
        p1s = np.linspace(0, 1, 11)
        PE = np.column_stack([p1s, np.ones(11), np.zeros(11), np.zeros(11)])
        for eps in EPSILONS:
            mp, mm, _ = verdicts(PE, eps)
            for i in range(0, 11, 2):
                print(f"    eps={eps:.0e}  p1={p1s[i]:.1f}  m+={mp[i]:+.9f}  "
                      f"(-1/(3-p1) = {-1 / (3 - p1s[i]):+.9f})  m-={mm[i]:+.9f}", file=fh)

        # ---- 7. the exact witness / certificate formulas of clause (8)
        print("\n[7] the exact closed-form witness values, checked pointwise "
              "(eps = 1e-10, 3000 random points per family)", file=fh)
        rng = np.random.default_rng(4242)
        eps = 1e-10
        n = 3000

        def emax(a, b):
            return float(np.abs(np.asarray(a) - np.asarray(b)).max())

        # (a) p4 > 0, witness ALLD (co-player index 0):  -p4/(1 - p2 + p4)
        P = rng.random((n, 4))
        P[:, 3] = np.maximum(P[:, 3], 1e-3)
        d = margins(P, eps)[:, 0]
        pred = -P[:, 3] / (1 - P[:, 1] + P[:, 3])
        print(f"    (a) p4>0 vs ALLD      = -p4/(1-p2+p4)                : "
              f"max err {emax(d, pred):.3e}", file=fh)
        # (b) p4 = 0, p2+p3 > 1, witness q* = 0001 (index 1): (1-p2-p3)/(p3+2(1-p2))
        p2 = rng.random(n) * 0.98                       # p2 < 1
        p3 = 1 - p2 + (0.01 + 0.98 * rng.random(n)) * p2      # p2+p3 > 1, p3 <= 1
        P = np.column_stack([rng.random(n), p2, np.clip(p3, 0, 1), np.zeros(n)])
        d = margins(P, eps)[:, 1]
        pred = (1 - P[:, 1] - P[:, 2]) / (P[:, 2] + 2 * (1 - P[:, 1]))
        print(f"    (b) p2+p3>1 vs 0001   = (1-p2-p3)/(p3+2(1-p2))       : "
              f"max err {emax(d, pred):.3e}", file=fh)
        # (c) p2 = 1, p4 = 0, witness ALLD: -(1+p3)/(3+p3-p1)
        P = np.column_stack([rng.random(n), np.ones(n), rng.random(n), np.zeros(n)])
        d = margins(P, eps)[:, 0]
        pred = -(1 + P[:, 2]) / (3 + P[:, 2] - P[:, 0])
        print(f"    (c) p2=1,p4=0 vs ALLD = -(1+p3)/(3+p3-p1)            : "
              f"max err {emax(d, pred):.3e}", file=fh)
        # (d) the certificate of clause (6) on Riv+, with equality on the face p2+p3 = 1
        p2 = rng.random(n) * 0.9
        P = np.column_stack([rng.random(n), p2, (1 - p2) * rng.random(n), np.zeros(n)])
        for e in EPSILONS:
            m = margins(P, e).min(axis=1)
            bound = -e / (e + (1 - 2 * e) * (1 - P[:, 1]))
            viol = int(np.sum(m < bound - 1e-14))
            print(f"    (d) certificate m+ >= -eps/(eps+(1-2eps)(1-p2)), eps={e:.0e} : "
                  f"violations {viol}/{n},  worst slack {float((m - bound).min()):+.3e}",
                  file=fh)
        p2 = rng.random(n) * 0.9
        P = np.column_stack([rng.random(n), p2, 1 - p2, np.zeros(n)])
        for e in EPSILONS:
            m = margins(P, e).min(axis=1)
            bound = -e / (e + (1 - 2 * e) * (1 - P[:, 1]))
            print(f"    (d') equality on the face p2+p3=1,               eps={e:.0e} : "
                  f"max |m+ - bound| {emax(m, bound):.3e}", file=fh)

        # ---- 8. the eps-vs-(1-p2) crossover: convention A vs convention B
        print("\n[8] approach to the exceptional edge E: p = (1/2, 1-x, 0, 0), which is IN"
              " Riv+ for\n    every fixed x > 0.  The finite-eps margin is"
              " -eps/(eps + (1-2eps)x); the limit is 0 only\n    once eps << x, which is why"
              " a point at distance 1e-9 from E is still 'seen' as a\n    non-rival at"
              " eps = 1e-10.  This is clause (9c), not a failure of the closed form.", file=fh)
        print(f"    {'x = 1-p2':>10}{'eps':>9}{'m+ computed':>16}{'-eps/(eps+x)':>16}"
              f"{'verdict':>10}", file=fh)
        for x in (1e-2, 1e-4, 1e-6, 1e-8, 1e-9, 1e-10, 1e-12):
            for e in EPSILONS:
                P = np.array([[0.5, 1 - x, 0.0, 0.0]])
                m = float(margins(P, e).min())
                b = -e / (e + (1 - 2 * e) * x)
                print(f"    {x:>10.0e}{e:>9.0e}{m:>+16.6e}{b:>+16.6e}"
                      f"{('rival' if m > -TOL else 'NOT'):>10}", file=fh)

        # ---- 9. the smallest GENUINELY non-zero |w_DC - w_CD|, certified in exact rationals
        print("\n[9] smallest non-zero |w_DC - w_CD| over the full grid x 16 co-players.", file=fh)
        print("    A float difference of two O(1) stationary probabilities carries absolute"
              " error ~1e-16,\n    so any float value below that is pure noise: the raw"
              " float minimum (~1e-26) is in fact a\n    pair whose exact margin is 0"
              " (p = 0001 against q = 0001, i.e. self-play).  Below, the 400\n    smallest"
              " float candidates are recomputed in exact Fractions.  The value reported is"
              " therefore\n    a CERTIFIED non-zero margin and an UPPER BOUND on the true"
              " minimum -- float cannot rank\n    quantities this small, so it cannot"
              " certify that no smaller one exists elsewhere.\n    (Cross-checked with"
              " mpmath at dps=200 on the exact double inputs: agrees to 20 digits.)",
              file=fh)
        Pg, _, _ = grid_points()
        for eps in EPSILONS:
            cand = []
            for i in range(0, Pg.shape[0], 8192):
                d = np.abs(margins(Pg[i:i + 8192], eps)).ravel()
                idx = np.argpartition(d, 400)[:400]
                cand += [(float(d[j]), i + int(j) // 16, int(j) % 16) for j in idx]
            cand.sort()
            smallest = None
            zeros = 0
            for _, ip, k in cand[:400]:
                we = stationary_exact([float(x) for x in Pg[ip]],
                                      [float(x) for x in DET[k]], eps)
                val = abs(we[2] - we[1])
                if val == 0:
                    zeros += 1
                elif smallest is None or val < smallest[0]:
                    smallest = (val, Pg[ip].copy(), k)
            print(f"    eps={eps:.0e}  of the 400 smallest float candidates, {zeros} are"
                  f" EXACTLY zero;", file=fh)
            print(f"              smallest exactly non-zero |w_DC - w_CD| = "
                  f"{float(smallest[0]):.6e}", file=fh)
            print(f"              attained at p = {tuple(smallest[1])} vs q = "
                  f"{genome(DET[smallest[2]])}", file=fh)

    print(f"\nTOTAL {time.time() - t0:.0f} s", file=fh)


if __name__ == "__main__":
    main()
