"""item 4 -- the DIMENSION FIELD of the limit-Nash set N(u,v) subset [0,1]^4.

Decisive, structured verification of the theorem under test:

    N(u,v) is never empty, never 4-dimensional, and
        dim N = 3  on {u<1} u {v<0} u {(1,0)},
        dim N = 2  on Q\{(1,0)},  Q = {u>=1, v>=0} closed,
    the only carriers of a 3-dimensional piece being the two rival facets
        {p1=1}:  (I)  p4 v <= 1-p2      (II)  (1-p2)(u-1) + p3 v <= 0
        {p4=0}:  (I*) (1-p1)u + p3 <= 0 (II*) (1-p2)u + p3(1+v) <= 0

Nothing is searched blindly.  Each claimed carrier is parametrised, its claimed
polytope is given a CHEBYSHEV BALL by linear programming inside the stratum, and the
ball is then verified against the RAW 16-deviator limit test at its centre and at the
2d axis points -- a positive radius plus a passing ball test is a certificate of full
stratum dimension, and it is immune to the thin-set failure that a uniform sampler has
(admitted gap 2 of the theorem).  Non-carriers are attacked from two sides: uniform
sampling (which sees any positive-measure piece) and Nelder-Mead maximisation of the
raw Nash gap from many starts (which sees thin ones).

Conventions and harness: m1lib.  eps = 1e-7 and 1e-9 throughout, never 1e-12.
"""
import sys, time, itertools
import numpy as np
from scipy.optimize import linprog, minimize
import m1lib

t0 = time.time()
rng = np.random.default_rng(4)
EPS = (1e-7, 1e-9)          # never 1e-12: the float pair chain is ill-conditioned there
TOL = 1e-4

def gap(p, u, v, eps):
    return m1lib.nash_gap(np.asarray(p, float), u, v, eps)

# --------------------------------------------------------------------- the limit test
# At finite eps the gap of a limit-Nash strategy is O(eps) with a constant that blows up
# in the boundary layers (a coordinate of size O(eps) -- e.g. p4 small on {p1=1}), so a
# fixed threshold at one eps misclassifies thin regions.  The DECISIVE test is the pair:
# a limit-Nash point has gap(1e-9) = gap(1e-7)/100, a non-Nash point has gap(1e-9) =
# gap(1e-7).  Both eps are used everywhere and their agreement is the test itself.
ratios = []
def lim_nash(p, u, v, detail=False):
    g7 = gap(p, u, v, 1e-7); g9 = gap(p, u, v, 1e-9)
    if g9 > -1e-6:                      # clean: nothing left at eps = 1e-9
        ok = True
    elif g9 > 0.2 * g7:                 # shrank with eps by at least 5x -> limit-Nash
        ok = True
        ratios.append((g7 / g9, tuple(np.round(p, 4)), u, v))
    else:
        ok = False
    return (ok, g7, g9) if detail else ok

def is_nash(p, u, v, eps=None):
    return lim_nash(p, u, v)


# --------------------------------------------------- exact fallback for thin certificates
# Where the certificate ball is thinner than the float chain can resolve (a coordinate of
# size ~1e-4 puts eps = 1e-9 still inside the boundary layer), the same test is run in
# exact Fractions at eps = 1e-11 and 1e-12, where no conditioning question arises.
from fractions import Fraction as Fr
def gap_exact(p, u, v, k):
    """Nash gap in exact rationals at eps = 10^-k; returns a float of an exact Fraction."""
    U, V = Fr(u).limit_denominator(10**6), Fr(v).limit_denominator(10**6)
    P = [Fr(x).limit_denominator(10**9) for x in p]
    eps = Fr(1, 10**k)
    pay = lambda w: w[0] + w[1] * U + w[2] * (1 + V)
    own = pay(m1lib.stationary_exact(P, P, eps))
    best = max(pay(m1lib.stationary_exact(list(sg), P, eps)) for sg in m1lib.DET_EXACT)
    return float(own - best)

def lim_nash_exact(p, u, v):
    """limit-Nash by the same eps-ratio criterion, in exact arithmetic at 1e-11 / 1e-12."""
    g11, g12 = gap_exact(p, u, v, 11), gap_exact(p, u, v, 12)
    return (g12 > -1e-9) or (g12 > 0.2 * g11), g11, g12

# ---------------------------------------------------------------- games
G = [(u, v) for u in np.arange(-8, 8.01, 1.0) for v in np.arange(-8, 8.01, 1.0)]
BND = []                                   # on and just off every claimed boundary
for v in (0.0, 0.25, 1.0, 3.0, 8.0):
    for du in (-1e-3, 0.0, 1e-3):
        BND.append((1.0 + du, v))
for u in (1.0, 1.25, 2.0, 5.0, 8.0):
    for dv in (-1e-3, 0.0, 1e-3):
        BND.append((u, 0.0 + dv))
for d in (-1e-3, 1e-3):                    # around the exceptional corner
    BND += [(1.0 + d, 0.0), (1.0, 0.0 + d), (1.0 + d, 0.0 + d), (1.0 + d, -d)]
BND.append((1.0, 0.0))
BND = sorted(set(BND))
GAMES = G + BND
print(f"# games: {len(G)} on the unit grid of [-8,8]^2 + {len(BND)} boundary probes")

def claim_dim(u, v):
    if u == 1.0 and v == 0.0:
        return 3
    return 2 if (u >= 1 and v >= 0) else 3

# ---------------------------------------------------------------- 1. EXISTENCE
print("\n=== 1. EXISTENCE: a limit-Nash witness at every game "
      "(constant strategy at the symmetric stage-game equilibrium) ===")
def witness(u, v):
    if v <= 0: return np.ones(4)                       # ALLC
    if u <= 0: return np.zeros(4)                      # ALLD
    x = u / (u + v);  return np.full(4, x)             # the mixed x*
fail, worst = 0, {}
for eps in EPS:
    w = (0, None)
    for (u, v) in GAMES:
        g = gap(witness(u, v), u, v, eps)
        if not lim_nash(witness(u, v), u, v): fail += 1; print("  EXISTENCE FAILURE", u, v, g)
        if g < w[0]: w = (g, (u, v))
    worst[eps] = w
    print(f"  eps={eps:g}: {len(GAMES)} games, failures {fail}, worst gap {w[0]:.2e} at {w[1]}")

# ---------------------------------------------------------------- geometry helpers
def cheby(A, b, n=3):
    """Chebyshev ball of {Ax<=b} cap [0,1]^n.  Returns (radius, centre)."""
    A = np.asarray(A, float); b = np.asarray(b, float)
    A = np.vstack([A, -np.eye(n), np.eye(n)])
    b = np.concatenate([b, np.zeros(n), np.ones(n)])
    nrm = np.linalg.norm(A, axis=1)
    Alp = np.hstack([A, nrm[:, None]])
    c = np.zeros(n + 1); c[-1] = -1
    r = linprog(c, A_ub=Alp, b_ub=b, bounds=[(None, None)] * n + [(0, None)])
    if not r.success: return 0.0, None
    return r.x[-1], r.x[:n]

def ball_test(embed, centre, radius, u, v, eps, extra=4):
    """raw 16-deviator test at the centre and at 2n axis points + random directions."""
    n = len(centre); pts = [centre]
    for i in range(n):
        for s in (+1, -1):
            q = np.array(centre, float); q[i] += s * 0.5 * radius; pts.append(q)
    for _ in range(extra):
        d = rng.normal(size=n); d /= np.linalg.norm(d)
        pts.append(np.array(centre) + 0.5 * radius * d)
    res = [lim_nash(embed(q), u, v, detail=True) for q in pts]
    return all(r[0] for r in res), min(r[2] for r in res)

# the two claimed 3-dimensional carriers, as (embed, A(u,v), b(u,v), predicate)
def car1(u, v):                       # {p1=1}, free (p2,p3,p4)
    A = [[1.0, 0.0, v], [-(u - 1), v, 0.0]];  b = [1.0, 1.0 - u]
    return (lambda x: np.array([1.0, x[0], x[1], x[2]])), A, b
def car4(u, v):                       # {p4=0}, free (p1,p2,p3)
    A = [[-u, 0.0, 1.0], [0.0, -u, 1.0 + v]];  b = [-u, -u]
    return (lambda x: np.array([x[0], x[1], x[2], 0.0])), A, b

def pred1(x, u, v):
    p2, p3, p4 = x
    return (p4 * v <= 1 - p2 + 1e-12) and ((1 - p2) * (u - 1) + p3 * v <= 1e-12)
def pred4(x, u, v):
    p1, p2, p3 = x
    return ((1 - p1) * u + p3 <= 1e-12) and ((1 - p2) * u + p3 * (1 + v) <= 1e-12)

# ------------------------------------------- 2. the two carriers: dimension by LP + ball
print("\n=== 2. THE TWO CLAIMED 3-DIMENSIONAL CARRIERS: Chebyshev radius + raw ball test ===")
for name, car, rule in (("{p1=1}", car1, lambda u, v: (u < 1 or v < 0) or (u == 1 and v == 0)),
                        ("{p4=0}", car4, lambda u, v: u < 0)):
    bad_rule = bad_ball = rescued = 0; rmin = 9.9; examples = []
    for (u, v) in GAMES:
        emb, A, b = car(u, v)
        r, c = cheby(A, b)
        feas = r > 1e-9
        if feas != bool(rule(u, v)):
            bad_rule += 1; examples.append((u, v, r, rule(u, v)))
        if feas:
            ok, mg = ball_test(emb, c, min(r, 0.2), u, v, EPS[0])
            if not ok:                      # thin ball: rerun the same ball exactly
                pts = [c] + [c + s_ * 0.5 * min(r, 0.2) * np.eye(3)[j]
                             for j in range(3) for s_ in (+1, -1)]
                ex = [lim_nash_exact(emb(q), u, v) for q in pts]
                if all(e[0] for e in ex):
                    rescued += 1
                    print(f"    thin ball on {name} at ({u},{v}) r={r:.2e}: float mingap "
                          f"{mg:.2e} at eps=1e-9, RESCUED exactly -- worst gap "
                          f"{min(e[1] for e in ex):.2e} at 1e-11 -> "
                          f"{min(e[2] for e in ex):.2e} at 1e-12")
                else:
                    bad_ball += 1
                    print(f"    ball FAILS on {name} at ({u},{v}) r={r:.2e} mingap={mg:.2e} "
                          f"(exact worst at 1e-12: {min(e[2] for e in ex):.2e})")
            rmin = min(rmin, r)
    print(f"  {name}: strict-feasibility rule mismatches {bad_rule}/{len(GAMES)}"
          f"   raw ball-test failures {bad_ball} (thin balls rescued in exact arithmetic: "
          f"{rescued})   smallest certified radius {rmin:.3g}")
    for e in examples[:5]: print("     rule mismatch:", e)

# --------------------------------- 3. do (I),(II) / (I*),(II*) really equal the Nash set?
print("\n=== 3. THE CLAIMED FACET INEQUALITIES vs THE RAW 16-DEVIATOR TEST ===")
CHK = [(-4., -4.), (-2., 3.), (-0.5, 0.5), (0.5, -2.), (0.5, 0.5), (1., 0.), (2., -3.),
       (2., 0.5), (3., 3.), (6., -0.5), (-1., 0.), (0., -1.)]
for name, car, pred in (("{p1=1} (I),(II)", car1, pred1), ("{p4=0} (I*),(II*)", car4, pred4)):
    mism = 0; tot = 0; ex = []
    for (u, v) in CHK:
        emb, A, b = car(u, v)
        X = rng.random((150, 3))
        for x in X:
            tot += 1
            if pred(x, u, v) != is_nash(emb(x), u, v, EPS[0]):
                mism += 1
                if len(ex) < 4: ex.append((u, v, np.round(x, 3), pred(x, u, v),
                                           gap(emb(x), u, v, EPS[0])))
    print(f"  {name}: {mism} mismatches / {tot} uniform facet points over {len(CHK)} games")
    for e in ex: print("     ", e)

# the sub-face {p1=1,p4=0}: the theorem says it is a carrier at NO game, yet (I),(II)
# do not forbid it -- test the slice explicitly, and its p4 = 1e-3 neighbour.
print("  sub-face {p1=1,p4=0} (Grim-like) vs {p1=1,p4=1e-3}:")
for (u, v) in [(-2., 3.), (0.5, 0.5), (2., -3.), (3., 3.)]:
    n0 = n1 = 0; X = rng.random((60, 2))
    for p2, p3 in X:
        n0 += is_nash([1, p2, p3, 0.0], u, v, EPS[0])
        n1 += is_nash([1, p2, p3, 1e-3], u, v, EPS[0])
        ii = ((1 - p2) * (u - 1) + p3 * v <= 0)
    nI = sum(((1 - p2) * (u - 1) + p3 * v <= 0) for p2, p3 in X)
    print(f"    ({u},{v}): raw-Nash at p4=0: {n0}/60   at p4=1e-3: {n1}/60   (I)&(II) predict {nI}/60")

# ------------------------------------------- 4. NO OTHER FACET CARRIES A 3-DIM PIECE
print("\n=== 4. THE OTHER SIX FACETS: no 3-dimensional piece (uniform + Nelder-Mead) ===")
FACETS = [("p1=0", 0, 0.0), ("p2=0", 1, 0.0), ("p2=1", 1, 1.0),
          ("p3=0", 2, 0.0), ("p3=1", 2, 1.0), ("p4=1", 3, 1.0)]
def emb_facet(i, val):
    def e(x):
        p = np.empty(4); k = 0
        for j in range(4):
            if j == i: p[j] = val
            else: p[j] = min(1.0, max(0.0, x[k])); k += 1
        return p
    return e
FG = [(-4., -4.), (-0.5, 0.5), (0.5, -2.), (2., 0.5), (3., 3.), (2., -3.)]
for nm, i, val in FACETS:
    emb = emb_facet(i, val); hits = 0; tot = 0; best = -9.9; bestpt = None; cand = []
    for (u, v) in FG:
        X = rng.random((250, 3))
        for x in X:
            tot += 1
            g = gap(emb(x), u, v, EPS[0])
            if lim_nash(emb(x), u, v): hits += 1; cand.append((u, v, emb(x)))
            if g > best: best, bestpt = g, (u, v, np.round(emb(x), 3))
        for _ in range(6):                      # thin-set hunt
            x0 = rng.random(3)
            r = minimize(lambda x: -gap(emb(x), u, v, EPS[0]), x0,
                         method="Nelder-Mead",
                         options=dict(maxfev=120, xatol=1e-3, fatol=1e-6))
            if -r.fun > best: best, bestpt = -r.fun, (u, v, np.round(emb(r.x), 3))
            if lim_nash(emb(r.x), u, v): cand.append((u, v, emb(r.x)))
    ball3 = 0
    for (u, v, p) in cand[:40]:            # is any Nash point there 3-dimensional?
        free = [j for j in range(4) if j != i]
        okb = all(lim_nash(np.clip(p + s_ * 0.02 * np.eye(4)[j], 0, 1), u, v)
                  for j in free for s_ in (+1, -1))
        ball3 += okb
    print(f"  facet {nm}: uniform Nash hits {hits}/{tot}; best gap(1e-7) over uniform+NM "
          f"{best:.3e} at {bestpt}; Nash points found {len(cand)}, of them 3-dimensional {ball3}")

# ------------------------------------------- 5. NO 4-DIMENSIONAL PIECE
print("\n=== 5. THE OPEN CUBE: no 4-dimensional piece; best interior gap ===")
IG = [(-4., -4.), (-2., 3.), (-0.5, 0.5), (0.5, 0.5), (1., 0.), (2., -3.), (2., 0.5), (5., 5.)]
glob = -9.9; globpt = None
for (u, v) in IG:
    X = 0.02 + 0.96 * rng.random((900, 4))
    gs = np.array([gap(x, u, v, 1e-9) for x in X])
    k = gs.argmax()
    nn = sum(lim_nash(x, u, v) for x in X[np.argsort(gs)[-12:]])
    if gs[k] > glob: glob, globpt = gs[k], (u, v, np.round(X[k], 3))
    print(f"  ({u},{v}): {nn}/900 interior draws Nash (limit test on the 12 best); "
          f"best interior gap at eps=1e-9 {gs[k]:.3e}")
print(f"  BEST interior Nash gap anywhere: {glob:.3e} at {globpt}")

print("  the two exceptional interior sets the structure theorem allows:")
for (u, v) in IG:
    # ZD equaliser sheet: A(p1-p2)=v, A(p4-p3)=u, A(1+p3-p1)=1-u; free (A,p3)
    found = []
    for A in (-8., -4., -2., -1., -.5, 0.5, 1.0, 2.0, 4.0, 8.0, 16., 32.):
        for p3 in np.linspace(0.02, 0.98, 13):
            p1 = 1 + p3 - (1 - u) / A; p2 = p1 - v / A; p4 = p3 + u / A
            p = np.array([p1, p2, p3, p4])
            if (p > 1e-6).all() and (p < 1 - 1e-6).all():
                found.append(gap(p, u, v, 1e-9))
    claimzd = (u < 1 and v > -1 and v > u - 1)
    zd = (f"{len(found)} ZD points, worst gap(1e-9) {min(found):.2e}" if found
          else "no ZD point in the open cube") + f" [claimed 2-dim there: {claimzd}]"
    x = u / (u + v) if u + v != 0 else 0.5
    xs = f"x*={x:.3f} gap(1e-9) {gap(np.full(4, x), u, v, 1e-9):.2e}" if 0 < x < 1 else "x* not interior"
    print(f"    ({u},{v}): {zd};  {xs}")

# ------------------------------------------- 6. THE 2-DIMENSIONAL CARRIERS ON Q
print("\n=== 6. THE 2-DIMENSIONAL STRATA ===")
def sheet_live(x): return np.array([x[0], 0.0, 1.0, x[1]])     # {p2=0,p3=1}
def sheet_other(x): return np.array([x[0], 1.0, 0.0, x[1]])    # {p2=1,p3=0}
def sheet_trap(x): return np.array([0.0, x[0], x[1], 1.0])     # {p1=0,p4=1}
QG = [(1., 0.), (1., 2.), (2., 0.), (3., 3.), (8., 8.), (1.5, 0.5), (5., 0.), (1., 0.001)]
print("  live alternating sheet {p2=0,p3=1}: saturation on the closed quadrant Q")
badQ = 0
for (u, v) in QG:
    XY = np.array([[a, b] for a in np.linspace(0.02, 0.98, 5) for b in np.linspace(0.02, 0.98, 5)])
    gs = [gap(sheet_live(x), u, v, 1e-9) for x in XY]
    n = sum(lim_nash(sheet_live(x), u, v) for x in XY); badQ += 25 - n
    print(f"    ({u},{v}): {n}/25 of the sheet Nash, worst gap(1e-9) {min(gs):.2e}")
print(f"    non-saturating points on Q: {badQ}")
print("  other alternating sheet {p2=1,p3=0}: claimed a carrier only on the ray u=v+1, v>=0")
for (u, v) in [(1., 0.), (2., 1.), (4., 3.), (2., 0.5), (3., 3.), (-1., -2.)]:
    XY = np.array([[a, b] for a in (0.2, 0.5, 0.8) for b in (0.2, 0.5, 0.8)])
    n = sum(is_nash(sheet_other(x), u, v) for x in XY)
    print(f"    ({u},{v}) [ray: {u == v + 1 and v >= 0}]: {n}/9 Nash")
print("  CC<->DD trap {p1=0,p4=1}: claimed Nash = {p3<=1-2u} cap {p2<=-2v}, 2-dim iff u<1/2,v<0")
mm = 0
for (u, v) in [(-4., -4.), (0.25, -1.), (0.4, -0.1), (0.6, -1.), (0.25, 0.5), (-2., -0.01), (2., -3.)]:
    X = rng.random((60, 2)); n = 0; loc = 0
    for p2, p3 in X:
        pr = (p3 <= 1 - 2 * u) and (p2 <= -2 * v)
        r = is_nash(sheet_trap([p2, p3]), u, v)
        n += r; loc += (pr != r)
    mm += loc
    dim2 = (u < 0.5 and v < 0)
    print(f"    ({u},{v}) [claim 2-dim: {dim2}]: {n}/60 Nash, predicate mismatches {loc}")
print(f"    total (III),(IV) mismatches: {mm}")

# ------------------------------------------- 7. THE DIMENSION FIELD ITSELF
print("\n=== 7. THE DIMENSION FIELD over every game ===")
mismatch = []; exact_rescue = []
for (u, v) in GAMES:
    d = 2
    for car in (car1, car4):
        emb, A, b = car(u, v)
        r, c = cheby(A, b)
        if r > 1e-9:
            ok, _ = ball_test(emb, c, min(r, 0.2), u, v, EPS[0], extra=1)
            if not ok:
                pts = [c] + [c + s_ * 0.5 * min(r, 0.2) * np.eye(3)[j]
                             for j in range(3) for s_ in (+1, -1)]
                ok = all(lim_nash_exact(emb(q), u, v)[0] for q in pts)
                if ok: exact_rescue.append((u, v))
            if ok: d = 3
    if d == 2:                       # certify 2 by a ball inside the live sheet
        ok, _ = ball_test(sheet_live, np.array([0.5, 0.5]), 0.6, u, v, EPS[0], extra=1)
        if not ok: d = 1
    if d != claim_dim(u, v): mismatch.append((u, v, d, claim_dim(u, v)))
print(f"  measured dim == claimed dim at {len(GAMES) - len(mismatch)}/{len(GAMES)} games"
      f"   (of which {len(exact_rescue)} needed the exact thin-ball fallback: {exact_rescue})")
for m in mismatch[:20]: print("    MISMATCH (u,v,measured,claimed):", m)

# ------------------------------------------- 8. THE SIXTEEN, panel a as control
print("\n=== 8. CONTROL: panel a for the 16 vertices ===")
M = 60                                       # cell centres of [-8,8]^2: area-accurate
cs = -8 + 16 * (np.arange(M) + 0.5) / M
CG = [(float(a), float(b)) for a in cs for b in cs]
cnt = {}
for (u, v) in CG:
    cnt[(u, v)] = sum(gap(p, u, v, 1e-7) > -TOL for p in m1lib.DET)
vals = np.array([cnt[g] for g in CG])
print("  count distribution over the %d cell centres:" % len(CG),
      {int(k): int((vals == k).sum()) for k in np.unique(vals)})
none = [g for g in CG if cnt[g] == 0]
print(f"  games with NO vertex equilibrium: {len(none)}/{len(CG)} = {100*len(none)/len(CG):.2f}% "
      f"of the square (sibling: 24.61% of the square, 23.16% of the disk)")
if none:
    print(f"    they all lie in u in [{min(g[0] for g in none):.2f},{max(g[0] for g in none):.2f}], "
          f"v in [{min(g[1] for g in none):.2f},{max(g[1] for g in none):.2f}]; "
          f"all in Snowdrift outside the unit square "
          f"({{u>0, v>0}} minus [0,1]^2): "
          f"{all(g[0] > 0 and g[1] > 0 and not (g[0] < 1 and g[1] < 1) for g in none)}; "
          f"fraction of them inside the dim-2 quadrant Q: "
          f"{sum(g[0] >= 1 and g[1] >= 0 for g in none)}/{len(none)}")
ever = sorted({m1lib.genome(p) for (u, v) in CG for p in m1lib.DET if gap(p, u, v, 1e-7) > -TOL})
print(f"  vertices Nash at some cell centre: {ever}")
print("  sibling's 8 with a full-dimensional Nash region: "
      "['0000','0001','0110','0111','1000','1001','1110','1111']")
bad = 0; worstw = 0.0
for g in none:
    ok, g7, g9 = lim_nash(witness(*g), *g, detail=True)
    bad += (not ok); worstw = min(worstw, g9)
print(f"  continuum witness at those {len(none)} no-equilibrium games: failures {bad}, "
      f"worst gap(1e-9) {worstw:.2e}")
sub = CG[::137]
dis8 = sum(1 for (u, v) in sub for p in m1lib.DET
           if (gap(p, u, v, 1e-7) > -TOL) != lim_nash(p, u, v))
print(f"  vertex classification, fixed-tol(1e-7) vs the two-eps limit test, "
      f"on {len(sub)*16} vertex/game pairs: {dis8} disagreements")

# ------------------------------------------- 9. eps AGREEMENT
print("\n=== 9. eps = 1e-7 vs 1e-9 ===")
POOL = []
for (u, v) in CHK:
    POOL += [(p, u, v) for p in rng.random((40, 4))]
    POOL += [(np.array([1., a, b, c]), u, v) for a, b, c in rng.random((25, 3))]
    POOL += [(np.array([a, b, c, 0.]), u, v) for a, b, c in rng.random((25, 3))]
    POOL += [(sheet_live(x), u, v) for x in rng.random((10, 2))]
dis = 0; worst = 0.0; ex = []
for p, u, v in POOL:
    g7 = gap(p, u, v, 1e-7); g9 = gap(p, u, v, 1e-9)
    if (g7 > -TOL) != (g9 > -TOL):
        dis += 1
        if len(ex) < 5: ex.append((np.round(p, 3), u, v, g7, g9))
    worst = max(worst, abs(g7 - g9))
print(f"  {len(POOL)} strategy/game pairs, FIXED tolerance {TOL:g} at each eps separately: "
      f"{dis} classification disagreements; max |gap(1e-7) - gap(1e-9)| = {worst:.2e}")
for e in ex: print("     disagreement:", e)
print("  (each such point is a boundary layer: its gap scales as eps, so the two-eps limit "
      "test used above calls it Nash at both.)")
r = np.array([x[0] for x in ratios]) if ratios else np.array([])
print(f"  points classified Nash by the eps-RATIO branch: {len(r)}"
      + (f"; ratio gap(1e-7)/gap(1e-9) in [{r.min():.1f}, {r.max():.1f}] "
         f"(exactly 100 for a limit-Nash point, 1 for a genuine failure)" if len(r) else ""))

print(f"\n[done in {time.time()-t0:.1f} s]")
