"""ROUTE B (item 2), EFFICIENCY ABOVE THE SWITCH LINE -- the ALTERNATORS.  All receipts.

THEOREM.  For p in [0,1]^4 under convention A, the self-play limit is W(p,p) = (0,1/2,1/2,0)
          -- equivalently pi(p,p) = Emax = (1+u+v)/2 for every game with u+v > 1 -- iff

    ALT = { p :  {p2,p3} = {0,1},   p1 < 1,   p4 > 0,   (p1,p4) != (0,1) }.

Six of the 81 strata; two unit squares in the (p1,p4)-plane (one over (p2,p3)=(0,1), one over
(1,0)), each minus the closed edges {p1=1} and {p4=0} and minus the further corner (p1,p4)=(0,1).

Run:  python3 alt_efficiency.py        (about 45 s)
"""
import itertools, random, sys, time
from fractions import Fraction

import numpy as np
import sympy as sp

sys.path.insert(0, '/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib
import alt_strata as A

SWAP = (0, 2, 1, 3)


def in_alt(p):
    p1, p2, p3, p4 = p
    return (((p2 == 0 and p3 == 1) or (p2 == 1 and p3 == 0))
            and p1 < 1 and p4 > 0 and not (p1 == 0 and p4 == 1))


def banner(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ---------------------------------------------------------------- 1. symbolic identities
banner("1.  EXACT SYMBOLIC IDENTITIES (all four p_i and eps symbolic)")
p1, p2, p3, p4, e = sp.symbols('p1 p2 p3 p4 epsilon')
hat = lambda x: e + (1 - 2 * e) * x


def statsym(p):
    P = [hat(x) for x in p]
    M = sp.zeros(4, 4)
    for s in range(4):
        a, b = P[s], P[SWAP[s]]
        M[s, 0], M[s, 1], M[s, 2], M[s, 3] = a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)
    w = sp.symbols('w0 w1 w2 w3')
    sol = sp.solve([sum(w[j] * M[j, i] for j in range(4)) - w[i] for i in range(3)]
                   + [sum(w) - 1], w, dict=True)[0]
    return [sp.cancel(sp.together(sol[x])) for x in w]


W = statsym([p1, p2, p3, p4])
print("  (a) self-play symmetry       w_CD - w_DC      =", sp.simplify(W[1] - W[2]))
print("  (b) p2 <-> p3 invariance     max_i |w_i - w_i(p2<->p3)| =",
      [sp.simplify(W[i] - W[i].subs({p2: p3, p3: p2}, simultaneous=True)) for i in range(4)])
a = 2 * hat(p1) * (1 - hat(p1)); b = 2 * hat(p4) * (1 - hat(p4))
al = hat(p2) * hat(p3); be = (1 - hat(p2)) * (1 - hat(p3))
print("  (c) cut balance  w_B(al+be) - (w_CC a + w_DD b) =",
      sp.simplify((W[1] + W[2]) * (al + be) - (W[0] * a + W[3] * b)))
L = sp.Matrix([[hat(p1)**2, a, (1 - hat(p1))**2],
               [al, 1 - al - be, be],
               [hat(p4)**2, b, (1 - hat(p4))**2]])
v = sp.symbols('v0 v1 v2')
s2 = sp.solve([sum(v[j] * L[j, i] for j in range(3)) - v[i] for i in range(2)] + [sum(v) - 1],
              v, dict=True)[0]
V = [sp.cancel(sp.together(s2[x])) for x in v]
print("  (d) lumping {CC},{CD,DC},{DD} exact at every eps:",
      [sp.simplify(x) for x in (W[0] - V[0], W[1] + W[2] - V[1], W[3] - V[2])])
h1, h4 = hat(p1), hat(p4); eta = e * (1 - e)
D = 2 * eta * (h4 + 1 - h1) + 2 * h4 * (1 - h1) * (1 + h1 - h4)
claim = [2 * eta * h4 / D, h4 * (1 - h1) * (1 + h1 - h4) / D,
         h4 * (1 - h1) * (1 + h1 - h4) / D, 2 * eta * (1 - h1) / D]
for c in [(0, 1), (1, 0)]:
    Wc = statsym([p1, sp.Integer(c[0]), sp.Integer(c[1]), p4])
    print("  (e) closed form at (p2,p3)=%s : exact minus claimed =" % (c,),
          [sp.simplify(Wc[i] - claim[i]) for i in range(4)])

# ---------------------------------------------------------------- 2. the 81 strata, two routes
banner("2.  ALL 81 STRATA, TWO INDEPENDENT ROUTES")
PS = [p1, p2, p3, p4]; X = sp.symbols('x1 x2 x3 x4', positive=True)
alt_tree, alt_lim, mism = [], [], 0
for st in itertools.product([0, 'i', 1], repeat=4):
    ok, lead, m = A.is_alternator(st)                       # leading-order tree theorem
    Wtree, _, _ = A.limit_W(st)
    sub = {PS[i]: (X[i] if st[i] == 'i' else sp.Integer(st[i])) for i in range(4)}
    Wlim = [sp.simplify(sp.limit(sp.cancel(x.subs(sub)), e, 0, '+')) for x in W]   # symbolic limit
    ok2 = (Wlim[0] == 0 and Wlim[3] == 0)
    if ok: alt_tree.append(st)
    if ok2: alt_lim.append(st)
    if ok != ok2 or any(sp.simplify(Wlim[i] - Wtree[i]) != 0 for i in range(4)): mism += 1
lab = lambda st: "".join(str(s) if s != 'i' else 'x' for s in st)
print("  disagreements between the two routes over all 81 strata (full W vector):", mism)
print("  alternator strata, tree theorem  :", [lab(s) for s in alt_tree])
print("  alternator strata, symbolic limit:", [lab(s) for s in alt_lim])
print("  closed form agrees on all 81 strata:",
      all((lab(st) in [lab(s) for s in alt_tree]) ==
          in_alt([0.5 if c == 'i' else c for c in st])
          for st in itertools.product([0, 'i', 1], repeat=4)))

# ---------------------------------------------------------------- 3. exact grids
banner("3.  EXACT-FRACTION GRIDS")
E4, E8 = Fraction(1, 10**4), Fraction(1, 10**8)
dev = lambda p, eps: (lambda w: w[0] + w[3])(m1lib.stationary_exact(list(p), list(p), eps))
G = [Fraction(i, 10) for i in range(11)]
t0 = time.time(); ma8 = 0.0; mn8 = 1.0; rat = []
for p in itertools.product(G, repeat=4):
    d8 = float(dev(p, E8))
    if in_alt(p): ma8 = max(ma8, d8); rat.append(float(dev(p, E4)) / max(d8, 1e-300))
    else: mn8 = min(mn8, d8)
print("  11^4 = 14641 points, eps=1e-8, %.0fs" % (time.time() - t0))
print("     predicted alternators : max (w_CC+w_DD) = %.3e" % ma8)
print("     predicted others      : min (w_CC+w_DD) = %.6f   (separation %.1e)" % (mn8, mn8 / ma8))
print("     dev(1e-4)/dev(1e-8) on alternators in [%.0f, %.0f]  (Theta(eps) => 1e4)"
      % (min(rat), max(rat)))
H = [Fraction(i, 20) for i in range(21)]; K = [Fraction(i, 10) for i in range(11)]
dis = 0
for q2 in H:
    for q3 in H:
        for q1 in K:
            for q4 in K:
                p = (q1, q2, q3, q4)
                if (float(dev(p, E8)) < 1e-4) != in_alt(p): dis += 1
print("  21^2 in (p2,p3) x 11^2 in (p1,p4) = 53361 points: disagreements =", dis)
rnd = random.Random(7); dis = 0
for _ in range(4000):
    p = tuple(min(Fraction(1), Fraction(rnd.randrange(0, 1000), 999)) if rnd.random() < 0.55
              else Fraction(rnd.choice([0, 1])) for _ in range(4))
    if (float(dev(p, E8)) < 1e-4) != in_alt(p): dis += 1
print("  4000 random mixed corner/interior rational points: disagreements =", dis)

# ---------------------------------------------------------------- 4. the sharp bound
banner("4.  THE SHARP NECESSITY BOUND   w_CC+w_DD >= 2s/(1+2s),  s = p2^p3 + (1-p2)^(1-p3) clamped")
worst = None
for p in itertools.product(G, repeat=4):
    h2, h3 = E8 + (1 - 2 * E8) * p[1], E8 + (1 - 2 * E8) * p[2]
    s = h2 * h3 + (1 - h2) * (1 - h3)
    sl = dev(p, E8) - 2 * s / (1 + 2 * s)
    if worst is None or sl < worst[0]: worst = (sl, p)
print("  11^4 exact: min slack = %.3e (>= 0 required), attained at p = %s"
      % (float(worst[0]), tuple(str(x) for x in worst[1])))

# ---------------------------------------------------------------- 5. the 16 vertices
banner("5.  THE 16 VERTICES -- the quantisation of ../../PartnersRivalsM1, exactly")
for pv in itertools.product([0, 1], repeat=4):
    Wv = statsym([sp.Integer(x) for x in pv])
    print("   %s  p2%sp3   w_CD(eps) = %-22s  limit %-4s %s"
          % ("".join(map(str, pv)), "=" if pv[1] == pv[2] else "!",
             sp.simplify(Wv[1]), sp.limit(Wv[1], e, 0, '+'),
             "ALT" if in_alt(pv) else ""))
print("  deterministic alternators:", [pv for pv in itertools.product([0, 1], repeat=4) if in_alt(pv)])

# ---------------------------------------------------------------- 6. friendly rivals above the line
banner("6.  FRIENDLY RIVALS ABOVE THE SWITCH LINE:  ALT cap Riv+ = ALT cap Riv- = empty")
bad = 0
for p in itertools.product([Fraction(i, 4) for i in range(5)], repeat=4):
    pf = [float(x) for x in p]
    if in_alt(p) and (m1lib.in_riv_plus(pf) or m1lib.in_riv_minus(pf)): bad += 1
print("  5^4 grid: points that are simultaneously alternator and rival:", bad)
print("  reason: Riv+ forces p4 = 0 but ALT forces p4 > 0;  Riv- forces p1 = 1 but ALT forces p1 < 1.")
