"""
The FULL-SPACE fixed-eps Nash layer: all 65536 residents, every cell.

    python3 ne16.py build [nproc]     ~5 min on 10 cores, once
    python3 ne16.py check             self-test

WHAT THIS ANSWERS.  For a game (u, v) and an error rate eps, which of the
65536 memory-two strategies are Nash equilibria against ALL 65536 deviations,
and which are STRICT ones?  `Figures/nashe4.py` answers the same question but
only for a handful of residents per cell (the top ten), because it runs the
exact test one Python integer at a time.  This module answers it for the whole
strategy space, which is what a MASS is: sum over sigma in N of pi(sigma)
needs a verdict for every sigma with pi(sigma) > 0, and at u = 1e-2 every one
of the 65536 has pi(sigma) > 0.

THE SIXTEEN ONE-FLIP DEVIATIONS ARE ENOUGH, AT eps > 0 -- FOR THE WEAK TEST.
CL/BinM2NE1 SI Proposition prop:sixteen (Final/PNAS/ms.tex:1571-1598) proves
that at 0 < eps < 1/2 the supremum of the deviator's long-run average payoff
over ALL behavioural deviations -- history-dependent, randomised, unbounded
memory, SUBJECT TO THE SAME IMPLEMENTATION ERROR -- equals the maximum over the
2^16 binary strategies, and that no behavioural deviation improves iff no
binary one does iff none of the sixteen one-flips does.  The mechanism: the
deviator faces a finite average-reward MDP on the sixteen states; at eps > 0
every policy's chain is irreducible, so the MDP is unichain and the Bellman
condition is state-wise.  This FAILS in the eps -> 0 limit, where the chain can
be multichain -- `Opt/exact.py` is the limit instrument and is not this one.

STRICTNESS IS A COROLLARY THE SI DOES NOT STATE.  The proposition is entirely
non-strict.  The extra step, written out because nothing we can cite writes it:
if all sixteen one-flip gains are strictly negative then, by the proposition's
own collapse g(q) - g(p) = sum_s nu_q(s) phi_p(s, q(s)) applied to each
one-flip, every Bellman slack phi_p(s, not p(s)) is strictly negative; and for
ANY stationary policy q != p the same identity gives a sum of terms
nu_q(s) phi_p(s, a) with nu_q(s) > 0 at EVERY state by irreducibility, so
g(q) < g(p).  Hence p is the unique average-optimal stationary policy and all
65535 other binary strategies do strictly worse.  Conversely one gain of
exactly 0 exhibits a distinct policy of equal payoff.  So

    weak  Nash   <=>  all sixteen  a + b u + d v <= 0
    STRICT Nash  <=>  all sixteen  a + b u + d v <  0

WHAT "STRICT" MEANS HERE, EXACTLY.  Strict against the 65536 eps-perturbed
binary strategies -- which is what panel d needs, and what the literature means
by a strict equilibrium of this game.  It does NOT extend to the behavioural
class the weak half of the proposition covers: a deviator that differs from p
on finitely many histories and then copies it forever earns exactly g(p), so
NO resident is ever strict against behavioural deviations, at any eps.  The
weak statement upgrades to "all behavioural"; the strict one does not.

THE TRIPLES DECIDE A SIGN AND NOTHING ELSE.  Clearing to primitive integers
divides each half-plane by its own positive factor, and those factors differ by
up to 1.5e23 WITHIN one resident, so a + b u + d v is a valid sign test and a
meaningless magnitude.  Nor does the one-flip maximum bound the full-space one:
an independent brute force over all 65535 deviations found the argmax is not a
one-flip in 86.8% of (resident, game) pairs, and at resident 22 on game
V_p000_PDW the largest one-flip gain is +1.17e-12 while a 13-flip mutant gains
+0.286.  The SIGNS agree -- that is the proposition, and it held at every one
of 175663 pairs tested -- but any float test with an absolute tolerance coarser
than ~1e-12 would call resident 22 an equilibrium there.  That is what the
exact arithmetic below is for.

THE TRIPLES ARE CL/BinM2NE1's, UNMODIFIED.  `exactcoef.exact16(p, eps)` returns
each half-plane cleared to a primitive integer triple.  They are big: measured over all
65536 x 16 x 3 of them at eps = 1/10000, the largest is 54 decimal digits (code
15060, flip 1), so int64 overflows and float64 cannot decide a tie -- though
every one converts to a finite double, the largest being 5.86e53, which is what
makes the float screen below legitimate.  None of the 1048576 triples is
identically zero, so every half-plane is a real constraint.  Building all 65536
costs ~42 ms each, hence the cache.

HOW A SIGN IS DECIDED, AND WHY IT IS STILL EXACT.  Evaluating 512 cells x
65536 residents x 16 half-planes in Python integers is 5.4e8 big-integer
operations and does not finish.  So each sign is first screened in float64,
and only the entries the screen cannot certify are redone in integers:

    g = a + b*u + d*v            in float64
    S = |a| + |b*u| + |d*v|      the same sum with every term positive

Each of the three products carries at most one rounding of the coefficient and
one of the product, and the two additions one each, so the absolute error is
below 4 * 2^-53 * S = 4.5e-16 * S.  The screen accepts a sign only when
|g| > 1e-12 * S, a margin of about 2250, and refers everything else -- ties,
near-ties, and any all-zero triple -- to exact integers.  Every (u, v) here is
an IEEE double and therefore a dyadic rational exactly, so

    a + b u + d v <= 0   <=>   a*D + b*un*(D/ud) + d*vn*(D/vd) <= 0,  D = ud*vd

with every quantity a Python int and D > 0.

BE HONEST ABOUT WHAT THE INTEGER BRANCH ACTUALLY DOES HERE.  On the 512
sunflower games it never runs: of all 512 x 65536 x 16 signs, ZERO fall inside
the screen, the smallest |g|/S anywhere being 5.27e-09 -- four orders of
magnitude clear of the 1e-12 gate and 1.2e7 times the error bound.  So every
verdict this module hands the figure is decided in float64, and it is exact
because the bound above certifies each sign, not because integers were used.
The integer branch is a guarantee for grids that do land on the tie stratum
(off-grid games do: see `_check`), and it was exercised separately -- the 415
tightest entries and 400000 random ones redone in integers, 0 sign
disagreements.  `nrefer()` reports the count so this stays visible.

NOR IS THE TIE STRATUM MERELY UNVISITED BY LUCK.  The games are exact multiples
of 1/1024, so an exact facet hit is possible rather than measure-zero.  There
is not one exactly-zero gain in the 536870912, and not one identically-zero
triple in the 1048576, so WEAK AND STRICT COINCIDE AS SETS in all 512 cells.
`verdicts` returns both anyway, because that is the cheap check that they do.
"""
import gzip
import os
import pickle
import sys
from fractions import Fraction as F
from math import gcd

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINM2NE1 = os.environ.get(
    "BINM2NE1_EXACT",
    "/Users/martin/Documents/CL/BinM2NE1/Final/Github/exact")
NC = 65536
NSH = 64                                  # 1024 residents per shard


def cachedir(en=1, ed=10000):
    """one directory per error rate, and DISTINCT for distinct rates.

    The first version named it by round(-log10(en/ed)), which sends 1/5000 and
    3/10000 both to `ne16_e4` and would silently serve coefficients computed at
    the wrong eps.  Only exact negative powers of ten keep the short name (so
    the built cache at 1/10000 is still `ne16_e4`); everything else is named by
    its reduced fraction.
    """
    g = gcd(en, ed)
    en, ed = en // g, ed // g
    k = 0
    d = ed
    while d % 10 == 0:
        d //= 10
        k += 1
    name = "e%d" % k if (en == 1 and d == 1) else "f%d_%d" % (en, ed)
    return os.path.join(HERE, os.pardir, "Data", "ne16_" + name)


def _shard(arg):
    lo, hi, en, ed, out = arg
    sys.path.insert(0, BINM2NE1)
    import exactcoef
    eps = F(en, ed)
    T = [[tuple(int(x) for x in t) for t in exactcoef.exact16(c, eps)]
         for c in range(lo, hi)]
    with gzip.open(out, "wb") as f:
        pickle.dump(T, f, protocol=4)
    return lo, hi


def build(nproc=10, en=1, ed=10000):
    from multiprocessing import Pool
    d = cachedir(en, ed)
    os.makedirs(d, exist_ok=True)
    step = NC // NSH
    jobs = []
    for i in range(NSH):
        out = os.path.join(d, "sh%03d.pkl.gz" % i)
        if os.path.exists(out):
            continue
        jobs.append((i * step, (i + 1) * step, en, ed, out))
    if not jobs:
        print("all %d shards present in %s" % (NSH, d))
        return
    print("building %d of %d shards on %d processes" % (len(jobs), NSH, nproc))
    with Pool(nproc) as p:
        for k, (lo, hi) in enumerate(p.imap_unordered(_shard, jobs)):
            print("  %3d/%d  codes %5d-%5d" % (k + 1, len(jobs), lo, hi - 1),
                  flush=True)


_C = {}


def load(en=1, ed=10000):
    """(A, B, D) object arrays (65536,16) of Python ints, and float64 copies"""
    key = (en, ed)
    if key in _C:
        return _C[key]
    d = cachedir(en, ed)
    step = NC // NSH
    A = np.empty((NC, 16), dtype=object)
    B = np.empty((NC, 16), dtype=object)
    D = np.empty((NC, 16), dtype=object)
    for i in range(NSH):
        with gzip.open(os.path.join(d, "sh%03d.pkl.gz" % i), "rb") as f:
            T = pickle.load(f)
        if len(T) != step:
            raise ValueError("shard %d holds %d codes, expected %d"
                             % (i, len(T), step))
        for j, tr in enumerate(T):
            c = i * step + j
            if len(tr) != 16:
                raise ValueError("code %d has %d half-planes" % (c, len(tr)))
            A[c] = [t[0] for t in tr]
            B[c] = [t[1] for t in tr]
            D[c] = [t[2] for t in tr]
    Af = np.array([[float(x) for x in r] for r in A], dtype=np.float64)
    Bf = np.array([[float(x) for x in r] for r in B], dtype=np.float64)
    Df = np.array([[float(x) for x in r] for r in D], dtype=np.float64)
    _C[key] = (A, B, D, Af, Bf, Df)
    return _C[key]


_NREF = [0]


def nrefer():
    """how many (resident, half-plane) signs the float screen could not
    certify and were redone in integers, since import"""
    return _NREF[0]


def verdicts(u, v, en=1, ed=10000, screen=1e-12):
    """(weak, strict) boolean arrays over all 65536 residents, exact.

    weak[c]   : every one of c's sixteen one-flip gains is <= 0 at (u, v)
    strict[c] : every one is < 0
    """
    A, B, D, Af, Bf, Df = load(en, ed)
    uf, vf = float(u), float(v)
    g = Af + Bf * uf + Df * vf
    S = np.abs(Af) + np.abs(Bf * uf) + np.abs(Df * vf)
    pos = g > screen * S                        # certainly > 0
    neg = g < -screen * S                       # certainly < 0
    todo = np.nonzero(~(pos | neg))             # ties, near-ties, zero triples
    if len(todo[0]):
        U, V = F(uf), F(vf)
        un, ud = U.numerator, U.denominator
        vn, vd = V.numerator, V.denominator
        DD = ud * vd
        cu, cv = un * (DD // ud), vn * (DD // vd)
        for c, k in zip(*todo):
            s = A[c][k] * DD + B[c][k] * cu + D[c][k] * cv
            pos[c, k] = s > 0
            neg[c, k] = s < 0
        _NREF[0] += len(todo[0])
    weak = ~pos.any(axis=1)
    strict = neg.all(axis=1)
    return weak, strict


def _check():
    import time
    A, B, D, Af, Bf, Df = load()
    # the numbers ../FiguresFinal/CODE_AND_DATA.md and F1_legend.tex quote
    # about the coefficients are asserted HERE, over all of them, so a rebuilt
    # cache that changed them would be caught rather than believed
    mx, arg, zero = 0, None, 0
    for c in range(NC):
        for k in range(16):
            a, b, d = A[c][k], B[c][k], D[c][k]
            if a == 0 and b == 0 and d == 0:
                zero += 1
            for x in (a, b, d):
                v = abs(int(x))
                if v > mx:
                    mx, arg = v, (c, k)
    print("loaded %d x 16 triples; largest |coefficient| %d digits (code %d, "
          "flip %d); as a double %.4g, finite: %s"
          % (NC, len(str(mx)), arg[0], arg[1], float(mx),
             bool(np.isfinite(Af).all() and np.isfinite(Bf).all()
                  and np.isfinite(Df).all())))
    print("identically-zero triples over all %d: %d" % (NC * 16, zero))

    sys.path.insert(0, os.path.join(HERE, os.pardir, "Figures"))
    import nashe4
    rng = np.random.default_rng(7)
    t = time.time()
    for u, v in ((-0.091796875, 0.083984375), (1.0, 1.0), (-3.5, 2.25),
                 (0.5, 0.5), (-27.0, 14.6875)):
        w, s = verdicts(u, v)
        codes = sorted(set(rng.integers(0, NC, 40).tolist())
                       | set(np.nonzero(w)[0][:15].tolist()))
        ref = nashe4.is_nash(u, v, codes)
        bad = [c for c, r in zip(codes, ref) if bool(w[c]) != bool(r)]
        print("  (%9.4f,%9.4f)  weak %5d  strict %5d   nashe4 agrees on "
              "%d/%d %s" % (u, v, w.sum(), s.sum(), len(codes) - len(bad),
                            len(codes), "" if not bad else "MISMATCH %s" % bad))
    print("%d signs referred to integers; %.1f s" % (nrefer(), time.time() - t))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "build":
        build(int(sys.argv[2]) if len(sys.argv) > 2 else 10)
    else:
        _check()
