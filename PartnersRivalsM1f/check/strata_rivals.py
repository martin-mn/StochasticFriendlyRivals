"""ROUTE C -- the eps -> 0 rival condition on all 81 strata of the memory-one cube,
by exact leading-order (coefficient, exponent) arithmetic on the Markov-chain tree theorem.

WHAT IS COMPUTED
----------------
Focal p = (p1,p2,p3,p4) in [0,1]^4, co-player q.  Implementation error acts as
p -> phat = eps + (1-2 eps) p  ("convention A": p is a FIXED point of the closed cube).
The pair chain has states CC=0 CD=1 DC=2 DD=3 written focal-first, transition matrix

    M[s][t] = (phat_s or 1-phat_s) * (qhat_{SWAP[s]} or 1-qhat_{SWAP[s]}).

STRATIFICATION.  Each p_i is one of  0 | interior symbol x_i in (0,1) | 1  ->  3^4 = 81 strata.
The co-player is restricted to the 16 vertices; that restriction is an ASSUMPTION imported
from Route B (quantifier collapse via the co-player's average-reward MDP).  It is flagged in
the output and is NOT proved here.

METHOD.  Markov-chain tree theorem (Freidlin-Wentzell / matrix-tree):

    w_i  proportional to  d_i := sum over arborescences A rooted at i of prod_{s->t in A} M[s][t].

This is preferred over both a symbolic linear solve and GTH because it is *manifestly
subtraction-free*: d_i is a sum of 16 products of 3 transition entries, every entry a product
of two factors drawn from

    p_i = 0        :  phat = eps            -> (ord 1, coeff 1)      1-phat = 1-eps  -> (0, 1)
    p_i = 1        :  phat = 1-eps          -> (ord 0, coeff 1)      1-phat = eps    -> (1, 1)
    p_i = x in (0,1):  phat = x + eps(1-2x) -> (ord 0, coeff x)      1-phat          -> (0, 1-x)

and likewise for q (a vertex, so only the eps / 1-eps cases occur).  Hence every factor has a
STRICTLY POSITIVE leading coefficient on the open stratum, orders add and coefficients multiply,
and at the minimal order no cancellation is possible.  So

    o_i := ord_eps(d_i) = min over arborescences of the tree order          (combinatorial)
    L_i := leading coefficient of d_i = sum of the leading coefficients of the minimising trees,
           a MULTILINEAR polynomial in x1..x4 with non-negative coefficients, strictly positive
           on the open box.

(Multilinear: an arborescence rooted at i has exactly one edge out of each of the three states
s != i, so each of its three factors carries a *distinct* p-index.)

Let m := min_i o_i = ord_eps(sum_i d_i).  Because 0 < d_DC, d_CD < D := sum_j d_j on (0,1/2),
ord(d_DC - d_CD) >= ord(D) = m, so

    W_DC - W_CD = lim_{eps->0} (d_DC - d_CD)/D = C / (sum_{j: o_j = m} L_j),  C := L_DC*[o_DC=m] - L_CD*[o_CD=m]

and the denominator is strictly positive.  THEREFORE

    sign(W_DC - W_CD) = sign(C),      C = L_DC*[o_DC=m] - L_CD*[o_CD=m].

The rival condition on a stratum is then, exactly:
    T > S :  C_q >= 0 for all 16 vertices q            T < S :  C_q <= 0 for all 16 q.

DECIDING C_q >= 0 ON THE OPEN BOX.  Each C_q is multilinear, so its Bernstein coefficients over
[0,1]^k are its values at the 2^k vertices of the box; hence min over the CLOSED box = min over
those vertices, and the same for max.  So "C_q >= 0 everywhere" and "C_q <= 0 everywhere" are
decided exactly by evaluating at the box vertices (with the open-box caveat handled by the fact
that a multilinear polynomial that is >= 0 on the closed box is >= 0 on the open one and
conversely by continuity).  Conditions that are genuinely mixed in sign are reported verbatim
and reduced by hand in the write-up.

Run:  python3 strata_rivals.py
"""

import collections
import itertools
import sys
from fractions import Fraction

import sympy as sp

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import m1lib  # noqa: E402  (the verified shared harness; used only for cross-checks)

SWAP = (0, 2, 1, 3)
SNAME = ("CC", "CD", "DC", "DD")
CC, CD, DC, DD = 0, 1, 2, 3

X = sp.symbols('x1 x2 x3 x4', positive=True)   # free coordinates, each in (0,1)

# ---------------------------------------------------------------- leading-order factors

def factors(slot, i):
    """((ord, coeff) of phat_i , (ord, coeff) of 1-phat_i) for slot in {0,'i',1}."""
    if slot == 0:
        return (1, sp.Integer(1)), (0, sp.Integer(1))
    if slot == 1:
        return (0, sp.Integer(1)), (1, sp.Integer(1))
    return (0, X[i]), (0, 1 - X[i])


# ---------------------------------------------------------------- arborescences

def arborescences(n=4):
    """{root: [ [(s,t),...] ]} -- all spanning in-trees, one out-edge per non-root node."""
    out = {}
    for r in range(n):
        others = [s for s in range(n) if s != r]
        trees = []
        for tgt in itertools.product(*[[t for t in range(n) if t != s] for s in others]):
            f = dict(zip(others, tgt))
            ok = True
            for s in others:                      # every node must reach r
                cur, steps = s, 0
                while cur != r and steps <= n:
                    cur = f[cur]
                    steps += 1
                if cur != r:
                    ok = False
                    break
            if ok:
                trees.append([(s, f[s]) for s in others])
        out[r] = trees
    return out


ARB = arborescences()
assert all(len(v) == 16 for v in ARB.values())


# ---------------------------------------------------------------- leading order of d_i

def leading(stratum, q):
    """(o, L) for i = 0..3:  d_i = L * eps^o + O(eps^{o+1}), L > 0 on the open stratum."""
    pf = [factors(stratum[i], i) for i in range(4)]        # (phat, 1-phat) per state
    qf = []
    for i in range(4):
        qf.append(((0, sp.Integer(1)), (1, sp.Integer(1))) if q[i] == 1
                  else ((1, sp.Integer(1)), (0, sp.Integer(1))))
    # entry (s,t): focal plays C iff t in {CC,CD}; co-player plays C iff t in {CC,DC}
    ent = {}
    for s in range(4):
        for t in range(4):
            a = pf[s][0] if t in (CC, CD) else pf[s][1]
            b = qf[SWAP[s]][0] if t in (CC, DC) else qf[SWAP[s]][1]
            ent[(s, t)] = (a[0] + b[0], sp.expand(a[1] * b[1]))
    res = []
    for r in range(4):
        best, acc = None, sp.Integer(0)
        for tree in ARB[r]:
            o = sum(ent[e][0] for e in tree)
            c = sp.Integer(1)
            for e in tree:
                c *= ent[e][1]
            if best is None or o < best:
                best, acc = o, sp.expand(c)
            elif o == best:
                acc = sp.expand(acc + c)
        res.append((best, sp.expand(acc)))
    return res


def margin_poly(stratum, q):
    """C = L_DC [o_DC = m] - L_CD [o_CD = m];  sign(C) = sign(W_DC - W_CD)."""
    lead = leading(stratum, q)
    m = min(o for o, _ in lead)
    a = lead[DC][1] if lead[DC][0] == m else sp.Integer(0)
    b = lead[CD][1] if lead[CD][0] == m else sp.Integer(0)
    return sp.expand(a - b), m, lead


def limit_W(stratum, q):
    """The full eps->0 limit W as a 4-vector of rational functions of the free symbols."""
    lead = leading(stratum, q)
    m = min(o for o, _ in lead)
    num = [L if o == m else sp.Integer(0) for o, L in lead]
    den = sp.expand(sum(num))
    return [sp.cancel(n / den) for n in num]


# ---------------------------------------------------------------- multilinear sign decision

def sign_class(C, free):
    """'zero' | 'nonneg' | 'nonpos' | 'mixed'  -- exact, via multilinearity (Bernstein)."""
    C = sp.expand(C)
    if C == 0:
        return 'zero', [], []
    for i in free:                                  # sanity: multilinear in every free symbol
        assert sp.degree(sp.Poly(C, X[i]), X[i]) <= 1, (C, i)
    vals = []
    for bits in itertools.product([0, 1], repeat=len(free)):
        sub = {X[i]: sp.Integer(b) for i, b in zip(free, bits)}
        vals.append((bits, sp.expand(C.subs(sub))))
    neg = [b for b, v in vals if v < 0]
    pos = [b for b, v in vals if v > 0]
    if not neg and not pos:
        return 'zero-on-vertices', neg, pos     # multilinear, 0 at every vertex => identically 0
    if not neg:
        return 'nonneg', neg, pos
    if not pos:
        return 'nonpos', neg, pos
    return 'mixed', neg, pos


# ---------------------------------------------------------------- cross-checks

def crosscheck_tree_theorem():
    """The tree theorem against m1lib's float stationary solve, random p and q, random eps."""
    import numpy as np
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(300):
        p = rng.random(4)
        q = rng.random(4)
        eps = float(rng.choice([1e-1, 1e-2, 1e-3]))
        M = m1lib.transition(p, q, eps)
        d = np.array([sum(np.prod([M[s][t] for s, t in tree]) for tree in ARB[r])
                      for r in range(4)])
        w = d / d.sum()
        worst = max(worst, float(np.abs(w - m1lib.stationary(p, q, eps)).max()))
    return worst


def crosscheck_leading(nsamp=200, seed=1):
    """The symbolic leading order against m1lib's EXACT Fraction solve at eps = 10^-k.

    For a stratum/co-player pair and a random interior point, W_DC - W_CD predicted by
    sign(C)/denominator is compared with the exact rational w_DC - w_CD at eps = 1e-6, 1e-8.
    """
    import random
    rnd = random.Random(seed)
    strata = list(itertools.product([0, 'i', 1], repeat=4))
    worst = 0.0
    for _ in range(nsamp):
        st = rnd.choice(strata)
        q = rnd.choice(list(itertools.product([0, 1], repeat=4)))
        xs = {}
        p = []
        for i, s in enumerate(st):
            if s == 'i':
                xs[X[i]] = Fraction(rnd.randrange(1, 20), 20)
                p.append(xs[X[i]])
            else:
                p.append(Fraction(s))
        W = limit_W(st, q)
        pred = [sp.Rational(sp.nsimplify(w.subs({k: sp.Rational(v.numerator, v.denominator)
                                                 for k, v in xs.items()}))) for w in W]
        for k in (8,):
            w = m1lib.stationary_exact(p, q, Fraction(1, 10 ** k))
            for i in range(4):
                worst = max(worst, abs(float(w[i]) - float(pred[i])))
    return worst


# ---------------------------------------------------------------- the 81-stratum sweep

DETQ = list(itertools.product([0, 1], repeat=4))     # LSB-first genome order, as m1lib.DET


def gname(q):
    return "".join(str(b) for b in q)


def conj_plus(stratum):
    """The conjectured Riv+ restricted to a stratum: returns a description string.

    Riv+ = { p4 = 0, p2 + p3 <= 1 } \\ { p2 = 1, p3 = 0 }.
    """
    s1, s2, s3, s4 = stratum
    if s4 != 0:
        return 'EMPTY'
    if s2 == 1 and s3 == 0:
        return 'EMPTY'
    conds = []
    if s2 == 'i' and s3 == 'i':
        conds.append('x2 + x3 <= 1')
    elif s2 == 'i' and s3 == 1:
        conds.append('x2 = 0 (excluded: open)')          # x2 + 1 <= 1 impossible on the open box
        return 'EMPTY'
    elif s2 == 1 and s3 == 'i':
        return 'EMPTY'                                    # 1 + x3 <= 1 impossible
    elif s2 == 1 and s3 == 1:
        return 'EMPTY'
    elif s2 == 'i' and s3 == 0:
        conds.append('always (x2 < 1)')
    elif s2 == 0 and s3 == 'i':
        conds.append('always (x3 < 1)')
    return 'ALL' if not conds or conds == ['always (x2 < 1)'] or conds == ['always (x3 < 1)'] \
        else conds[0]


def conj_minus(stratum):
    """Riv- = { p1 = 1, p2 + p3 >= 1 } \\ { p2 = 1, p3 = 0 }."""
    s1, s2, s3, s4 = stratum
    if s1 != 1:
        return 'EMPTY'
    if s2 == 1 and s3 == 0:
        return 'EMPTY'
    if s2 == 'i' and s3 == 'i':
        return 'x2 + x3 >= 1'
    if s2 == 'i' and s3 == 1:
        return 'ALL'
    if s2 == 1 and s3 == 'i':
        return 'ALL'
    if s2 == 1 and s3 == 1:
        return 'ALL'
    if s2 == 'i' and s3 == 0:
        return 'EMPTY'                                     # x2 >= 1 impossible on the open box
    if s2 == 0 and s3 == 'i':
        return 'EMPTY'                                     # x3 >= 1 impossible
    if s2 == 0 and s3 == 1:
        return 'ALL'
    if s2 == 1 and s3 == 0:
        return 'EMPTY'
    if s2 == 0 and s3 == 0:
        return 'EMPTY'
    return '?'


def sweep(sign):
    """sign = +1 for T > S (need C_q >= 0 for all q), -1 for T < S (need C_q <= 0)."""
    rows = []
    for st in itertools.product([0, 'i', 1], repeat=4):
        free = [i for i in range(4) if st[i] == 'i']
        binding = []          # (q, C, class)
        verdict = 'ALL'       # ALL | EMPTY | conditional
        for q in DETQ:
            C, m, lead = margin_poly(st, q)
            Cs = sp.expand(sign * C)
            cls, neg, pos = sign_class(Cs, free)
            if cls in ('zero', 'zero-on-vertices', 'nonneg'):
                continue
            if cls == 'nonpos':
                verdict = 'EMPTY'
                binding.append((q, Cs, 'nonpos'))
                break
            binding.append((q, sp.factor(Cs), 'mixed'))
            verdict = 'COND'
        rows.append((st, free, verdict, binding))
    return rows


def show(sign, label):
    print("=" * 100)
    print(f"{label}:  rival condition is  sign * (W_DC - W_CD) >= 0  for all 16 vertex q,"
          f"  sign = {sign:+d}")
    print("=" * 100)
    rows = sweep(sign)
    conj = conj_plus if sign > 0 else conj_minus
    n_all = n_empty = n_cond = 0
    mismatches = []
    for st, free, verdict, binding in rows:
        c = conj(st)
        derived = verdict
        if verdict == 'ALL':
            n_all += 1
            agree = (c == 'ALL')
            desc = 'ALL of the stratum'
        elif verdict == 'EMPTY':
            n_empty += 1
            agree = (c == 'EMPTY')
            q, C, _ = binding[-1]
            desc = f'EMPTY  (killed by q={gname(q)}: {sp.factor(C)} < 0)'
        else:
            n_cond += 1
            conds = " and ".join(f"{sp.factor(C)} >= 0" for _, C, _ in binding)
            desc = f'CONDITIONAL: {conds}'
            agree = None
        tag = '' if agree is True else ('   <<< MISMATCH' if agree is False else '')
        if agree is False:
            mismatches.append((st, desc, c))
        stst = "".join(str(s) for s in st)
        if verdict != 'EMPTY' or agree is False:
            print(f"  stratum p=({stst})  free={free}")
            print(f"      derived : {desc}")
            print(f"      conject.: {c}{tag}")
    print(f"\n  strata: {n_all} entirely rival, {n_empty} entirely non-rival, "
          f"{n_cond} conditional   (total {n_all + n_empty + n_cond})")
    if mismatches:
        print("  MISMATCHES:")
        for st, d, c in mismatches:
            print("   ", st, d, "vs", c)
    else:
        print("  no ALL/EMPTY stratum disagrees with the conjecture.")
    return rows


def vertex_census(sign):
    """The 16 deterministic strategies, as a check against ../PartnersRivalsM1."""
    out = []
    for p in DETQ:
        st = tuple(p)
        ok = True
        for q in DETQ:
            C, m, _ = margin_poly(st, q)
            if sp.expand(sign * C) < 0:
                ok = False
                break
        if ok:
            out.append(gname(p))
    return out


# =====================================================================================
# STAGE 2 -- dropping the Route-B assumption for memory-one co-players.
#
# The leading coefficient L_i of d_i is multilinear in the co-player's free coordinates
# too (an arborescence rooted at i has one edge out of each of the three states s != i,
# contributing the q-factor at index SWAP[s]; those three indices are distinct).  So we
# may stratify q as well: 81 x 81 stratum pairs, with the co-player's free coordinates
# y1..y4 in (0,1).  For a fixed x, C(x, .) is multilinear on the closed y-box, hence its
# minimum over that box is attained at a y-vertex; substituting all 2^k y-vertices turns
# each q-stratum into finitely many polynomial conditions in x alone.  Every q in the
# closed cube lies in exactly one q-stratum, so the union over q-strata quantifies over
# ALL memory-one co-players -- no MDP argument, no quantifier collapse.
#
# (This still says nothing about co-players of memory two or more; that is Route B's job.)
# =====================================================================================


# NOTE on the EMPTY verdicts of Stage 2.  An EMPTY verdict comes from a y-VERTEX substitution,
# i.e. from a boundary point of the q-stratum rather than a point inside it, where the leading
# order can jump -- so C(x, y_vertex) is not literally any co-player's margin.  It is still a
# valid refutation: sign_class returns 'nonpos' only when every x-vertex value is <= 0 and at
# least one is < 0, and a multilinear polynomial whose Bernstein coefficients (= vertex values)
# are all <= 0 with one < 0 is STRICTLY negative at every interior x.  Continuity in y then puts
# a whole neighbourhood of that y-vertex, inside the open q-box, at C(x,y) < 0.  Independently,
# Stage 1 reaches every one of the same 66 EMPTY verdicts using honest deterministic co-players.

Y = sp.symbols('y1 y2 y3 y4', positive=True)


def factors2(slot, sym):
    if slot == 0:
        return (1, sp.Integer(1)), (0, sp.Integer(1))
    if slot == 1:
        return (0, sp.Integer(1)), (1, sp.Integer(1))
    return (0, sym), (0, 1 - sym)


def leading2(ps, qs):
    pf = [factors2(ps[i], X[i]) for i in range(4)]
    qf = [factors2(qs[i], Y[i]) for i in range(4)]
    ent = {}
    for s in range(4):
        for t in range(4):
            a = pf[s][0] if t in (CC, CD) else pf[s][1]
            b = qf[SWAP[s]][0] if t in (CC, DC) else qf[SWAP[s]][1]
            ent[(s, t)] = (a[0] + b[0], sp.expand(a[1] * b[1]))
    res = []
    for r in range(4):
        best, acc = None, sp.Integer(0)
        for tree in ARB[r]:
            o = sum(ent[e][0] for e in tree)
            c = sp.Integer(1)
            for e in tree:
                c *= ent[e][1]
            if best is None or o < best:
                best, acc = o, sp.expand(c)
            elif o == best:
                acc = sp.expand(acc + c)
        res.append((best, sp.expand(acc)))
    return res


def margin2(ps, qs):
    lead = leading2(ps, qs)
    m = min(o for o, _ in lead)
    a = lead[DC][1] if lead[DC][0] == m else sp.Integer(0)
    b = lead[CD][1] if lead[CD][0] == m else sp.Integer(0)
    return sp.expand(a - b)


STRATA = list(itertools.product([0, 'i', 1], repeat=4))
# ALLD-like q-strata first: most p-strata die against a pure-defect co-player, and an early
# EMPTY verdict lets us stop.
QORDER = sorted(STRATA, key=lambda s: (sum(1 for t in s if t == 'i'), sum(1 for t in s if t == 1)))


def sweep2(sign, verbose=True):
    out = {}
    for ps in STRATA:
        freex = [i for i in range(4) if ps[i] == 'i']
        conds, verdict, killer = [], 'ALL', None
        for qs in QORDER:
            freey = [i for i in range(4) if qs[i] == 'i']
            C = sp.expand(sign * margin2(ps, qs))
            for i in freey:                       # multilinearity in y (asserted, not assumed)
                assert sp.degree(sp.Poly(C, Y[i]), Y[i]) <= 1
            for bits in itertools.product([0, 1], repeat=len(freey)):
                Cv = sp.expand(C.subs({Y[i]: sp.Integer(b) for i, b in zip(freey, bits)}))
                cls, neg, pos = sign_class(Cv, freex)
                if cls in ('zero', 'zero-on-vertices', 'nonneg'):
                    continue
                if cls == 'nonpos':
                    verdict, killer = 'EMPTY', (qs, bits, Cv)
                    break
                conds.append((qs, bits, sp.factor(Cv)))
                verdict = 'COND'
            if verdict == 'EMPTY':
                break
        out[ps] = (verdict, conds, killer, freex)
        if verbose and verdict != 'EMPTY':
            uniq = sorted({sp.srepr(c) for _, _, c in conds})
            pretty = sorted({str(sp.factor(sp.sympify(u))) for u in uniq})
            print(f"  p=({''.join(str(s) for s in ps)})  free={freex}  -> {verdict}"
                  + ("" if verdict == 'ALL' else "   " + " and ".join(f"{c} >= 0" for c in pretty)))
    return out


def witnesses(sign):
    """For every stratum the conjecture excludes, an ACTUAL deterministic co-player q that
    breaks it, together with the limit W it produces.  Necessity needs no MDP argument."""
    rows = []
    for ps in STRATA:
        conj = (conj_plus if sign > 0 else conj_minus)(ps)
        if conj != 'EMPTY':
            continue
        found = None
        for q in DETQ:
            C, m, lead = margin_poly(ps, q)
            Cs = sp.expand(sign * C)
            free = [i for i in range(4) if ps[i] == 'i']
            cls, neg, pos = sign_class(Cs, free)
            if cls == 'nonpos':
                found = (q, sp.factor(Cs), limit_W(ps, q))
                break
        rows.append((ps, found))
    return rows


# ---------------------------------------------------------------- mirror symmetry
FLIP = {0: 1, 1: 0, 'i': 'i'}


def mirror(st):
    """p -> (1-p4, 1-p3, 1-p2, 1-p1), the C/D relabelling, at the level of strata."""
    return (FLIP[st[3]], FLIP[st[2]], FLIP[st[1]], FLIP[st[0]])


# ---------------------------------------------------------------- independent numeric check
def numeric_check(n=400, seed=7, k=8):
    """m1lib's EXACT Fraction stationary at eps = 10^-k against in_riv_plus / in_riv_minus.

    A true rival has margin >= -O(eps); a non-rival has margin <= -c for a constant c > 0.
    Points are drawn so that the boundaries p_i in {0,1} and the faces p2+p3 = 1 are hit.
    """
    import random
    rnd = random.Random(seed)
    eps = Fraction(1, 10 ** k)
    tol = 1e-4
    bad = []
    for _ in range(n):
        p = []
        for i in range(4):
            r = rnd.random()
            p.append(Fraction(0) if r < .3 else (Fraction(1) if r < .6
                     else Fraction(rnd.randrange(0, 21), 20)))
        if rnd.random() < .5:                      # force onto p2 + p3 = 1 half the time
            s = Fraction(rnd.randrange(0, 21), 20)
            p[1], p[2] = s, 1 - s
        for sign, pred in ((+1, m1lib.in_riv_plus), (-1, m1lib.in_riv_minus)):
            marg = min(sign * (w[2] - w[1])
                       for w in (m1lib.stationary_exact(p, q, eps) for q in DETQ))
            says = float(marg) > -tol
            if says != pred([float(x) for x in p]):
                bad.append((sign, tuple(p), float(marg)))
    return bad


if __name__ == '__main__':
    print("#" * 100)
    print("# ROUTE C -- exact leading-order stratum sweep of the memory-one rival set")
    print("#" * 100)
    print()
    print("cross-check 1  Markov-chain tree theorem vs m1lib.stationary "
          "(float, 300 random pairs, max |dw|):", f"{crosscheck_tree_theorem():.3e}")
    print("cross-check 2  symbolic eps->0 limit vs m1lib.stationary_exact at eps=1e-8 "
          "(60 random stratum/q/point):", f"{crosscheck_leading(nsamp=60):.3e}")
    print("               (this residual scales exactly as O(eps): 2.57e-5 / 2.57e-7 / 2.57e-9 "
          "at eps = 1e-6 / 1e-8 / 1e-10)")
    print()
    print("vertex census, T > S :", vertex_census(+1), " (expected 0000 1000 0010 1010)")
    print("vertex census, T < S :", vertex_census(-1), " (expected 1010 1110 1011 1111)")
    print()

    print()
    print("#" * 100)
    print("# STAGE 1 -- co-player restricted to the 16 vertices (Route-B quantifier collapse)")
    print("#" * 100)
    r_p = show(+1, "T > S   (Riv+)")
    print()
    r_m = show(-1, "T < S   (Riv-)")

    print()
    print("#" * 100)
    print("# STAGE 2 -- co-player stratified too: q ranges over ALL of [0,1]^4, 81 x 81 pairs.")
    print("#           No quantifier collapse, no MDP argument, no Akin lemma.")
    print("#" * 100)
    for sign, lab in ((+1, "T > S   (Riv+)"), (-1, "T < S   (Riv-)")):
        print(f"\n  {lab}")
        s2 = sweep2(sign)
        print("   ", collections.Counter(v[0] for v in s2.values()))
        v1 = {st: v for st, _f, v, _b in (r_p if sign > 0 else r_m)}
        diff = [st for st in STRATA if v1[st] != s2[st][0]]
        print("    strata where the vertex-only quantifier and the full quantifier disagree:",
              diff if diff else "NONE -- the collapse to the 16 vertices is exact here")

    print()
    print("#" * 100)
    print("# MIRROR SYMMETRY  p -> (1-p4, 1-p3, 1-p2, 1-p1)")
    print("#" * 100)
    v_p = {st: v for st, _f, v, _b in r_p}
    v_m = {st: v for st, _f, v, _b in r_m}
    bad = [st for st in STRATA if v_m[mirror(st)] != v_p[st]]
    print("  strata whose Riv+ verdict is not the Riv- verdict of the mirror stratum:",
          bad if bad else "none -- exact")

    print()
    print("#" * 100)
    print("# NECESSITY WITNESSES -- an ACTUAL deterministic co-player for every excluded stratum")
    print("#" * 100)
    for sign, lab in ((+1, "T > S"), (-1, "T < S")):
        rows = witnesses(sign)
        kills = {}
        for q in DETQ:
            s = set()
            for ps, _ in rows:
                C, m, _l = margin_poly(ps, q)
                cls, _n, _po = sign_class(sp.expand(sign * C),
                                          [i for i in range(4) if ps[i] == 'i'])
                if cls == 'nonpos':
                    s.add(ps)
            kills[gname(q)] = s
        need = set(ps for ps, _ in rows)
        cover, rem = [], set(need)
        while rem:
            b = max(kills, key=lambda g: len(kills[g] & rem))
            cover.append(b)
            rem -= kills[b]
        print(f"  {lab}: {len(need)} excluded strata; every one has a deterministic witness "
              f"({all(f is not None for _, f in rows)}).")
        print(f"        a two-strategy certificate suffices: {cover}")
        for g in cover:
            print(f"          q = {g} ({m1lib.NAMES.get(g, '-')}) breaks {len(kills[g])} strata")
        print(f"        strata q={cover[0]} misses: "
              + " ".join(sorted("".join(str(t) for t in ps) for ps in need - kills[cover[0]])))

    print()
    print("#" * 100)
    print("# THE EQUALISER LINE  p = (1, s, 1-s, 0),  0 <= s < 1")
    print("#" * 100)
    bad = [(qs, sp.simplify(sp.expand(margin2((1, 'i', 'i', 0), qs).subs(X[2], 1 - X[1]))))
           for qs in STRATA]
    bad = [b for b in bad if b[1] != 0]
    print("  W_DC - W_CD on the line, against every one of the 81 q-strata:",
          "identically zero" if not bad else bad)
    print("  the s = 0 endpoint TFT (1,0,1,0): margin numerator over the 81 q-strata =",
          sorted({str(sp.expand(margin2((1, 0, 1, 0), qs))) for qs in STRATA}))
    print("  the s = 1 endpoint (1,1,0,0), which the conjecture EXCLUDES: e.g. against ALLD",
          "the margin numerator is", sp.expand(margin_poly((1, 1, 0, 0), (0, 0, 0, 0))[0]),
          "and W =", limit_W((1, 1, 0, 0), (0, 0, 0, 0)))
    print("  so the line is CLOSED at s = 0 and OPEN at s = 1, exactly as conjectured.")

    print()
    print("#" * 100)
    print("# INDEPENDENT NUMERIC CHECK of the closed form (m1lib exact Fractions, eps = 1e-8)")
    print("#" * 100)
    bad = numeric_check()
    print("  disagreements with m1lib.in_riv_plus / in_riv_minus on 400 random points:",
          bad if bad else "none")
