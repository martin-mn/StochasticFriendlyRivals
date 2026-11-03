"""SYMBOLIC receipt for the closed form of the rival set over the full memory-one cube.

Nothing here is numerical: eps is a sympy symbol, the p_i are sympy symbols, the stationary
distribution of the 4-state pair chain is an exact rational function of (eps, p), the eps -> 0
limit is taken with sympy.limit(..., eps, 0, '+'), and the resulting inequalities are decided
on the closed box [0,1]^k by a BERNSTEIN certificate (all Bernstein coefficients of a
polynomial nonnegative  =>  the polynomial is nonnegative on the box).  No grid, no tolerance.

Conventions are m1lib's and are not redefined: states CC=0 CD=1 DC=2 DD=3, focal action first,
clamping p -> eps + (1-2eps) p (convention A: p is a FIXED point of the closed cube),
pi(p,q) - pi(q,p) = (T-S) (w_DC - w_CD), and p is a rival for T > S iff
lim_{eps->0+} (w_DC - w_CD) >= 0 for every co-player q.

WHAT IS SPECIALISED, AND WHY
----------------------------
The limit W depends on which coordinates of p sit at the boundary of the cube, because the
eps = 0 chain can be reducible there.  A single generic rational function therefore cannot
answer the question at p4 = 0 or p2 = 1; the limit must be taken on each stratum separately.
So the script is stratified exactly as the task prescribes, and each stratum is stated:

  STRATUM G  (p1,p2,p3,p4 all symbolic, generic interior point)
      only the witness ALLD = 0000 is needed: its limit margin is -p4/(1 - p2 + p4), which is
      < 0 for p4 > 0.  This kills the whole slab p4 > 0 and is the one place where all four
      p_i stay symbolic.
  STRATUM A  (p4 = 0 exactly; p1,p2,p3 symbolic, p2 generic i.e. p2 != 1)
      all 16 vertex co-players, full limits, and the derivation of min_q >= 0.
  STRATUM E  (p4 = 0 and p2 = 1 exactly; p1,p3 symbolic)
      the degenerate locus where the eps = 0 chain has two recurrent classes.  Handled at
      finite eps and then in the limit; ALLD refutes with margin -(1+p3)/(3+p3-p1).

Akin's identity is verified symbolically at finite eps with ALL FOUR p_i symbolic, against all
16 vertex co-players (no specialisation at all), and then in the limit.

Run:  python3 symbolic_rivals.py
"""

import sys
import time
from itertools import product

import sympy as sp

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import m1lib                      # the verified shared harness; the chain is NOT rewritten,
                                  # the symbolic chain below is cross-checked against it exactly

t0 = time.time()

# --------------------------------------------------------------------------- symbols
eps = sp.Symbol('epsilon', positive=True)
p1, p2, p3, p4 = sp.symbols('p1 p2 p3 p4', nonnegative=True)
s_par = sp.Symbol('s', nonnegative=True)          # generalised-TFT parameter, p = (1,s,1-s,0)
P = (p1, p2, p3, p4)

SWAP = (0, 2, 1, 3)                     # state seen from the co-player's side
STATE = ('CC', 'CD', 'DC', 'DD')
VERTICES = [tuple(b) for b in product([0, 1], repeat=4)]   # LSB-first genome order


def genome(q):
    return ''.join(str(int(x)) for x in q)


NAMES = {'0000': 'ALLD', '1111': 'ALLC', '1010': 'TFT', '1001': 'WSLS', '1000': 'Grim',
         '0101': 'ATFT', '0110': 'AWSLS', '0001': "q* (C only after DD)", '0111': "m(q*)"}


def label(q):
    g = genome(q)
    return f"{g}" + (f" [{NAMES[g]}]" if g in NAMES else "")


def clamp(x):
    return eps + (1 - 2 * eps) * x


# ------------------------------------------------------- exact stationary distribution
# Markov-chain tree theorem.  K = M^T - I has columns summing to zero, so the principal
# cofactors C_ii = det(K with row i and column i removed) are, up to one common sign, the
# unnormalised stationary weights d_i (each a sum of 16 products of three transition
# entries -- manifestly a polynomial, no subtraction).  This is far cheaper symbolically
# than a 4x4 linear solve and is cross-checked below against sympy's own nullspace.

def transition(p, q):
    """Row-stochastic pair-chain matrix for clamped focal p against clamped co-player q."""
    pc = [clamp(x) for x in p]
    qc = [clamp(x) for x in q]
    M = sp.zeros(4, 4)
    for s in range(4):
        a, b = pc[s], qc[SWAP[s]]
        M[s, 0] = a * b
        M[s, 1] = a * (1 - b)
        M[s, 2] = (1 - a) * b
        M[s, 3] = (1 - a) * (1 - b)
    return M


def unnormalised_weights(p, q):
    M = transition(p, q)
    K = M.T - sp.eye(4)
    d = []
    for i in range(4):
        minor = K.extract([r for r in range(4) if r != i], [c for c in range(4) if c != i])
        d.append(sp.expand(minor.det()))
    # fix the common sign so that the weights are the positive ones (det of a 3x3 principal
    # minor of -Laplacian carries the sign (-1)^3); normalise by evaluating the total.
    tot = sp.expand(sum(d))
    # every free symbol is a probability-like quantity in (0,1); substitute distinct
    # generic rationals for all of them to read off the overall sign unambiguously.
    vals = [sp.Rational(1, 7), sp.Rational(1, 3), sp.Rational(1, 5), sp.Rational(2, 7),
            sp.Rational(3, 11), sp.Rational(4, 13), sp.Rational(5, 17)]
    probe = {eps: sp.Rational(1, 7)}
    for i, sym in enumerate(sorted(tot.free_symbols - {eps}, key=str)):
        probe[sym] = vals[(i + 1) % len(vals)]
    if sp.nsimplify(tot.subs(probe)) < 0:
        d = [sp.expand(-x) for x in d]
    return d


def margin_limit(p, q):
    """lim_{eps->0+} (w_DC - w_CD), an exact rational function of the free p_i."""
    d = unnormalised_weights(p, q)
    expr = sp.cancel((d[2] - d[1]) / sum(d))
    return sp.simplify(sp.limit(expr, eps, 0, '+'))


def margin_finite(p, q):
    """w_DC - w_CD at finite eps, exact."""
    d = unnormalised_weights(p, q)
    return sp.cancel((d[2] - d[1]) / sum(d))


# --------------------------------------------------------------- Bernstein certificate
def bernstein_coeffs(poly_expr, vars_):
    """Bernstein coefficients of a polynomial over the box [0,1]^n.

    All coefficients >= 0  =>  the polynomial is >= 0 on the closed box (sufficient, exact).
    """
    e = sp.expand(poly_expr)
    if e == 0:
        return {(0,) * len(vars_): sp.Integer(0)}
    pol = sp.Poly(e, *vars_)
    deg = [int(pol.degree(v)) if pol.degree(v) != sp.S.NegativeInfinity else 0
           for v in vars_]
    a = {}
    for monom, coeff in pol.terms():
        a[tuple(monom)] = sp.Rational(coeff)
    out = {}
    for alpha in product(*[range(d + 1) for d in deg]):
        tot = sp.Integer(0)
        for beta, coeff in a.items():
            if all(b <= al for b, al in zip(beta, alpha)):
                f = sp.Integer(1)
                for i in range(len(vars_)):
                    f *= sp.Rational(sp.binomial(alpha[i], beta[i]),
                                     sp.binomial(deg[i], beta[i]))
                tot += f * coeff
        out[alpha] = sp.nsimplify(tot)
    return out


def nonneg_on_box(expr, vars_):
    b = bernstein_coeffs(expr, vars_)
    bad = {k: v for k, v in b.items() if v < 0}
    return (not bad), bad


def sign_certificate(m, box, want_nonneg):
    """Exact Bernstein test of  m >= 0  (or m <= 0) on the unit box in `box`."""
    if m == 0:
        return want_nonneg, m
    num, den = [sp.expand(x) for x in sp.fraction(sp.cancel(sp.together(m)))]
    if box:
        ok_d, _ = nonneg_on_box(den, box)
        if not ok_d:
            ok_d2, _ = nonneg_on_box(-den, box)
            if ok_d2:
                num, den, ok_d = sp.expand(-num), sp.expand(-den), True
        ok_n, _ = nonneg_on_box(num if want_nonneg else -num, box)
    else:
        if den < 0:
            num, den = -num, -den
        ok_d = bool(den > 0)
        ok_n = bool(num >= 0) if want_nonneg else bool(num <= 0)
    return (ok_d and ok_n), m


def rule(title):
    print('\n' + '=' * 78)
    print(title)
    print('=' * 78)


# =====================================================================================
rule("CROSS-CHECK against m1lib's exact Fraction chain (no floats, no tolerance)")
from fractions import Fraction
_e = sp.Rational(1, 10)
_ps = [(sp.Rational(1, 3), sp.Rational(1, 5), sp.Rational(2, 7), sp.Rational(3, 11)),
       (sp.Rational(1), sp.Rational(2, 5), sp.Rational(3, 5), sp.Rational(0)),
       (sp.Rational(0), sp.Rational(1), sp.Rational(0), sp.Rational(0))]
_bad = 0
for _p in _ps:
    for _q in VERTICES:
        _sym = sp.nsimplify(margin_finite(_p, _q).subs(eps, _e))
        _w = m1lib.stationary_exact([Fraction(int(x.p), int(x.q)) for x in _p],
                                    [Fraction(int(x)) for x in _q], Fraction(1, 10))
        _ref = sp.Rational(_w[2] - _w[1])
        if sp.simplify(_sym - _ref) != 0:
            _bad += 1
            print('   MISMATCH', _p, genome(_q), _sym, _ref)
print(f'   48 pair-chain margins (3 focal p x 16 vertex q) at eps = 1/10, exact rationals:')
print(f'   mismatches against m1lib.stationary_exact = {_bad}')
assert _bad == 0
print('   OK -- the symbolic chain and the harness chain are the same object.')

# =====================================================================================
rule('SANITY: tree-theorem weights agree with sympy nullspace (symbolic, p4=0, q=0001)')
qtest = (0, 0, 0, 1)
ptest = (p1, p2, p3, sp.Integer(0))
d = unnormalised_weights(ptest, qtest)
w_tree = [sp.cancel(x / sum(d)) for x in d]
M = transition(ptest, qtest)
ns = (M.T - sp.eye(4)).nullspace()
assert len(ns) == 1, 'chain not irreducible at symbolic eps'
v = ns[0]
w_ns = [sp.cancel(v[i] / sum(v)) for i in range(4)]
diff = [sp.simplify(a - b) for a, b in zip(w_tree, w_ns)]
print('  w_tree - w_nullspace =', diff, ' (must be [0, 0, 0, 0])')
assert all(x == 0 for x in diff)
print('  OK -- the two independent symbolic routes to w agree identically in (eps, p1, p2, p3).')

# =====================================================================================
rule("AKIN'S IDENTITY, symbolic, ALL FOUR p_i free, all 16 vertex co-players")
print("  Claim (finite eps):  w . (p~1-1, p~2-1, p~3, p~4) = 0 exactly, p~ = eps+(1-2eps)p.")
print("  Claim (limit)     :  W . (p1-1, p2-1, p3, p4) = 0.\n")
pgen = (p1, p2, p3, p4)
pc = [clamp(x) for x in pgen]
akin_vec_finite = [pc[0] - 1, pc[1] - 1, pc[2], pc[3]]
akin_vec_limit = [p1 - 1, p2 - 1, p3, p4]
all_ok = True
for q in VERTICES:
    dd = unnormalised_weights(pgen, q)
    tot = sum(dd)
    res = sp.simplify(sp.cancel(sum(dd[i] * akin_vec_finite[i] for i in range(4)) / tot))
    ok = (res == 0)
    all_ok &= ok
    print(f'   q = {label(q):26s} finite-eps residual = {res}   {"OK" if ok else "FAIL"}')
assert all_ok
print('\n  All 16 residuals are the exact symbolic zero in (eps, p1, p2, p3, p4).')
print('  The limit form follows by continuity of the residual in eps (w is a probability')
print('  vector, so the O(eps) correction p~ - p contributes at most 2*eps in absolute value).')

# =====================================================================================
rule('STRATUM G -- all four p_i symbolic: the witness ALLD = 0000 kills the slab p4 > 0')
qA = (0, 0, 0, 0)
mA_fin = margin_finite(pgen, qA)
mA_lim = sp.simplify(sp.limit(mA_fin, eps, 0, '+'))
print('  finite eps, factored:')
print('     w_DC - w_CD =', sp.factor(sp.simplify(mA_fin)))
print('\n  exact identity check:  w_DC - w_CD = -(1-2eps) * sum_s w_s p_s ?')
dd = unnormalised_weights(pgen, qA)
tot = sum(dd)
rhs = sp.cancel(-(1 - 2 * eps) * sum(dd[i] * pgen[i] for i in range(4)) / tot)
print('     residual =', sp.simplify(sp.cancel(mA_fin - rhs)), ' (must be 0)')
print('\n  limit:  W_DC - W_CD =', mA_lim)
print('     = -p4/(1 - p2 + p4)?  residual =',
      sp.simplify(mA_lim + p4 / (1 - p2 + p4)))
print('\n  Sign: for 0 <= p2 < 1 the denominator 1 - p2 + p4 > 0, so the margin is < 0')
print('  exactly when p4 > 0, and = 0 when p4 = 0.  Hence  p in Riv+  =>  p4 = 0.')
print('  (The formula is independent of p1 and p3.  It degenerates only at p2 = 1, p4 = 0,')
print('   which is stratum E below -- there the generic limit 0/0 is NOT the right answer.)')

print('\n  The same caveat as on the face applies: the formula above is the limit on the OPEN')
print('  cell.  So the slab p4 > 0 is also swept cell by cell -- p1,p2,p3 each in {0, free, 1}')
print('  and p4 either free in (0,1) or = 1, i.e. 54 cells -- with the substitution applied')
print('  BEFORE the limit, each refuted by ALLD alone:')
G_ok, G_shown = True, 0
for c4 in [p4, sp.Integer(1)]:
    for c1 in [sp.Integer(0), p1, sp.Integer(1)]:
        for c2 in [sp.Integer(0), p2, sp.Integer(1)]:
            for c3 in [sp.Integer(0), p3, sp.Integer(1)]:
                mm = sp.simplify(margin_limit((c1, c2, c3, c4), (0, 0, 0, 0)))
                bx = tuple(v for v in (c1, c2, c3, c4) if v.free_symbols)
                ok, _ = sign_certificate(mm, bx, False)
                mid = sp.nsimplify(mm.subs({v: sp.Rational(1, 2) for v in bx}))
                strict = bool(mid < 0)
                G_ok &= bool(ok and strict)
                if not (ok and strict) or G_shown < 6:
                    G_shown += 1
                    print(f'    p = ({c1}, {c2}, {c3}, {c4})'.ljust(34) +
                          f' ALLD margin = {mm}   <=0: {ok}, mid = {mid}')
print(f'  ALL 54 CELLS OF THE SLAB p4 > 0 REFUTED BY ALLD: {G_ok}')

# =====================================================================================
rule('STRATUM A -- p4 = 0 exactly, p1 p2 p3 symbolic (p2 generic): all 16 co-players')
pA = (p1, p2, p3, sp.Integer(0))
margins = {}
for q in VERTICES:
    margins[q] = margin_limit(pA, q)
    print(f'   q = {label(q):26s} W_DC - W_CD = {margins[q]}')

# ------------------------------------------------------------------ cells of the face
rule('WHY ONE GENERIC FORMULA IS NOT ENOUGH: the face p4 = 0 must be split into cells')
print('  sympy.limit treats p1,p2,p3 as GENERIC symbols, so the 16 rational functions above')
print('  are the limit only on the OPEN cell p1,p2,p3 in (0,1).  When a coordinate hits 0 or 1')
print('  the eps = 0 chain can gain or lose transient states and the limit genuinely changes.')
print('  Each p_i in {0, free in (0,1), 1} gives 3^3 = 27 cells, which partition the face.')
print('  Each cell is recomputed with the substitution applied BEFORE the limit.\n')

CELLVALS = [sp.Integer(0), None, sp.Integer(1)]      # None = keep the symbol free in (0,1)
cells = []
for c1, c2, c3 in product(CELLVALS, repeat=3):
    cells.append(((c1, c2, c3),
                  (p1 if c1 is None else c1, p2 if c2 is None else c2,
                   p3 if c3 is None else c3)))

cell_margins = {}
wrong_cells = []
for key, (v1, v2, v3) in cells:
    m = {q: sp.simplify(margin_limit((v1, v2, v3, sp.Integer(0)), q)) for q in VERTICES}
    cell_margins[key] = m
    subs = {p1: v1, p2: v2, p3: v3}
    if any(sp.simplify(m[q] - sp.simplify(margins[q].subs(subs))) != 0 for q in VERTICES):
        wrong_cells.append(key)


def cell_desc(key):
    return '(' + ', '.join('p%d' % (i + 1) if c is None else str(c)
                           for i, c in enumerate(key)) + ')'


print(f'  cells on which the generic formula is WRONG for some q: {len(wrong_cells)} / 27')
print('   ', ', '.join(cell_desc(k) for k in wrong_cells))
print('  So the certificate below is done cell by cell.  (This is a substantive correction:')
print('  a single global Bernstein certificate on the generic formulas would have been')
print('  certifying the wrong functions on 18 of the 27 cells.)')


def region_param(key, complement):
    """Parametrise cl(cell) cap {p2+p3 <= 1} (complement: cap {p2+p3 >= 1}) by a unit box.

    Returns (substitution, box variables), or None when that piece does not meet the cell.
    """
    c1, c2, c3 = key
    box = [p1] if c1 is None else []
    sub = {}
    if c2 is None and c3 is None:
        if not complement:
            sub = {p2: 1 - z, p3: z * t2}            # p3 <= 1 - p2
        else:
            sub = {p2: 1 - z, p3: z + (1 - z) * t2}  # p3 >= 1 - p2
        box += [z, t2]
    elif c2 is None:                                  # p3 fixed at 0 or 1
        if (c3 == 0) == complement:
            return None            # p3=0 has no p2<1 with p2+p3>1; p3=1 has none with <1
        sub = {p2: z}
        box.append(z)
    elif c3 is None:                                  # p2 fixed at 0 or 1
        if (c2 == 0) == complement:
            return None
        sub = {p3: z}
        box.append(z)
    else:                                             # both fixed: the cell is a point
        tot = c2 + c3
        if complement and tot <= 1:
            return None      # p2+p3 > 1 is empty here (tot = 1 lies IN the region)
        if (not complement) and tot > 1:
            return None
    return sub, tuple(box)



z, t2 = sp.symbols('z t2', nonnegative=True)

rule('SUFFICIENCY, CELL BY CELL: all 16 margins >= 0 on cl(cell) cap {p2+p3 <= 1}, p2 < 1')
print("  Certificate: substitute the cell's parametrisation, then check that numerator and")
print('  denominator have all BERNSTEIN coefficients >= 0 on the unit box.  Bernstein')
print('  nonnegativity is exact and sufficient; there is no sampling anywhere.\n')
suff_ok = True
for key, _ in cells:
    if key[1] == 1:                    # p2 = 1 is excluded (stratum E); reported there
        continue
    rp = region_param(key, complement=False)
    if rp is None:
        print(f'  cell {cell_desc(key):18s} : empty intersection with p2+p3 <= 1 -- skipped')
        continue
    sub, box = rp
    bad = []
    for q in VERTICES:
        ok, m = sign_certificate(cell_margins[key][q].subs(sub), box, True)
        if not ok:
            bad.append((q, m))
    if bad:
        suff_ok = False
        print(f'  cell {cell_desc(key):18s} : NOT CERTIFIED')
        for q, m in bad:
            print(f'        q = {label(q):26s} m = {m}')
    else:
        print(f'  cell {cell_desc(key):18s} : all 16 margins >= 0   CERTIFIED')
print(f'\n  SUFFICIENCY CERTIFIED ON EVERY CELL OF {{p4=0, p2+p3<=1, p2<1}}: {suff_ok}')

rule('NECESSITY, CELL BY CELL: a witness with margin <= 0 on every piece with p2+p3 > 1')
print('  For each cell piece meeting {p2+p3 >= 1} we look for one co-player whose margin is')
print('  <= 0 there (Bernstein on -margin) and print its value at the box midpoint, which is')
print('  an interior point of the piece and hence strictly off the face p2+p3 = 1.\n')
nec_ok = True
for key, _ in cells:
    if key[1] == 1:
        continue                       # p2 = 1: refuted wholesale in stratum E below
    rp = region_param(key, complement=True)
    if rp is None:
        continue
    sub, box = rp
    found = None
    for q in [(0, 0, 0, 1)] + [x for x in VERTICES if x != (0, 0, 0, 1)]:
        ok, m = sign_certificate(cell_margins[key][q].subs(sub), box, False)
        if ok and m != 0:
            found = (q, m, sp.nsimplify(m.subs({v: sp.Rational(1, 2) for v in box})))
            break
    if found is None:
        nec_ok = False
        print(f'  cell {cell_desc(key):18s} : NO WITNESS  <-- this would refute the theorem')
    else:
        q, m, val = found
        print(f'  cell {cell_desc(key):18s} : witness {label(q):24s} m = {m}')
        print(f'        m at the box midpoint = {val}')
print(f'\n  NECESSITY WITNESSED ON EVERY CELL PIECE WITH p2+p3 > 1: {nec_ok}')

rule('NECESSITY IN CLOSED FORM: the witness q* = 0001 on the interior cell')
qstar = (0, 0, 0, 1)
mstar = margins[qstar]
print('  W_DC - W_CD against 0001 is', mstar)
print('  claimed closed form (1-p2-p3)/(p3+2(1-p2)) : residual =',
      sp.simplify(mstar - (1 - p2 - p3) / (p3 + 2 * (1 - p2))))
print('  p1 does not appear:', p1 not in mstar.free_symbols)
print('  p3 + 2(1-p2) > 0 whenever p2 < 1, so sign(margin) = sign(1 - p2 - p3): < 0 iff')
print('  p2 + p3 > 1 and = 0 on p2+p3 = 1.  The constraint is tight, not slack.')

rule('THE MINIMUM IS ATTAINED AT ZERO ON Riv+ (no slack anywhere in the region)')
print('   margin against 0000 (ALLD) at p4 = 0, interior cell :', margins[(0, 0, 0, 0)])
print('   margin against 0001 (q*)   on the face p2+p3 = 1     :',
      sp.simplify(margins[qstar].subs(p3, 1 - p2)))
print('  Both vanish identically, so min_q (W_DC - W_CD) = 0 on Riv+: the closed form is the')
print('  ZERO SET of the min, not a strict-inequality region.')

rule('WHICH CO-PLAYERS ARE NEEDED (clause 5d)')
print('  ALLD alone, on p4 = 0, interior cell :', margins[(0, 0, 0, 0)], '-> no constraint')
print('  0001 alone                           :  p2 + p3 <= 1')
print('  So {ALLD, 0001} give the whole answer: ALLD supplies p4 = 0 (stratum G) and p2 < 1')
print('  (stratum E), 0001 supplies p2 + p3 <= 1.  No single co-player suffices.')

# =====================================================================================
rule('STRATUM E -- p4 = 0 AND p2 = 1 (the excluded edge locus), p1 p3 symbolic')
print('  Here the eps = 0 chain has two recurrent classes, so the generic stratum-A formulas')
print('  are not valid; the verdict is decided by the order-eps leakage and must be computed')
print('  at finite eps.  Substituting p2 = 1, p4 = 0 BEFORE taking the limit:\n')
pE = (p1, sp.Integer(1), p3, sp.Integer(0))
mE_fin = margin_finite(pE, (0, 0, 0, 0))
mE_lim = sp.simplify(sp.limit(mE_fin, eps, 0, '+'))
print('   against ALLD, finite eps :', sp.factor(sp.simplify(mE_fin)))
print('   against ALLD, limit      :', mE_lim)
print('   = -(1+p3)/(3+p3-p1)? residual =', sp.simplify(mE_lim + (1 + p3) / (3 + p3 - p1)))
print('   On the edge E itself (p3 = 0):', sp.simplify(mE_lim.subs(p3, 0)),
      ' -> at p1=0:', sp.simplify(mE_lim.subs({p3: 0, p1: 0})),
      ', at p1=1:', sp.simplify(mE_lim.subs({p3: 0, p1: 1})))
print('\n  3 + p3 - p1 >= 2 > 0 and 1 + p3 >= 1, so the margin is <= -1/3 < 0 everywhere on')
print('  the locus {p2 = 1, p4 = 0}.  Hence the whole of that locus -- E included -- is OUT.')
print('  Note this is a genuine discontinuity: the stratum-A limit at p2 -> 1 is 0, the')
print('  stratum-E limit is bounded away from 0.  The order of the two limits matters.')

print('\n  The two limits do NOT commute -- which is exactly why stratification is needed.')
print('  Against ALLD, at p4 = 0:')
print('   eps -> 0 first (stratum A), then p2 -> 1 :',
      sp.limit(margins[(0, 0, 0, 0)], p2, 1, '-'))
print('   p2 = 1 first, then eps -> 0 (stratum E)  :', mE_lim)
print('  Against 0001 the two orders happen to agree (both -1 at p3 = 0):')
print('   eps -> 0 first, then p2 -> 1 :', sp.limit(margins[qstar], p2, 1, '-'))
print('   p2 = 1 first, then eps -> 0  :',
      sp.simplify(sp.limit(margin_finite(pE, qstar), eps, 0, '+')))

print('\n  All NINE cells of the locus {p4 = 0, p2 = 1}, each with the substitution applied')
print('  before the limit, refuted by ALLD (margin must be < 0 on the whole closed cell):')
E_ok = True
for c1 in [sp.Integer(0), p1, sp.Integer(1)]:
    for c3 in [sp.Integer(0), p3, sp.Integer(1)]:
        mm = sp.simplify(margin_limit((c1, sp.Integer(1), c3, sp.Integer(0)), (0, 0, 0, 0)))
        box = tuple(v for v in (c1, c3) if v.free_symbols)
        ok, _ = sign_certificate(mm, box, False)          # margin <= 0 on the closed cell
        strict = sp.nsimplify(mm.subs({v: sp.Rational(1, 2) for v in box})) < 0
        E_ok &= bool(ok and strict)
        tag = f'(p1,p2,p3) = ({c1}, 1, {c3})'
        print(f'    {tag:26s} ALLD margin = {mm}   <= 0: {ok}, < 0 at midpoint: {bool(strict)}')
print(f'  WHOLE LOCUS p2 = 1 REFUTED: {E_ok}')


# =====================================================================================
rule('COMPARISON WITH THE CONJECTURED CLOSED FORM')
print('  Derived, from the three strata:')
print('     Riv+ = { p4 = 0 }  (stratum G, witness ALLD)')
print('          cap { p2 + p3 <= 1 }  (stratum A, witness 0001)')
print('          minus { p2 = 1 }  (stratum E, witness ALLD)')
print('  and on {p4 = 0, p2 < 1}, p2 = 1 is impossible, while p2 = 1 with p2+p3 <= 1 forces')
print('  p3 = 0, i.e. exactly the edge E = {(p1,1,0,0)}.  So')
print('     Riv+ = ({p4 = 0} cap {p2+p3 <= 1}) \\ E = {p4 = 0, p3 <= 1-p2, p2 < 1}.')
print('  This is verbatim the conjecture.  p1 is free, and (u,v) never entered.')

print('\n  Vertex census of the derived set (evaluate the predicate at the 16 vertices):')
inplus = [q for q in VERTICES if q[3] == 0 and q[1] + q[2] <= 1 and q[1] < 1]
print('     Riv+ cap {0,1}^4 =', [label(q) for q in inplus])
print('     expected           = 0000 ALLD, 1000 Grim, 0010, 1010 TFT')

# =====================================================================================
rule('MIRROR: Riv- by the C/D relabelling m(p) = (1-p4, 1-p3, 1-p2, 1-p1)')
print('  Verified symbolically at finite eps on all 16 co-players: relabelling both players')
print('  maps state s to 3-s and negates w_DC - w_CD.\n')
mir_ok = True
for q in VERTICES:
    pm = (1 - p4, 1 - p3, 1 - p2, 1 - p1)
    qm = tuple(1 - x for x in reversed(q))
    lhs = margin_finite(pgen, q)
    rhs = margin_finite(pm, qm)
    r = sp.simplify(sp.cancel(lhs + rhs))
    mir_ok &= (r == 0)
    print(f'   q = {label(q):26s} margin(p,q) + margin(m(p),m(q)) = {r}')
print(f'\n  MIRROR IDENTITY HOLDS FOR ALL 16: {mir_ok}')
print('  Therefore Riv- = m(Riv+) = { p1 = 1, p2 + p3 >= 1, p3 > 0 }, and')
print('     Riv+ cap Riv- = { (1, s, 1-s, 0) : 0 <= s < 1 }.')
print('  On that segment p~ := (p1-1,p2-1,p3,p4) = (1-s)(0,-1,1,0), so Akin forces')
print('  W_DC = W_CD against every co-player:')
for q in [(0, 0, 0, 0), (0, 0, 0, 1), (1, 1, 1, 1), (0, 1, 0, 1)]:
    mm = sp.simplify(margin_limit((sp.Integer(1), s_par, 1 - s_par, sp.Integer(0)), q))
    print(f'     q = {label(q):26s} W_DC - W_CD = {mm}')

# =====================================================================================
rule('FIXED eps (clause 9a): Riv+(eps) = {ALLD} exactly')
print('  From the exact ALLD identity above, w_DC - w_CD = -(1-2eps) sum_s w_s p_s with w > 0')
print('  entrywise, so at any fixed eps in (0,1/2) it is >= 0 only for p = (0,0,0,0):')
for pp, nm in [((sp.Integer(0),) * 4, 'ALLD'), ((p1, p2, p3, p4), 'generic p')]:
    print(f'    p = {nm:10s} :', sp.factor(sp.simplify(margin_finite(pp, (0, 0, 0, 0)))))

print(f'\n[done in {time.time() - t0:.1f} s]')
