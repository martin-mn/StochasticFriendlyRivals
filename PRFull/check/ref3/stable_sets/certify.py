"""Symbolic certificates, stratum by stratum, for the memory-one stability statements of SI 'Partners' and
'The dimension of the stable set' other than Eq. (3) (proved in check/ref3/eq3_boundary) and the cycle family.

Claims (domain = the strata listed; games as stated):
  B  {p_DD=0, p_CC<1} off (p_CD,p_DC) in {(0,1),(1,0)}; all games:
        stable  <=>  g1 = (1-p_CC) u + p_DC <= 0  and  g2 = (1-p_CD) u + p_DC (1+v) <= 0
  C  alternators (x,0,1,y), x<1, y>0, (x,y) != (0,1); games u+v >= 1:
        stable  <=>  h1 = (1+v-u) y - (1+u+v) <= 0  and  h2 = (1+v-u) x - 2v <= 0
  C1 (1,0,1,y), y>0, on u+v = 1:  stable <=> y(u-1) >= -1
  D  frozen (x,1,0,y), x<1, y>0, (x,y) != (0,1); games u+v >= 1:  stable <=> u - v = 1
  D1 (1,1,0,y), y>0, on u+v = 1:  stable <=> (u,v) = (1,0)

Method.  The eps->0 limits W(sigma,p) on each open stratum are rational functions of the free coordinates
(normlead.py: the normaliser lead has one sign on every open stratum, for all 17 chains).  Each gain
G_sigma = A u + B v + C is then written
  (=>) as  sum_i lambda_i * (claimed constraint or domain constraint)  -  mu,  lambda_i, mu >= 0
       (Farkas: no deviator gains on the claimed region), trying every basis of at most three of the
       multipliers; and
  (<=) the claimed constraints are themselves positive multiples of some deviators' gains (or, on a
       line domain, the gains restricted to the line).
Every sign is certified on the whole open stratum by substituting x = t/(1+t) and checking that all
coefficients of the cleared numerator have one sign (factor by factor if needed).
Run: python3 certify.py   (12 processes, about a minute)
"""
import itertools
import os
import sys
from multiprocessing import Pool

import sympy as sp

sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'farkas_display'))
import m1rows  # noqa: E402

u, v = sp.symbols('u v')
X = sp.symbols('x0:4', positive=True)
T = sp.symbols('t0:4', positive=True)
BITS = list(itertools.product([0, 1], repeat=4))
NAMES = [''.join(map(str, b)) for b in BITS]


def sign_of(expr, free):
    """+1, -1, 0 if expr has that sign on the whole open stratum; None if not certified."""
    expr = sp.factor(sp.together(expr))
    if expr == 0:
        return 0
    num, den = sp.fraction(expr)
    s = 1
    for part in (num, den):
        fl = sp.factor_list(part)
        s *= 1 if fl[0] > 0 else -1
        for fac, mult in fl[1]:
            if mult % 2 == 0:
                # even power: positive unless the factor vanishes inside; check it does not
                if one_sign(fac, free) is None:
                    return None
                continue
            sg = one_sign(fac, free)
            if sg is None:
                return None
            s *= sg
    return s


def one_sign(fac, free):
    if not fac.free_symbols:
        return 1 if fac > 0 else -1
    sub = {X[i]: T[i] / (1 + T[i]) for i in free}
    num, den = sp.fraction(sp.together(fac.subs(sub)))
    # den is a product of (1+t)^k, positive
    co = sp.Poly(sp.expand(num), *[T[i] for i in free]).coeffs()
    if all(c > 0 for c in co):
        return 1
    if all(c < 0 for c in co):
        return -1
    return None


def gains(p):
    V = m1rows.stationary_limit(p, p)
    out = []
    for b in BITS:
        W = m1rows.stationary_limit(tuple(sp.Integer(t) for t in b), p)
        d = [sp.together(W[k] - V[k]) for k in range(4)]
        G = sp.together(d[1] * (u - 1) + d[2] * v - d[3])
        num, den = sp.fraction(G)
        P = sp.Poly(sp.expand(num), u, v)
        A = sp.factor(P.coeff_monomial(u) / den)
        B = sp.factor(P.coeff_monomial(v) / den)
        C = sp.factor(P.coeff_monomial(1) / den)
        out.append((A, B, C))
    return V, out


def affine(expr):
    P = sp.Poly(sp.expand(expr), u, v)
    return (P.coeff_monomial(u), P.coeff_monomial(v), P.coeff_monomial(1))


def farkas(G, cons, free):
    """Find lambda_i >= 0 (i over cons) and mu >= 0 with G = sum lambda_i cons_i - mu (as affine forms),
    certified on the whole stratum, trying every basis of at most three multipliers (an LP with three
    equality constraints has a basic solution with at most three non-zeros).  Returns the multipliers or None."""
    n = len(cons)
    lam = sp.symbols('l0:%d' % n)
    mu = sp.Symbol('mu')
    allv = list(lam) + [mu]
    if all(sp.simplify(z) == 0 for z in G):
        return {str(w): 0 for w in allv}
    for k in (1, 2, 3):
        for basis in itertools.combinations(allv, k):
            sub0 = {w: 0 for w in allv if w not in basis}
            eq = [sp.together((sum(lam[i] * cons[i][j] for i in range(n)) - (mu if j == 2 else 0) - G[j]).subs(sub0))
                  for j in range(3)]
            sol = sp.solve(eq, list(basis), dict=True)
            if len(sol) != 1 or any(w not in sol[0] for w in basis):
                continue
            vals = dict(sub0)
            vals.update(sol[0])
            if any((vals[w].free_symbols & {u, v}) for w in basis):
                continue
            ok = all(sp.simplify(sum(vals[lam[i]] * cons[i][j] for i in range(n)) - (vals[mu] if j == 2 else 0)
                                 - G[j]) == 0 for j in range(3))
            if not ok:
                continue
            signs = [sign_of(vals[w], free) for w in allv]
            if all(sg is not None and sg >= 0 for sg in signs):
                return {str(w): sp.factor(vals[w]) for w in allv}
    return None


def strata_of(claim):
    out = []
    for st in itertools.product('0o1', repeat=4):
        cv = [None if c == 'o' else int(c) for c in st]
        pc, a, b, c = cv
        if claim == 'B':
            if c == 0 and pc != 1 and (a, b) not in ((0, 1), (1, 0)):
                out.append(st)
        elif claim in ('C', 'D'):
            ab = (0, 1) if claim == 'C' else (1, 0)
            if (a, b) == ab and pc != 1 and c != 0 and (pc, c) != (0, 1):
                out.append(st)
        elif claim == 'C1':
            if (pc, a, b) == (1, 0, 1) and c != 0:
                out.append(st)
        elif claim == 'D1':
            if (pc, a, b) == (1, 1, 0) and c != 0:
                out.append(st)
    return out


def work(job):
    claim, st = job
    free = [i for i, ch in enumerate(st) if ch == 'o']
    p = tuple(sp.Integer(0) if ch == '0' else sp.Integer(1) if ch == '1' else X[i] for i, ch in enumerate(st))
    pc, a, b, c = p
    V, G = gains(p)
    rep = {'claim': claim, 'stratum': ''.join(st), 'V': [str(w) for w in V]}
    problems = []
    if claim == 'B':
        g1 = affine((1 - pc) * u + b)
        g2 = affine((1 - a) * u + b * (1 + v))
        cons = [g1, g2]
        for name, Gs in zip(NAMES, G):
            if farkas(Gs, cons, free) is None:
                problems.append(('=>', name, [str(z) for z in Gs]))
        # converse: each claimed constraint is a positive multiple of some gain
        for k, g in enumerate(cons):
            found = None
            for name, Gs in zip(NAMES, G):
                # Gs = r * g with r > 0 ?
                r = None
                okk = True
                for j in range(3):
                    if g[j] == 0:
                        if sp.simplify(Gs[j]) != 0:
                            okk = False
                    else:
                        rj = sp.simplify(Gs[j] / g[j])
                        if r is None:
                            r = rj
                        elif sp.simplify(rj - r) != 0:
                            okk = False
                if okk and r is not None and sign_of(r, free) == 1:
                    found = (name, str(sp.factor(r)))
                    break
            if found is None:
                problems.append(('<=', 'g%d' % (k + 1), 'no deviator gain is a positive multiple'))
            else:
                rep['witness_g%d' % (k + 1)] = found
    elif claim == 'C':
        h1 = affine((1 + v - u) * c - (1 + u + v))
        h2 = affine((1 + v - u) * pc - 2 * v)
        dom = affine(1 - u - v)
        for name, Gs in zip(NAMES, G):
            if farkas(Gs, [h1, h2, dom], free) is None:
                problems.append(('=>', name, [str(z) for z in Gs]))
        for k, g in enumerate([h1, h2]):
            found = None
            for name, Gs in zip(NAMES, G):
                rs = [sp.simplify(Gs[j] / g[j]) for j in range(3) if g[j] != 0]
                zero_ok = all(sp.simplify(Gs[j]) == 0 for j in range(3) if g[j] == 0)
                if zero_ok and rs and all(sp.simplify(r - rs[0]) == 0 for r in rs) and sign_of(rs[0], free) == 1:
                    found = (name, str(sp.factor(rs[0])))
                    break
            if found is None:
                problems.append(('<=', 'h%d' % (k + 1), 'no witness'))
            else:
                rep['witness_h%d' % (k + 1)] = found
    elif claim == 'D':
        # on u = 1+v, v >= 0 every gain must be <= 0; off the line ALLD or ALLC gains
        s = sp.Symbol('s', nonnegative=True)
        for name, Gs in zip(NAMES, G):
            val0 = sp.simplify(Gs[0] * 1 + Gs[2])           # v = 0, u = 1
            slope = sp.simplify(Gs[0] + Gs[1])              # d/dv along u = 1+v
            s0, s1 = sign_of(val0, free), sign_of(slope, free)
            if s0 is None or s1 is None or s0 > 0 or s1 > 0:
                problems.append(('=>', name, str(val0), str(slope)))
        gd = [sp.simplify(z) for z in G[0]]       # ALLD
        gc = [sp.simplify(z) for z in G[15]]      # ALLC
        rep['ALLD'] = [str(z) for z in gd]
        rep['ALLC'] = [str(z) for z in gc]
        if not (sp.simplify(gd[0] + sp.Rational(1, 2)) == 0 and sp.simplify(gd[1] - sp.Rational(1, 2)) == 0
                and sp.simplify(gd[2] - sp.Rational(1, 2)) == 0):
            problems.append(('<=', 'ALLD gain is not (1+v-u)/2', [str(z) for z in gd]))
        if not (sp.simplify(gc[0] - sp.Rational(1, 2)) == 0 and sp.simplify(gc[1] + sp.Rational(1, 2)) == 0
                and sp.simplify(gc[2] + sp.Rational(1, 2)) == 0):
            problems.append(('<=', 'ALLC gain is not (u-v-1)/2', [str(z) for z in gc]))
    elif claim in ('C1', 'D1'):
        # restrict to v = 1-u: gain = (A - B) u + (B + C)
        lin = [(name, sp.simplify(Gs[0] - Gs[1]), sp.simplify(Gs[1] + Gs[2])) for name, Gs in zip(NAMES, G)]
        rep['line_gains'] = [(n, str(sp.factor(k)), str(sp.factor(m))) for n, k, m in lin]
        if claim == 'C1':
            # claimed: y(u-1) + 1 >= 0  <=>  -y u + (y - 1) <= 0 ; Farkas in one variable
            cons = [(-c, sp.Integer(0), c - 1)]
            for n, k, m in lin:
                if farkas((k, sp.Integer(0), m), cons, free) is None:
                    problems.append(('=>', n, str(k), str(m)))
            # converse: some gain on the line is a positive multiple of -y u + y - 1
            ok = any(sp.simplify(m * (-c) - k * (c - 1)) == 0 and sign_of(k / (-c), free) == 1 for n, k, m in lin
                     if k != 0)
            if not ok:
                problems.append(('<=', 'no witness on the line'))
        else:
            cons = [(sp.Integer(1), sp.Integer(0), sp.Integer(-1)), (sp.Integer(-1), sp.Integer(0), sp.Integer(1))]
            for n, k, m in lin:
                if farkas((k, sp.Integer(0), m), cons, free) is None:
                    problems.append(('=>', n, str(k), str(m)))
            pos = any(sp.simplify(k + m) == 0 and sign_of(k, free) == 1 for n, k, m in lin)
            neg = any(sp.simplify(k + m) == 0 and sign_of(k, free) == -1 for n, k, m in lin)
            if not (pos and neg):
                problems.append(('<=', 'no witnesses u<=1 and u>=1 on the line'))
    rep['problems'] = problems
    return rep


if __name__ == '__main__':
    jobs = [(cl, st) for cl in ('B', 'C', 'C1', 'D', 'D1') for st in strata_of(cl)]
    with Pool(int(os.environ.get('CERT_MAX_THREADS') or 12)) as pool:
        reps = pool.map(work, jobs, chunksize=1)
    nprob = 0
    for r in reps:
        tag = 'OK ' if not r['problems'] else 'FAIL'
        nprob += bool(r['problems'])
        extra = {k: r[k] for k in r if k.startswith('witness')}
        print(tag, r['claim'], r['stratum'], 'V =', r['V'], extra, r['problems'] if r['problems'] else '')
    print('jobs', len(reps), 'with problems', nprob)
