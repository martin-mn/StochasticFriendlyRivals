#!/usr/bin/env python3
"""Prover P1: exact proof of the partner conditions I<=0, II<=0 on the 13 boundary strata
of R = {p = (1,a,b,c): c > 0} minus the segments (a,b) = (0,1), (1,0), plus a re-derivation
of the open-facet display of the SI 'Partners' paragraph.

    python3 prove_strata.py            (writes prove_strata.log and prove_strata.json here)

Self-contained: it imports only sympy.  Nothing is floating point; no tolerance anywhere.

SETTING.  States CC=0, CD=1, DC=2, DD=3 from the deviator's (focal) view.  The deviator
plays sigma[s]; the resident answers p[M[s]] with M swapping CD and DC.  Errors clamp every
cooperation probability x to eps + (1-2 eps) x.  For 0 < eps < 1/2 every transition
probability is positive, so the chain is irreducible, and by the Markov chain tree theorem
(equivalently: adj(I-P) has rank one with rows proportional to the stationary vector) the
stationary distribution is w_X = N_X(eps)/Ntot(eps) with N_X the principal 3x3 minor of
I - P with row/column X deleted.  The code checks sum_X N_X (I-P)_{XY} == 0 identically.

LIMIT, EXACTLY AND AT EVERY POINT OF A STRATUM.  N_X and Ntot are polynomials in eps whose
coefficients are polynomials in the free stratum parameters.  Let k0 be the lowest eps-order
whose coefficient L(params) of Ntot is not the zero polynomial.  (i) we check that the
eps^k coefficients of every N_X vanish IDENTICALLY for k < k0 (then also those of Ntot);
(ii) we certify that L(params) > 0 at EVERY point of the open stratum (sign certificate
below).  Then at every point  W_X = lim w_X = [eps^k0]N_X / L, a rational function whose
denominator never vanishes on the stratum.

SIGN CERTIFICATES.  A rational function is factored over Q; every irreducible factor f is
certified sign-definite on the open box (0,1)^k of the free parameters by
  'manifest'  : f is one of x, 1-x (x a free parameter) or a non-zero rational constant;
  'subst'     : x -> t/(1+t) for every free x, multiplied by (1+t)^{deg_x f}; the resulting
                polynomial in t > 0 is non-zero and all its coefficients have one strict sign
                (these are binomial-scaled Bernstein coefficients of f on [0,1]^k), so f has
                that sign at every interior point.
A claim 'q >= 0' is proved either because q == 0 identically or because q > 0 is certified.

GAIN.  g_sigma = (u-1) W_CD + v W_DC - W_DD  (needs W(p,p) = (1,0,0,0), which is proved).
Decomposition g = l1 I + l2 II - mu, I = c v - (1-a), II = (1-a)(u-1) + b v:
  a < 1 : l2 = W_CD/(1-a),  l1 = (W_DC - b l2)/c,  mu = W_DD - (1-a) l1   (unique);
  a = 1 : I = c v, II = b v; need W_CD == 0, and then l1 = W_DC/c, l2 = 0, mu = W_DD.
Sufficiency: l1, l2, mu >= 0 for all 16 sigma.
Necessity : a sigma with g = lam*I, lam > 0 (l2 == mu == 0), and one with g = lam'*II,
            lam' > 0 (l1 == mu == 0).  For a = 1, I > 0 or II > 0 iff v > 0 (b, c > 0),
            and a sigma with W_CD == W_DD == 0 < W_DC gives g = W_DC v > 0.
"""
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from sympy import QQ
from sympy.polys.rings import ring

HERE = os.path.dirname(os.path.abspath(__file__))
R, E, A, B, C = ring('e,a,b,c', QQ)
a_s, b_s, c_s, e_s = sp.symbols('a b c e')
u, v = sp.symbols('u v')
SYM = {'a': a_s, 'b': b_s, 'c': c_s}
RG = {'a': A, 'b': B, 'c': C}
M = (0, 2, 1, 3)
NAMES = {'0000': 'ALLD', '1111': 'ALLC', '1010': 'TFT', '1001': 'WSLS', '1000': 'Grim'}
SIGMAS = [''.join(map(str, bits)) for bits in itertools.product([0, 1], repeat=4)]

LOG = []


def log(*args):
    s = ' '.join(str(x) for x in args)
    LOG.append(s)
    print(s)
    sys.stdout.flush()


# ----------------------------------------------------------------- the chain and its limit

def cofactors(sig, p):
    """sig, p: 4-tuples of ring elements.  Returns [N_CC, N_CD, N_DC, N_DD] (ring elements)."""
    L = [[None] * 4 for _ in range(4)]
    for s in range(4):
        x = E + (1 - 2 * E) * sig[s]
        y = E + (1 - 2 * E) * p[M[s]]
        row = (x * y, x * (1 - y), (1 - x) * y, (1 - x) * (1 - y))
        for j in range(4):
            L[s][j] = (R(1) if s == j else R(0)) - row[j]
    # rows of I-P sum to zero (stochasticity), identically
    for s in range(4):
        assert sum(L[s]) == 0
    N = []
    for k in range(4):
        idx = [i for i in range(4) if i != k]
        m = [[L[i][j] for j in idx] for i in idx]
        N.append(m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                 - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                 + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    # stationarity N (I-P) == 0, identically in eps and the parameters
    for j in range(4):
        assert sum((N[i] * L[i][j] for i in range(4)), R(0)) == 0
    return N


def by_order(poly):
    """{k: coefficient of e^k as a ring element free of e}"""
    out = {}
    for mon, cf in poly.terms():
        k = mon[0]
        out[k] = out.get(k, R(0)) + R.from_dict({(0,) + tuple(mon[1:]): cf})
    return out


def to_expr(rel):
    return rel.as_expr(e_s, a_s, b_s, c_s)


# ----------------------------------------------------------------- sign certificates

def factor_cert(f, free):
    """f: irreducible sympy polynomial factor.  Returns (sign, certificate string)."""
    fv = [x for x in free if f.has(x)]
    if not fv:
        val = sp.Rational(f)
        assert val != 0
        return (1 if val > 0 else -1), 'const'
    for x in fv:
        if sp.expand(f - x) == 0:
            return 1, 'manifest:%s' % x
        if sp.expand(f - (1 - x)) == 0:
            return 1, 'manifest:1-%s' % x
        if sp.expand(f + (1 - x)) == 0:
            return -1, 'manifest:-(1-%s)' % x
    ts = sp.symbols('t1:%d' % (len(fv) + 1), positive=True)
    g = f
    for x, t in zip(fv, ts):
        d = sp.degree(f, x)
        g = sp.together(g.subs(x, t / (1 + t))) * (1 + t) ** d
    g = sp.expand(sp.cancel(g))
    P = sp.Poly(g, *ts)
    cfs = [cf for cf in P.coeffs() if cf != 0]
    assert cfs, 'zero after substitution'
    if all(cf > 0 for cf in cfs):
        return 1, 'subst(%d coeffs>0)' % len(cfs)
    if all(cf < 0 for cf in cfs):
        return -1, 'subst(%d coeffs<0)' % len(cfs)
    return None, 'FAILED subst on %s' % f


def sign_cert(expr, free):
    """Certified strict sign of a rational function on the open box of `free`.
    Returns (sign, certs) with sign in {1,-1} or (0, ['identically 0']) or (None, ...)."""
    expr = sp.cancel(sp.together(expr))
    if expr == 0:
        return 0, ['identically 0']
    num, den = sp.fraction(expr)
    sign = 1
    certs = []
    for part, tag in ((num, 'num'), (den, 'den')):
        cst, facs = sp.factor_list(sp.expand(part), *free) if free else (part, [])
        # the factorisation is re-multiplied and compared with the input (no trust in factor_list)
        assert sp.expand(cst * sp.Mul(*[f ** k for f, k in facs]) - part) == 0
        cst = sp.Rational(cst)
        assert cst != 0
        if cst < 0:
            sign = -sign
        for f, k in facs:
            s, cert = factor_cert(f, free)
            if s is None:
                return None, certs + ['%s factor %s: %s' % (tag, f, cert)]
            if s < 0 and k % 2 == 1:
                sign = -sign
            certs.append('%s[%s]^%d:%s' % (tag, sp.sstr(f), k, cert))
    return sign, certs


def nonneg(expr, free):
    s, certs = sign_cert(expr, free)
    return (s is not None and s >= 0), s, certs


def positive(expr, free):
    s, certs = sign_cert(expr, free)
    return s == 1, certs


# ----------------------------------------------------------------- limits on a stratum

def limit_W(sig, p, free):
    """Exact eps->0 limit of the stationary distribution, valid at every point of the open
    stratum.  Returns (W list of sympy exprs, info dict)."""
    N = cofactors(sig, p)
    tot = N[0] + N[1] + N[2] + N[3]
    T = by_order(tot)
    k0 = min(k for k, cf in T.items() if cf != 0)
    lead = to_expr(T[k0])
    # (i) lower orders of every numerator vanish identically
    for X in range(4):
        for k, cf in by_order(N[X]).items():
            if k < k0:
                assert cf == 0, 'numerator %d has a non-zero eps^%d term' % (X, k)
    # (ii) leading coefficient of the normaliser sign-definite (positive) on the stratum
    ok, certs = positive(lead, free)
    nums = [to_expr(by_order(N[X]).get(k0, R(0))) for X in range(4)]
    W = [sp.factor(sp.cancel(n / lead)) for n in nums]
    assert sp.simplify(sum(W) - 1) == 0
    return W, {'k0': k0, 'lead': sp.factor(lead), 'lead_ok': ok, 'lead_certs': certs}


def point_limit(sig, p_vals):
    """Independent cross-check: exact limit at one rational point (univariate in eps)."""
    pr = tuple(R(sp.Rational(x)) if not isinstance(x, type(A)) else x for x in p_vals)
    N = cofactors(sig, pr)
    tot = N[0] + N[1] + N[2] + N[3]
    T = by_order(tot)
    k0 = min(k for k, cf in T.items() if cf != 0)
    lead = to_expr(T[k0])
    return [sp.Rational(to_expr(by_order(N[X]).get(k0, R(0)))) / sp.Rational(lead) for X in range(4)]


# ----------------------------------------------------------------- strata

def strata():
    out = []
    for ca in '01*':
        for cb in '01*':
            for cc in '*1':
                if (ca, cb, cc) == ('*', '*', '*'):
                    continue
                if (ca, cb) in (('0', '1'), ('1', '0')):
                    continue
                out.append((ca, cb, cc))
    return out


def setup(code):
    ca, cb, cc = code
    vals_r, vals_s, free = [], [], []
    for name, cd in zip('abc', code):
        if cd == '*':
            vals_r.append(RG[name]); vals_s.append(SYM[name]); free.append(SYM[name])
        else:
            vals_r.append(R(int(cd))); vals_s.append(sp.Integer(int(cd)))
    return vals_r, vals_s, free


def sig_tuple(s):
    return tuple(R(int(ch)) for ch in s)


def label(code):
    return '(1,%s,%s,%s)' % code


def prove_stratum(code, rng):
    vals_r, (a, b, c), free = setup(code)
    p = (R(1),) + tuple(vals_r)
    I = c * v - (1 - a)
    II = (1 - a) * (u - 1) + b * v
    res = {'stratum': label(code), 'free': [str(x) for x in free]}
    log('=' * 100)
    log('STRATUM %s   free parameters %s in (0,1)' % (label(code), [str(x) for x in free]))
    # self-play
    Wpp, inf = limit_W(p, p, free)
    sp_ok = [sp.simplify(w) for w in Wpp] == [1, 0, 0, 0] and inf['lead_ok']
    log('  self-play W(p,p) = %s ; eps-order %d, normaliser lead %s, certified>0: %s %s'
        % (Wpp, inf['k0'], inf['lead'], inf['lead_ok'], inf['lead_certs']))
    res['selfplay_cc'] = bool(sp_ok)
    limit_ok = inf['lead_ok']
    suff_ok = True
    rows = {}
    for s in SIGMAS:
        W, inf = limit_W(sig_tuple(s), p, free)
        limit_ok = limit_ok and inf['lead_ok']
        g = sp.expand((u - 1) * W[1] + v * W[2] - W[3])
        row = {'W': [str(w) for w in W], 'k0': inf['k0'], 'lead': str(inf['lead']),
               'lead_ok': inf['lead_ok'], 'lead_certs': inf['lead_certs']}
        if a != 1:
            l2 = sp.factor(sp.cancel(W[1] / (1 - a)))
            l1 = sp.factor(sp.cancel((W[2] - b * l2) / c))
            mu = sp.factor(sp.cancel(W[3] - (1 - a) * l1))
            ok_cd = True
        else:
            ok_cd = sp.simplify(W[1]) == 0
            l2 = sp.Integer(0)
            l1 = sp.factor(sp.cancel(W[2] / c))
            mu = sp.factor(W[3])
        # exact identity g == l1 I + l2 II - mu
        ident = sp.simplify(sp.cancel(g - (l1 * I + l2 * II - mu))) == 0
        certs = {}
        okl = ok_cd and ident
        for nm, q in (('l1', l1), ('l2', l2), ('mu', mu)):
            ok, sgn, cs = nonneg(q, free)
            certs[nm] = cs
            okl = okl and ok
        suff_ok = suff_ok and okl
        row.update({'l1': str(l1), 'l2': str(l2), 'mu': str(mu), 'identity': ident,
                    'W_CD_zero(a=1)': ok_cd if a == 1 else None, 'certs': certs, 'ok': okl})
        rows[s] = row
        log('  %s %-5s k0=%d lead=%s [%s]  W=%s' % (s, NAMES.get(s, ''), inf['k0'], inf['lead'],
                                                   'pos:' + ';'.join(inf['lead_certs']) if inf['lead_ok'] else 'LEAD NOT CERTIFIED',
                                                   W))
        log('        g = %s*I + %s*II - %s   identity:%s  ok:%s' % (l1, l2, mu, ident, okl))
        log('        certs l1:%s | l2:%s | mu:%s' % (certs['l1'], certs['l2'], certs['mu']))
        # point cross-check at random rational interior points
        for _ in range(2):
            pt = {x: sp.Rational(rng.randint(1, 96), 97) for x in free}
            pv = [x.subs(pt) if hasattr(x, 'subs') else x for x in (sp.Integer(1), a, b, c)]
            Wp = point_limit(sig_tuple(s), pv)
            assert Wp == [sp.Rational(w.subs(pt)) for w in W], ('point mismatch', s, pt)
    res['limit_valid'] = bool(limit_ok)
    res['sufficiency'] = bool(suff_ok)
    # necessity
    witnesses = []
    if a != 1:
        wI = [s for s in SIGMAS if rows[s]['ok'] and sp.sympify(rows[s]['l2']) == 0
              and sp.sympify(rows[s]['mu']) == 0 and positive(sp.sympify(rows[s]['l1']), free)[0]]
        wII = [s for s in SIGMAS if rows[s]['ok'] and sp.sympify(rows[s]['l1']) == 0
               and sp.sympify(rows[s]['mu']) == 0 and positive(sp.sympify(rows[s]['l2']), free)[0]]
        nec = bool(wI) and bool(wII)
        witnesses = {'I': wI, 'II': wII}
        log('  NECESSITY: g = lam*I, lam>0 certified, for %s ; g = lam*II, lam>0 for %s' % (wI, wII))
        if '0000' in wI:
            log('     ALLD: g = %s * I' % rows['0000']['l1'])
        if '0111' in wII:
            log('     0111: g = %s * II' % rows['0111']['l2'])
    else:
        wv = [s for s in SIGMAS if sp.simplify(sp.sympify(rows[s]['W'][1])) == 0
              and sp.simplify(sp.sympify(rows[s]['W'][3])) == 0
              and positive(sp.sympify(rows[s]['W'][2]), free)[0]]
        nec = bool(wv)
        witnesses = {'v': wv}
        log('  NECESSITY (a=1, I>0 or II>0 iff v>0): g = W_DC * v with W_DC>0, W_CD=W_DD=0 for %s' % wv)
    res['necessity'] = nec
    res['witnesses'] = witnesses
    res['rows'] = rows
    log('  RESULT %s: selfplay_cc=%s limit_valid=%s sufficiency=%s necessity=%s'
        % (label(code), res['selfplay_cc'], res['limit_valid'], res['sufficiency'], res['necessity']))
    return res


# ----------------------------------------------------------------- open facet display

def open_facet(rng):
    a, b, c = a_s, b_s, c_s
    free = [a, b, c]
    p = (R(1), A, B, C)
    I = c * v - (1 - a)
    II = (1 - a) * (u - 1) + b * v
    disp = {
        '1001': (0, 0, 1), '1010': (0, 0, 1), '1011': (0, 0, 1), '1101': (0, 0, 1),
        '1110': (0, 0, 1), '1111': (0, 0, 1),
        '0000': (1, 0, 1 - a + c), '0100': (1, 0, 1 - a + c),
        '1000': (2, 0, 2 - 2 * a + 3 * c),
        '1100': (1, 0, 2 - 2 * a + 2 * c),
        '0110': (0, 1, 1 - a + 2 * b), '0111': (0, 1, 1 - a + 2 * b),
        '0001': (1, 1 - c, 2 - 2 * a + b + c - b * c),
        '0010': (1 - b, c, (1 - a) * (1 - b) + 2 * c),
        '0011': (1 - b, 1, 2 * (1 - a + c) + 2 * b * (a - c)),
        '0101': (b, 1 - c, (1 - a) * (1 - c) + b * (3 - 2 * a)),
    }
    log('=' * 100)
    log('OPEN FACET 0<a,b,c<1: re-derivation of the SI display')
    Wpp, inf = limit_W(p, p, free)
    selfplay = [sp.simplify(w) for w in Wpp] == [1, 0, 0, 0] and inf['lead_ok']
    log('  self-play W(p,p) = %s ; lead %s certified: %s %s' % (Wpp, inf['lead'], inf['lead_ok'], inf['lead_certs']))
    all_ok = selfplay
    rows = {}
    for s in SIGMAS:
        W, inf = limit_W(sig_tuple(s), p, free)
        g = (u - 1) * W[1] + v * W[2] - W[3]
        k1, k2, den = disp[s]
        match = sp.simplify(sp.cancel(g - (k1 * I + k2 * II) / den)) == 0
        okd, cd = positive(den, free)
        o1, s1, c1 = nonneg(k1, free)
        o2, s2, c2 = nonneg(k2, free)
        ok = match and okd and o1 and o2 and inf['lead_ok']
        all_ok = all_ok and ok
        rows[s] = {'match': match, 'den_pos': okd, 'k1_nonneg': o1, 'k2_nonneg': o2,
                   'lead_ok': inf['lead_ok']}
        log('  %s %-5s lead=%s [%s] gain = [(%s) I + (%s) II]/(%s)  match:%s den>0:%s %s  coeffs>=0:%s,%s'
            % (s, NAMES.get(s, ''), inf['lead'], ';'.join(inf['lead_certs']), k1, k2, den, match, okd,
               cd, o1, o2))
        for _ in range(2):
            pt = {x: sp.Rational(rng.randint(1, 96), 97) for x in free}
            Wp = point_limit(sig_tuple(s), [1, pt[a], pt[b], pt[c]])
            assert Wp == [sp.Rational(w.subs(pt)) for w in W]
    log('  OPEN FACET: all 16 displayed gains confirmed (mu = 0, denominators > 0, coefficients >= 0, limits valid): %s' % all_ok)
    return all_ok, rows


def main():
    t0 = time.time()
    rng = random.Random(20261024)
    sts = strata()
    log('13 boundary strata:', ', '.join(label(s) for s in sts), '  (count %d)' % len(sts))
    assert len(sts) == 13
    results = [prove_stratum(s, rng) for s in sts]
    of_ok, of_rows = open_facet(rng)
    log('=' * 100)
    log('SUMMARY')
    for r in results:
        log('  %-14s free=%-12s selfplay_cc=%s limit_valid=%s sufficiency=%s necessity=%s witnesses=%s'
            % (r['stratum'], ','.join(r['free']) or '-', r['selfplay_cc'], r['limit_valid'],
               r['sufficiency'], r['necessity'], r['witnesses']))
    overall = all(r['selfplay_cc'] and r['limit_valid'] and r['sufficiency'] and r['necessity']
                  for r in results) and of_ok
    log('  open facet display confirmed: %s' % of_ok)
    log('  OVERALL PROVED: %s   (%.1f s)' % (overall, time.time() - t0))
    with open(os.path.join(HERE, 'prove_strata.log'), 'w') as fh:
        fh.write('\n'.join(LOG) + '\n')
    with open(os.path.join(HERE, 'prove_strata.json'), 'w') as fh:
        json.dump({'strata': results, 'open_facet_ok': of_ok, 'open_facet_rows': of_rows,
                   'overall': overall}, fh, indent=1, default=str)


if __name__ == '__main__':
    main()
