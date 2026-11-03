"""Prover 2 -- independent search for parameter values where the formulas could break.

Takes every numerator/denominator polynomial recorded by p2_proof.py (denominators of mu_i,
of absorption probabilities, of Q entries, of w, of W, the determinants of the stationary and
transient linear systems, and numerators/denominators of l1, l2, mu) and searches the OPEN
stratum for real zeros, without using the x = t/(1+t) certificates:
  * 1 free parameter: exact Sturm count of real roots in the open interval (0,1);
  * 2 free parameters: exact Sturm counts on 2 x 101 axis-parallel lines (including the lines
    at 10^-6 and 1-10^-6), plus a check that the polynomial has no zero on the critical set
    {f = 0, df/dx = 0} inside the box (resultant in x, Sturm in y), which would be needed for
    a zero curve to avoid every line;
  * 3 free parameters (open facet): Sturm counts on 3 x 21^2 axis-parallel lines.
Run: python3 p2_rootsearch.py   (~1 min)
"""
import os, pickle, sys, time
from fractions import Fraction as F
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from p2lib import A, B, C, stratum_values

X = sp.Symbol('x')


def open_roots_1d(f, x):
    """Number of distinct real roots of the univariate polynomial f in the OPEN interval (0,1)."""
    P = sp.Poly(f, x)
    if P.is_zero:
        return float('inf')
    if P.degree() <= 0:
        return 0
    P = sp.Poly(sp.sqf_part(P.as_expr()), x)
    n = int(P.count_roots(0, 1))
    n -= int(P.eval(0) == 0) + int(P.eval(1) == 0)
    return n


def lines(N):
    vals = [sp.Rational(k, N) for k in range(1, N)] + [sp.Rational(1, 10 ** 6), 1 - sp.Rational(1, 10 ** 6)]
    return vals


def search(code, f, free):
    f = sp.expand(f)
    if not free or not f.free_symbols:
        return 0 if f != 0 else 1, 'constant %s' % f
    if len(free) == 1:
        n = open_roots_1d(f, free[0])
        return n, 'Sturm on (0,1): %d roots' % n
    if len(free) == 2:
        x, y = free
        bad = 0
        for val in lines(100):
            bad += open_roots_1d(f.subs(x, val), y) > 0 and 1 or 0
            bad += open_roots_1d(f.subs(y, val), x) > 0 and 1 or 0
        # critical set: zeros of f with df/dx = 0; project to y by the resultant
        g = sp.diff(f, x)
        crit = 0
        if sp.Poly(f, x).degree() >= 1 and g != 0:
            r = sp.resultant(f, g, x)
            if r != 0 and sp.Poly(r, y).degree() >= 1:
                for yr in sp.Poly(r, y).real_roots():
                    if 0 < yr < 1:
                        # does f(., yr) vanish in (0,1) at a point where df/dx = 0 ?
                        fy = sp.Poly(sp.expand(f.subs(y, yr)), x)
                        for xr in fy.real_roots():
                            if 0 < xr < 1:
                                crit += 1
        return bad + crit, 'lines with a root: %d; critical zeros in box: %d' % (bad, crit)
    if len(free) == 3:
        bad = 0
        vals = [sp.Rational(k, 20) for k in range(1, 20)] + [sp.Rational(1, 10 ** 6), 1 - sp.Rational(1, 10 ** 6)]
        for i in range(3):
            others = [free[j] for j in range(3) if j != i]
            for v1 in vals:
                for v2 in vals:
                    bad += open_roots_1d(f.subs({others[0]: v1, others[1]: v2}), free[i]) > 0 and 1 or 0
        return bad, 'lines with a root: %d' % bad
    raise ValueError


if __name__ == '__main__':
    t0 = time.time()
    data = pickle.load(open(os.path.join(HERE, 'p2_polys.pkl'), 'rb'))
    seen = {}
    for code, sg, role, srep in data:
        f = sp.sympify(sp.parse_expr(srep) if False else eval(srep, sp.__dict__.copy()))
        f = f.subs({sp.Symbol('a', positive=True): A, sp.Symbol('b', positive=True): B,
                    sp.Symbol('c', positive=True): C})
        f = sp.expand(f)
        if f.is_number:
            continue
        # search irreducible factors (a zero of a product is a zero of a factor)
        _, free = stratum_values(code)
        for fac, _m in sp.factor_list(f, *free)[1]:
            key = (code, sp.srepr(sp.expand(fac)))
            seen.setdefault(key, []).append((sg, role))
    print('distinct irreducible factors to search:', len(seen))
    total_bad = 0
    bycode = {}
    for (code, srep), uses in sorted(seen.items()):
        _, free = stratum_values(code)
        f = eval(srep, sp.__dict__.copy())
        n, how = search(code, f, free)
        bycode.setdefault(code, []).append((f, n, how))
        if n:
            total_bad += 1
            print('ZERO FOUND in open stratum', code, f, how, uses[:4])
    for code, lst in bycode.items():
        print('%s: %d distinct factors searched, zeros in open stratum: %d   e.g. %s'
              % (code, len(lst), sum(1 for _, n, _ in lst if n), [str(f) for f, _, _ in lst[:6]]))
    print('TOTAL factors with a zero in the open stratum:', total_bad, ' time %.1fs' % (time.time() - t0))
