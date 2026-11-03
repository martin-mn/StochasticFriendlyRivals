#!/usr/bin/env python3
"""Independent check of the TYPESET box: parse every gain in si_box3_strata.tex back from LaTeX, and compare it,
and the open-facet display for every deviator the box does not list, with the exact limits of
PartnersRivalsM1f/check/exact_partners.py (a separate code base: Fractions, tree weights) at random rational
points of each stratum and random games.  Also checks L > 0 at those points against the self-play limit (1,0,0,0).

    python3 -I check_box_vs_exact_partners.py si_box3_strata.tex
    (or si_box3_strata_body.tex, the box without its title and legend, in release 1.0.0)
"""
import random
import re
import os
import sys
from fractions import Fraction as F

import sympy as sp
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application)

# exact_partners.py of the memory-one project: PartnersRivalsM1f/check, found from PRFull/check/ref3/eq3_boundary/<dir>/
_M1_CAND = [os.path.join(os.path.dirname(os.path.abspath(__file__)), *(['..'] * 5), 'PartnersRivalsM1f', 'check'),
            '/Users/martin/Documents/CL1/PartnersRivalsM1f/check']
sys.path.insert(0, next(d for d in _M1_CAND if os.path.isfile(os.path.join(d, 'exact_partners.py'))))
sys.dont_write_bytecode = True
import exact_partners as ep                                                   # noqa: E402

a, b, c, u, v, Is, Js = sp.symbols('a b c u v I K')
LOC = {'a': a, 'b': b, 'c': c, 'u': u, 'v': v, 'I': Is, 'K': Js}
TR = standard_transformations + (implicit_multiplication_application,)
NAMES = {r'\ALLD{}': '0000', r'\ALLC{}': '1111', r'\TFT{}': '1010', r'\WSLS{}': '1001', r'\Grim{}': '1000'}
DISP = {
    '1001': (0, 0, 1), '1010': (0, 0, 1), '1011': (0, 0, 1), '1101': (0, 0, 1),
    '1110': (0, 0, 1), '1111': (0, 0, 1),
    '0000': (1, 0, 1 - a + c), '0100': (1, 0, 1 - a + c),
    '1000': (2, 0, 2 - 2 * a + 3 * c), '1100': (1, 0, 2 - 2 * a + 2 * c),
    '0110': (0, 1, 1 - a + 2 * b), '0111': (0, 1, 1 - a + 2 * b),
    '0001': (1, 1 - c, 2 - 2 * a + b + c - b * c),
    '0010': (1 - b, c, (1 - a) * (1 - b) + 2 * c),
    '0011': (1 - b, 1, 2 * (1 - a + c) + 2 * b * (a - c)),
    '0101': (b, 1 - c, (1 - a) * (1 - c) + b * (3 - 2 * a)),
}


def parse_formula(t):
    # II becomes K, not J: '2J' would be a Python complex literal
    t = t.replace(r'I\!I', 'K').replace('[', '(').replace(']', ')')
    t = re.sub(r'\\tfrac\{(\d+)\}\{(\d+)\}', r'(\1/\2)', t)
    t = re.sub(r'\^(\d+)', r'**\1', t)
    return parse_expr(t, local_dict=LOC, transformations=TR)


def parse_L(cell):
    parts = re.findall(r'\\mbox\{\$(.*?)\$\}', cell)
    expr = sp.Integer(1)
    for k, ptxt in enumerate(parts):
        if k == 0:
            ptxt = ptxt.split('=', 1)[1]
        expr *= parse_formula(ptxt)
    return expr


def main():
    tex = open(sys.argv[1]).read()
    rng = random.Random(20261024)
    rows = re.findall(r'^\$\(1,([01abc]),([01abc]),([01abc])\)\$\\newline (.*?) & (.*?) & (.*?)\\\\$', tex, re.M)
    assert len(rows) == 13, len(rows)
    n_listed = n_disp = n_L = 0
    for ca, cb, cc, Lcell, gcell, _ in rows:
        sub = {}
        for sym, cd in ((a, ca), (b, cb), (c, cc)):
            if cd in '01':
                sub[sym] = sp.Integer(int(cd))
        listed = {}
        if gcell != 'none':
            for ent in gcell.split('; '):
                name, form = ent.split(':~', 1)
                name = name.replace(r'$^\dagger$', '')
                key = NAMES.get(name, name.strip('$'))
                listed[key] = parse_formula(form.strip('$'))
        Lexpr = parse_L(Lcell)
        Iexpr = (c * v - (1 - a)).subs(sub)
        Jexpr = ((1 - a) * (u - 1) + b * v).subs(sub)
        for trial in range(6):
            pt = {}
            for sym in (a, b, c):
                pt[sym] = sub[sym] if sym in sub else sp.Rational(rng.randint(1, 96), 97)
            uu, vv = sp.Rational(rng.randint(-300, 300), 97), sp.Rational(rng.randint(-300, 300), 97)
            p = tuple(F(int(sp.numer(x)), int(sp.denom(x))) for x in (1, pt[a], pt[b], pt[c]))
            V, _, ok = ep.limitW(p, p, 1)
            assert ok and V == (1, 0, 0, 0)
            assert Lexpr.subs(pt) > 0
            n_L += 1
            for q in ep.QS:
                s = ''.join(map(str, q))
                W, _, ok = ep.limitW(tuple(F(x) for x in q), p, 1)
                assert ok
                gain = W[0] + W[1] * F(int(sp.numer(uu)), int(sp.denom(uu))) \
                    + W[2] * (1 + F(int(sp.numer(vv)), int(sp.denom(vv)))) - 1
                gain = sp.Rational(gain.numerator, gain.denominator)
                if s in listed:
                    val = listed[s].subs({Is: Iexpr, Js: Jexpr}).subs(pt).subs({u: uu, v: vv})
                    n_listed += 1
                else:
                    k1, k2, den = DISP[s]
                    den = sp.sympify(den).subs(sub)
                    val = ((sp.sympify(k1) * Iexpr + sp.sympify(k2) * Jexpr) / den).subs(pt).subs({u: uu, v: vv})
                    n_disp += 1
                assert not val.has(sp.Float), val
                val = sp.simplify(val)
                assert val.is_Rational and val == gain, (ca, cb, cc, s, pt, uu, vv, val, gain, type(val))
    print('typeset box agrees with exact_partners.limitW: %d listed gains, %d display gains, %d L>0 checks; 0 mismatches'
          % (n_listed, n_disp, n_L))


if __name__ == '__main__':
    main()
