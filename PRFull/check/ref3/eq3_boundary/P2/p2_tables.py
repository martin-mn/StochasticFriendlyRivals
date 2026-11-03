"""Compact per-stratum table from p2_summary.json: closed-class counts, multipliers of ALLD and
0111, and the distinct irreducible factors that occur in all certified quantities."""
import json, os, sys
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from p2lib import A, B, C, stratum_values
LOC = {'a': A, 'b': B, 'c': C}
S = json.load(open(os.path.join(HERE, 'p2_summary.json')))
for code, d in S.items():
    _, free = stratum_values(code)
    ks = [r['k'] for r in d['rows'].values()]
    facs = set()
    for r in d['rows'].values():
        for key in ('l1', 'l2', 'mu'):
            e = sp.factor(sp.sympify(r[key], locals=LOC))
            n, dd = sp.fraction(e)
            for part in (n, dd):
                if free and not part.is_number:
                    for f, m in sp.factor_list(part, *free)[1]:
                        facs.add(str(f))
        for w in r['W']:
            e = sp.factor(sp.sympify(w, locals=LOC))
            for part in sp.fraction(e):
                if free and not part.is_number:
                    for f, m in sp.factor_list(part, *free)[1]:
                        facs.add(str(f))
    g = lambda s: (d['rows'][s]['l1'], d['rows'][s]['l2'], d['rows'][s]['mu'])
    print('(1,%s,%s,%s): #sigma with 1/2/3 closed classes of M0: %d/%d/%d; mu==0 for all: %s'
          % (code[0], code[1], code[2], ks.count(1), ks.count(2), ks.count(3), not d['mu_nonzero']))
    print('     ALLD (l1,l2,mu)=%s ; 0111 (l1,l2,mu)=%s' % (g('0000'), g('0111')))
    print('     irreducible factors in W, l1, l2, mu: %s' % sorted(facs, key=len))
