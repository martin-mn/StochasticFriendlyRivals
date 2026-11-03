"""Local validation of the stratcert kit (run on the Mac, nice 10, 3 cores):
 (1) for units already run by the kit (OUTDIR/t_0.tsv), a DIRECT call of the ORIGINAL library
     PRFull/check/cert: ce.certify_point(witness, wedge=locus) for every locus of the unit, compared field by
     field with the kit's line (certified, counts by route, min margin, exact self-play);
 (2) the fallback ladder: a unit whose first candidate is a deliberately non-competitive point of the same
     stratum must reject it (negative margin) and certify the next candidate.
usage: python3 prep/validate.py OUTDIR"""
import sys, os
from fractions import Fraction
KIT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/check')          # the ORIGINAL library
from cert import engine as en, exact as ex, certify as ce
assert os.path.dirname(ce.__file__) == '/Users/martin/Documents/CL1/PRFull/check/cert', ce.__file__
out = sys.argv[1]
rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(out, 't_0.tsv'))]
H = rows[0]; H[0] = H[0].lstrip('#'); rows = [dict(zip(H, r)) for r in rows[1:]]
nbad = 0; n = 0
for r in rows:
    if r['certified'] != '1': continue
    p = [Fraction(x) for x in r['witness'].split(',')]
    for loc in r['loci'].split(','):
        d = ce.certify_point(p, wedge=loc, nthreads=3); n += 1
        n_pos_exact = d['n_unresolved_by_interval'] - d['n_exact_zero']
        mine = dict(certified='1' if d['certified'] else '0', n_pos_interval=str(d['n_pos'] - n_pos_exact),
                    n_pos_exact=str(n_pos_exact), n_struct_zero=str(d['n_struct_zero']), n_exact_zero=str(d['n_exact_zero']),
                    n_neg=str(d['n_neg']), n_unres_interval=str(d['n_unresolved_by_interval']),
                    min_margin=str(Fraction(d['min_margin'])), selfplay=','.join(str(Fraction(w)) for w in d['selfplay']))
        diff = {k: (r[k], v) for k, v in mine.items() if r[k] != v}
        ok = (not diff) and d['efficiency_class'] == r['class'] == en.efficiency_class(r['pattern'])
        nbad += not ok
        print('uid %5s %s %-4s dim %s: direct certify_point(wedge=%s) %s %s' % (r['uid'], r['pattern'], loc, r['dim'], loc,
              'AGREES' if ok else 'DIFFERS', diff if diff else ''))
print('(1) %d direct calls, %d disagreements' % (n, nbad))

# (2) fallback ladder, with the kit's own module
sys.path.insert(0, KIT)
import importlib; sc = importlib.import_module('stratcert')
U = sc.load_units(); P = sc.load_points()
for uid in (7646, 8229):          # W dim 4, E dim 4
    u = U[uid]; pts = P[uid]
    # a non-competitive point of the same stratum: float scan of dyadic draws until one has a margin < -1e-6
    import random, numpy as np
    rng = random.Random(7); bad = None
    for _ in range(2000):
        q = [Fraction(int(c)) if c != '*' else Fraction(rng.randint(1, 1023), 1024) for c in u['pattern']]
        S = en.f_scan(np.array([float(x) for x in q]), None, nthreads=3)
        if (u['sign']*(S['W'][:, 2] - S['W'][:, 1])).min() < -1e-6: bad = q; break
    assert bad is not None
    res = sc.certify_unit((u, [('BAD', bad)] + pts, [], 0))
    print('(2) uid %d: bad point first -> certified %s via %s (tried %s)' % (uid, res['certified'], res['route'], res['n_tried']))
    assert res['certified'] == 1 and res['route'] == 'census:1' and res['n_tried'] == 2
    res = sc.certify_unit((u, [('BAD', bad)], [], 0))
    print('    only the bad point, no random search -> certified %s, reason: %s' % (res['certified'], res['reason']))
    assert res['certified'] == 0
    res = sc.certify_unit((u, [('BAD', bad)], [], 4096))
    print('    only the bad point, random search on -> certified %s via %s, %s s' % (res['certified'], res['route'], res['seconds']))
