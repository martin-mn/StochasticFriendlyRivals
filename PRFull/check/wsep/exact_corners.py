"""Exact (rational, no floating point) re-check of the corner tables of a random sample of the
certificates in out/<cls>.tsv: every corner c of the stratum cube is evaluated with ex.lead_exact
against the certifying q; E_c >= E must hold; corners with E_c = E must have exact margin <= 0 and at
least one must be < 0; and the counts (negative, zero, E_c > E) must equal those the interval
engine recorded.  usage: nice -n 10 python3 exact_corners.py [per-class sample] [max d]
"""
import sys, os, random, itertools, json, time
from fractions import Fraction
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from cert import exact as ex

def job(row):
    pat, E, q, nneg, nz, nhigh = row
    fr = [j for j, ch in enumerate(pat) if ch == '*']
    base = [Fraction(1, 2) if ch == '*' else Fraction(int(ch)) for ch in pat]
    Q = ex.g2p(q)
    if ex.lead_exact(base, Q)['E'] != E: return pat, 'E mismatch'
    a = b = c = 0
    for bits in itertools.product((0, 1), repeat=len(fr)):
        x = list(base)
        for j, v in zip(fr, bits): x[j] = Fraction(v)
        L = ex.lead_exact(x, Q)
        if L['E'] < E: return pat, 'E_c < E'
        if L['E'] > E: c += 1; continue
        m = L['W'][2] - L['W'][1]
        if m > 0: return pat, 'positive corner'
        if m < 0: a += 1
        else: b += 1
    if a < 1: return pat, 'no negative corner'
    if (a, b, c) != (nneg, nz, nhigh): return pat, 'counts differ %s %s' % ((a, b, c), (nneg, nz, nhigh))
    return pat, 'ok'

def main():
    ns = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    rng = random.Random(77); t0 = time.time(); rep = {}
    for name in ['c6', 'c8', 'c68', 'c9a', 'c9b']:
        R = []
        for l in open(os.path.join(HERE, 'out', name + '.tsv')):
            a = l.split('\t')
            if int(a[2]) <= dmax: R.append((a[1], int(a[3]), int(a[4]), int(a[6]), int(a[7]) + int(a[8]), int(a[9])))
        S = rng.sample(R, min(ns, len(R)))
        with Pool(3) as P:
            res = P.map(job, S, chunksize=8)
        bad = [r for r in res if r[1] != 'ok']
        ncorn = sum(1 << r[0].count('*') for r in S)
        rep[name] = dict(sampled=len(S), corners=ncorn, failures=len(bad), examples=bad[:3], max_d=max(r[0].count('*') for r in S))
        print(name, rep[name], round(time.time() - t0, 1), 's', flush=True)
    json.dump(rep, open(os.path.join(HERE, 'out', 'exact_corners.json'), 'w'), indent=1)

if __name__ == '__main__':
    main()
