"""Fully rational re-verification of EVERY corner of EVERY certificate (no floating point).

The certificates use 42.45 M corner evaluations but only 583,817 distinct (binary corner, q) pairs
(out/corner_pairs.npy, from corner_pairs.py).  Each pair is evaluated once with ex.lead_exact
(gmpy2 rationals): its own minimum tree order E_c and the exact sign of w_DC - w_CD.  Then every
row of out/*.tsv and out/extra/*.tsv is re-checked against this table:
  every corner has E_c >= E;  every corner with E_c = E has exact margin <= 0;  at least one is < 0;
  and (#negative, #zero, #E_c>E) equal the counts the interval engine recorded.
E is the stratum's order at the centre, recorded by wsep and proved equal to the exact value for every
row by verify.py.  Together with verify.py's exact centre margins, every number the certificates use
is then checked in exact rational arithmetic.
usage: nice -n 10 python3 exact_table.py
"""
import os, sys, glob, json, time
import numpy as np
from fractions import Fraction
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'out')
sys.path.insert(0, os.path.join(HERE, '..'))
from cert import exact as ex

def ev(chunk):
    out = []
    for c, q in chunk:
        L = ex.lead_exact(ex.g2p(int(c)), ex.g2p(int(q)))
        m = L['W'][2] - L['W'][1]
        out.append((L['E'], (m > 0) - (m < 0)))
    return out

def main():
    t0 = time.time()
    P = np.load(os.path.join(OUT, 'corner_pairs.npy')).astype(np.int64)
    chunks = [P[i:i + 2000] for i in range(0, len(P), 2000)]
    with Pool(3) as pool:
        R = [x for ch in pool.map(ev, chunks) for x in ch]
    Ec = np.array([r[0] for r in R], dtype=np.int64); sg = np.array([r[1] for r in R], dtype=np.int8)
    key = P[:, 0]*65536 + P[:, 1]; order = np.argsort(key); key = key[order]; Ec = Ec[order]; sg = sg[order]
    np.savez_compressed(os.path.join(OUT, 'corner_table_exact.npz'), key=key, Ec=Ec, sign=sg)
    t1 = time.time()
    print('exact table: %d pairs, %.1f s (%.3f ms per pair per core)' % (len(P), t1 - t0, 3e3*(t1 - t0)/len(P)), flush=True)
    SUB = {}
    def subsets(fr):
        k = tuple(fr)
        if k not in SUB:
            w = np.array([1 << j for j in fr], dtype=np.int64); m = np.arange(1 << len(fr))
            SUB[k] = ((m[:, None] >> np.arange(len(fr))) & 1) @ w if fr else np.zeros(1, dtype=np.int64)
        return SUB[k]
    rep = {}
    for p in sorted(glob.glob(os.path.join(OUT, '*.tsv')) + glob.glob(os.path.join(OUT, 'extra', '*.tsv'))):
        if p.endswith('exact.tsv'): continue
        name = os.path.relpath(p, OUT)[:-4]; n = bad = 0; examples = []
        for l in open(p):
            a = l.split('\t'); pat = a[1]; E = int(a[3]); q = int(a[4])
            nneg, nz, nhigh = int(a[6]), int(a[7]) + int(a[8]), int(a[9])
            base = sum(1 << j for j, ch in enumerate(pat) if ch == '1')
            cg = base + subsets([j for j, ch in enumerate(pat) if ch == '*'])
            kk = cg*65536 + q; idx = np.searchsorted(key, kk)
            assert (key[idx] == kk).all()
            e = Ec[idx]; s = sg[idx]; same = e == E
            ok = (e >= E).all() and (s[same] <= 0).all() and (s[same] < 0).any() and \
                 (int((s[same] < 0).sum()), int((s[same] == 0).sum()), int((e > E).sum())) == (nneg, nz, nhigh)
            n += 1
            if not ok:
                bad += 1
                if len(examples) < 5: examples.append(pat)
        rep[name] = dict(rows=n, rows_failing=bad, examples=examples)
        print(name, rep[name], flush=True)
    rep['_pairs'] = len(P); rep['_seconds'] = round(time.time() - t0, 1)
    json.dump(rep, open(os.path.join(OUT, 'exact_table.json'), 'w'), indent=1)

if __name__ == '__main__':
    main()
