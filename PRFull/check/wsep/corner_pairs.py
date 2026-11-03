"""The distinct (binary corner, co-player) pairs used by all certificates in out/*.tsv and
out/extra/*.tsv, written to out/corner_pairs.npy as int32 pairs (corner genome, q)."""
import os, glob, numpy as np, collections, time
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'out')
t0 = time.time(); per_q = collections.defaultdict(list); nrows = 0; ncorn = 0
SUB = {}
def subsets(fr):
    key = tuple(fr)
    if key not in SUB:
        w = np.array([1 << j for j in fr], dtype=np.int64); m = np.arange(1 << len(fr))
        SUB[key] = ((m[:, None] >> np.arange(len(fr))) & 1) @ w if fr else np.zeros(1, dtype=np.int64)
    return SUB[key]
for p in sorted(glob.glob(os.path.join(OUT, '*.tsv')) + glob.glob(os.path.join(OUT, 'extra', '*.tsv'))):
    if p.endswith('exact.tsv'): continue
    for l in open(p):
        a = l.split('\t'); pat = a[1]; q = int(a[4])
        base = sum(1 << j for j, ch in enumerate(pat) if ch == '1')
        fr = [j for j, ch in enumerate(pat) if ch == '*']
        per_q[q].append(base + subsets(fr)); nrows += 1; ncorn += 1 << len(fr)
pairs = []
for q, L in per_q.items():
    u = np.unique(np.concatenate(L)); pairs.append(np.stack([u, np.full(len(u), q)], 1))
P = np.concatenate(pairs).astype(np.int32)
np.save(os.path.join(OUT, 'corner_pairs.npy'), P)
print('rows', nrows, 'corner evaluations', ncorn, 'distinct co-players', len(per_q), 'distinct (corner, q) pairs', len(P),
      'seconds', round(time.time() - t0, 1))
