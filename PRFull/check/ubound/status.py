"""status of an exclusion run: per locus, the efficient strata of each dimension not yet excluded"""
import sys, numpy as np
N3 = 3**16
_cache = {}
def base():
    if not _cache:
        cls = np.fromfile('cls.bin', dtype=np.uint8); dim = np.load('dim.npy')
        _cache.update(cls=cls, dim=dim)
    return _cache['cls'], _cache['dim']
def digit(idx, j): return (idx // 3**j) % 3
def load(prefix):
    idx = np.fromfile(prefix + '.idx', dtype=np.uint32).astype(np.int64)
    qp = np.fromfile(prefix + '.qp', dtype=np.int32); qm = np.fromfile(prefix + '.qm', dtype=np.int32)
    return idx, qp, qm
# locus: (class set, sign, bound = smallest dimension that must be excluded, extra filter)
LOCI = {
    'W':     ({1}, +1, 8, None),
    'N':     ({2}, +1, 9, None),
    'E':     ({2}, -1, 9, None),
    'ray+':  ({1, 2, 3}, +1, 9, None),
    'W68':   ({1}, +1, 6, 'p6=p8=0'),
    'S':     ({1}, -1, 16, None),
    'ray-':  ({1, 2, 3}, -1, 16, None),
}
import os, json
for _k, _b in json.loads(os.environ.get('UB_BOUNDS', '{}')).items():     # e.g. UB_BOUNDS='{"W": 0}'
    LOCI[_k] = LOCI[_k][:2] + (_b,) + LOCI[_k][3:]
def masks(idx):
    cls, dim = base(); c = cls[idx]; d = dim[idx]
    out = {}
    for k, (cs, sg, b, extra) in LOCI.items():
        m = np.isin(c, list(cs))
        if extra == 'p6=p8=0': m &= (digit(idx, 6) == 0) & (digit(idx, 8) == 0)
        out[k] = (m, sg, b)
    return out, d
def pattern(i):
    s = []
    for j in range(16): s.append('01*'[i % 3]); i //= 3
    return ''.join(s)
def report(prefix, show=True):
    idx, qp, qm = load(prefix); M, d = masks(idx); res = {}
    for k, (m, sg, b) in M.items():
        q = qp if sg > 0 else qm
        left = m & (q < 0)
        hist_all = np.bincount(d[m], minlength=17); hist_left = np.bincount(d[left], minlength=17)
        tgt = m & (d >= b); tleft = tgt & (q < 0)
        res[k] = dict(targets=int(tgt.sum()), left=int(tleft.sum()), hist_all=hist_all.tolist(), hist_left=hist_left.tolist(),
                      maxdim_left=int(d[left].max()) if left.any() else -1)
        if show:
            print('%-5s sign %+d bound %2d: targets %8d  not excluded %8d | not-excluded efficient strata by dim %s'
                  % (k, sg, b, tgt.sum(), tleft.sum(), hist_left.tolist()))
    return res
if __name__ == '__main__':
    report(sys.argv[1])
