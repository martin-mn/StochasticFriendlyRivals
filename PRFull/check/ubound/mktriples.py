"""the (stratum, sign, q) triples that make up the upper-bound certificate, for `ub verify`"""
import sys, numpy as np, status as st
prefix, out = sys.argv[1], sys.argv[2]; sub = int(sys.argv[3]) if len(sys.argv) > 3 else 0
idx, qp, qm = st.load(prefix); M, d = st.masks(idx)
need = {+1: np.zeros(len(idx), bool), -1: np.zeros(len(idx), bool)}
for k, (m, sg, b) in M.items():
    if b > 15: continue
    need[sg] |= m & (d >= b)
T = []
for sg, q in ((+1, qp), (-1, qm)):
    t = np.where(need[sg])[0]; assert (q[t] >= 0).all(), 'uncovered target'
    T.append(np.stack([idx[t], np.full(len(t), sg), q[t].astype(np.int64)], 1))
T = np.concatenate(T).astype(np.int64)
rng = np.random.default_rng(5); T = T[rng.permutation(len(T))]
if sub: T = T[:sub]
T.tofile(out); print(out, len(T), 'triples; corner evaluations', int((2.0**d[np.searchsorted(idx, T[:,0])]).sum()))
