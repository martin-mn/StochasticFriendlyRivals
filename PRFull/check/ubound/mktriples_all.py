"""by-product: the (stratum, sign, q) triples of EVERY excluded efficient stratum of W, N, E and the
T>S ray, all dimensions (run with UB_BOUNDS='{"W":0,"W68":0,"N":0,"E":0,"ray+":0}'), in two halves"""
import sys, numpy as np, status as st
prefix, out = sys.argv[1], sys.argv[2]
idx, qp, qm = st.load(prefix); M, d = st.masks(idx)
need = {1: np.zeros(len(idx), bool), -1: np.zeros(len(idx), bool)}
for k, (m, sg, b) in M.items():
    if b > 15: continue
    need[sg] |= m & (d >= b) & ((qp if sg > 0 else qm) >= 0)
T = []
for sg, q in ((1, qp), (-1, qm)):
    t = np.where(need[sg])[0]; T.append(np.stack([idx[t], np.full(len(t), sg), q[t].astype(np.int64)], 1))
T = np.concatenate(T).astype(np.int64); T = T[np.random.default_rng(13).permutation(len(T))]
h = len(T)//2; T[:h].tofile(out + '_A.bin'); T[h:].tofile(out + '_B.bin'); print(len(T), 'triples in two halves')
