"""Greedy set cover of the separation strata (c6 + c8 + c68 survivors of the ALLD test) by the
co-players in out/cover/<tag>.q, from the success masks written by cover.c.  Every stratum must be
covered.  Writes out/cover/<tag>_cover.json: the chosen co-players in order and the strata each adds."""
import os, json, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, 'out', 'cover')
tag = sys.argv[1] if len(sys.argv) > 1 else 'sep'
Q = [int(l) for l in open(os.path.join(D, tag + '.q')) if l.strip()]
rows = [l.split('\t') for l in open(os.path.join(D, tag + '.mask'))]
W = np.array([[int(w, 16) for w in r[1:5]] for r in rows], dtype=np.uint64)
B = np.zeros((len(rows), len(Q)), dtype=bool)
for k in range(len(Q)):
    B[:, k] = (W[:, k >> 6] >> np.uint64(k & 63)) & np.uint64(1)
left = np.ones(len(rows), dtype=bool); chosen = []
while left.any():
    gains = B[left].sum(0); k = int(gains.argmax())
    if gains[k] == 0: break
    chosen.append((Q[k], int(gains[k]))); left &= ~B[:, k]
out = dict(strata=len(rows), coplayers_tested=len(Q), uncoverable_by_list=int((~B.any(1)).sum()),
           cover_size=len(chosen), left=int(left.sum()), greedy_cover=chosen,
           strata_per_coplayer_any=dict(zip(map(str, Q), map(int, B.sum(0)))))
json.dump(out, open(os.path.join(D, tag + '_cover.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'strata_per_coplayer_any'}))
