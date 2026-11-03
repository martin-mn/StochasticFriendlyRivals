"""greedy minimum cover of the certificate triples by co-players (from `ub sets`)"""
import numpy as np, sys
qs = [int(x) for x in open('runs/tried.txt').read().split()]
T = np.fromfile('runs/tri_full.bin', dtype=np.int64).reshape(-1, 3); n = len(T); nb = (n + 7)//8
B = np.unpackbits(np.fromfile('runs/sets.bin', dtype=np.uint8).reshape(len(qs), nb), axis=1, bitorder='little')[:, :n].astype(bool)
assert B.any(0).all(), 'some triple not covered by any tried q'
left = np.ones(n, bool); chosen = []
while left.any():
    gain = (B & left).sum(1); i = int(gain.argmax()); chosen.append((qs[i], int(gain[i]))); left &= ~B[i]
print('greedy cover:', len(chosen), 'co-players'); print(chosen)
# per sign
for sg in (1, -1):
    m = T[:, 1] == sg; l = m.copy(); ch = []
    while l.any():
        gain = (B & l).sum(1); i = int(gain.argmax()); ch.append(qs[i]); l &= ~B[i]
    print('sign %+d: %d targets, greedy cover %d co-players: %s' % (sg, m.sum(), len(ch), ch))
open('runs/qstar.txt', 'w').write('\n'.join(str(q) for q, _ in chosen) + '\n')
# assign each triple its first chosen covering q
Q = np.array([q for q, _ in chosen]); idxs = [qs.index(q) for q in Q]
assign = np.full(n, -1, np.int64)
for q, i in zip(Q, idxs):
    s = (assign < 0) & B[i]; assign[s] = q
T2 = T.copy(); T2[:, 2] = assign; T2.tofile('runs/tri_star.bin'); print('wrote runs/qstar.txt and runs/tri_star.bin')
