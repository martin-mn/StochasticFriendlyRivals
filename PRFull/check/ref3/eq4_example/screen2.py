import numpy as np
from we_core import *
Q = all_binary()
names=['CC','CD','DC','DD']
for x in [0.39, 0.40, 0.6, 0.9]:
    p = list(P_EX); p[9] = x
    m, w = margins(p, Q, 1e-12)
    mn = m.min()
    ties = np.where(m < mn + 1e-9)[0]
    print(f"x={x} min={mn:.12f} n_ties={len(ties)}")
    bits = Q[ties].astype(int)
    # which coordinates are fixed across ties
    fixed = ''.join(str(bits[0,s]) if (bits[:,s]==bits[0,s]).all() else '*' for s in range(16))
    print("  pattern over ties (co-player's own view, s=0..15):", fixed)
    print("  q6 (co-player at its <CD,DC>, focal offer state) values:", set(bits[:,6]))
