import numpy as np, sys, time
from we_core import *
Q = all_binary()
t=time.time()
for x in [0.25, 5/13, 0.40, 0.45, 0.5, 0.0, 1.0]:
    p = list(P_EX); p[9] = x
    for eps in [1e-9, 1e-12]:
        m, w = margins(p, Q, eps)
        i = np.argsort(m)
        print(f"x={x:.5f} eps={eps:g} min={m[i[0]]:.12f}  worst codes:", [(int(c), round(m[c],9)) for c in i[:8]])
    # self-play
    T = transition_batch(p, np.array([p]), 1e-12); pi = gth_b(T); print("   selfplay w:", pi.reshape(4,4).sum(1).round(9))
print(f"{time.time()-t:.1f} s")
