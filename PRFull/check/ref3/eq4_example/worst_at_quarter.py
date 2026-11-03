import numpy as np
from we_core import *
Q = all_binary()
p = list(P_EX)
m, w = margins(p, Q, 1e-12)
mn = m.min(); print("min margin at p9=1/4 (eps=1e-12):", mn)
ties = np.where(m < 1e-9)[0]
print("co-players within 1e-9 of 0:", len(ties))
q6 = Q[ties, 6].astype(int)
print("  of these, q6=0 (defect against offer):", (q6==0).sum(), " q6=1:", (q6==1).sum())
fam = [52620, 36236, 56716, 40332, 52612, 36228, 56708, 40324]
print("family margins at 1/4:", m[fam], " 7/174 =", 7/174)
# rank of family at 1/4
order = np.sort(m)
print("number of co-players with margin < 7/174 - 1e-9:", (m < 7/174 - 1e-9).sum())
# examples of minimizers
print("examples q6=1:", ties[q6==1][:5], " q6=0:", ties[q6==0][:5])
np.save('ties_quarter.npy', ties)
# what fraction of minimizers end in mutual cooperation
print("minimizers with w_CC>0.999:", (w[ties,0]>0.999).sum())
