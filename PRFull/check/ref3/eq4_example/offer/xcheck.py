"""Cross-checks with the earlier, independent code in check/ref3/eq4_example:
(1) float screen (we_core, eps=1e-12) just above 5/13: how many co-players have a negative margin;
(2) sympy exact limits over Q(x, eps) (exact.py) for three binders that are not among the eight."""
import sys, time
sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/check/ref3/eq4_example')
import numpy as np
from we_core import all_binary, margins, P_EX
Q = all_binary()
for x9 in [5/13 - 1e-3, 5/13 + 1e-3, 0.4]:
    p = list(P_EX); p[9] = x9
    m, w = margins(p, Q, 1e-12)
    print('float: x = %.6f  co-players with margin < -1e-9: %d   min %.3e' % (x9, (m < -1e-9).sum(), m.min()))
import sympy as sp
from exact import EX, code2q as c2q, x, margin_limit
p = list(EX); p[9] = x
t = time.time()
for c in [33153, 33920, 49536]:
    mm, w = margin_limit(p, c2q(c))
    mm = sp.factor(mm)
    print('sympy exact: co-player', c, 'margin(x) =', mm, ' root:', sp.solve(sp.numer(sp.together(mm)), x))
print('%.1fs' % (time.time() - t))
