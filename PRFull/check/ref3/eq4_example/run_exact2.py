import sympy as sp
from exact import *
R = sp.Rational
p = list(EX)  # p9 = 1/4
for c in [64, 65, 0, 52620]:
    m, w = margin_limit(p, code2q(c))
    print(f"p9=1/4 co-player {c} q6={code2q(c)[6]}: margin={m}, w={w}")
for xv in [0, 1]:
    pp = list(EX); pp[9] = xv
    m, w = margin_limit(pp, pp)
    print(f"p9={xv} self-play w = {w}")
    for c in [52620, 58764]:
        m, w = margin_limit(pp, code2q(c))
        print(f"   p9={xv} vs {c}: margin {m}")
