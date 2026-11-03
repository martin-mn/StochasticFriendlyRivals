"""Eq. (3) on the efficient edges with p_DD = 0, {p_CC=1, p_DD=0, max(p_CD,p_DC)=1} off the segments, i.e. the
strata (1,a,1,0), (1,1,b,0), 1110: symbolic Farkas certificate (=>) and witnesses (<=), as certify.py."""
import sympy as sp
from certify import gains, farkas, affine, sign_of, NAMES, X, u, v

for st in ['1o10', '11o0', '1110']:
    free = [i for i, ch in enumerate(st) if ch == 'o']
    p = tuple(sp.Integer(0) if ch == '0' else sp.Integer(1) if ch == '1' else X[i] for i, ch in enumerate(st))
    pc, a, b, c = p
    V, G = gains(p)
    I = affine(c * v - (1 - a))
    II = affine((1 - a) * (u - 1) + b * v)
    cons = [z for z in (I, II) if any(sp.simplify(w) != 0 for w in z[:2])]
    probs = []
    for name, Gs in zip(NAMES, G):
        if farkas(Gs, cons, free) is None:
            probs.append(name)
    wit = []
    for g in cons:
        for name, Gs in zip(NAMES, G):
            rs = [sp.simplify(Gs[j] / g[j]) for j in range(3) if g[j] != 0]
            if all(sp.simplify(Gs[j]) == 0 for j in range(3) if g[j] == 0) and rs and \
                    all(sp.simplify(r - rs[0]) == 0 for r in rs) and sign_of(rs[0], free) == 1:
                wit.append((name, str(sp.factor(rs[0]))))
                break
    print(st, 'V', V, 'binding constraints', cons, 'Farkas failures', probs, 'witnesses', wit,
          'I constant sign', sign_of(sp.Add(*[I[2]]), free) if I[:2] == (0, 0) else 'n/a')
