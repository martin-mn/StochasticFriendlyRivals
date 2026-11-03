"""Excursion lengths from s* (in the closed class of s*) for all 400 co-players that bind at 5/13."""
from fractions import Fraction as F
from lo import code2q, pstar
import renewal
from renewal import closed_classes, excursion, SSTAR
def load(name):
    return {int(c): F(m) for c, m in (l.split() for l in open(name))}
m513, m14 = load('margins_5_13.txt'), load('margins_1_4.txt')
bind = sorted(c for c in m513 if m513[c] == 0 and m14[c] != 0)
p = pstar(F(1, 4))
L = {}
for c in bind:
    q = code2q(c)
    cs = next(cl for cl in closed_classes(p, q) if SSTAR in cl)
    orig = renewal.closed_classes; renewal.closed_classes = lambda a, b: [cs]
    _, res = excursion(p, q); renewal.closed_classes = orig
    d = {lab.split(',')[0]: (g, l) for pr, lab, g, l in res}
    assert d['focal C'][0] == 1 and d['focal D'][0] == F(-5, 8)
    L[c] = (d['focal C'][1], d['focal D'][1])
mo = min(v[0] for v in L.values()); mn = min(v[1] for v in L.values())
print(len(L), 'binders; per-branch gains +1 and -5/8 for all; min offer length', mo, 'min no-offer length', mn)
print('attaining both minima:', sorted(c for c, v in L.items() if v == (mo, mn)))
print('attaining offer min:', len([c for c, v in L.items() if v[0] == mo]), '; no-offer min:', sorted(c for c, v in L.items() if v[1] == mn))
