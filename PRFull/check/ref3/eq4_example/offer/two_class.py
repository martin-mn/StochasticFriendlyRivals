"""The binders whose pair chain at eps=0 has two closed classes: accounting in the class of s*."""
from fractions import Fraction as F
from lo import code2q, pstar, limit
from renewal import closed_classes, step0, SSTAR, st, NAMES, excursion
import renewal

def load(name):
    return {int(c): F(m) for c, m in (l.split() for l in open(name))}
m513, m14 = load('margins_5_13.txt'), load('margins_1_4.txt')
bind = sorted(c for c in m513 if m513[c] == 0 and m14[c] != 0)
p = pstar(F(1, 4))
summary = {}
for c in bind:
    q = code2q(c)
    cls = closed_classes(p, q)
    if len(cls) == 1:
        continue
    cs = next(cl for cl in cls if SSTAR in cl)
    other = [cl for cl in cls if SSTAR not in cl]
    # excursion within the class of s*: reuse renewal.excursion by restricting (closed_classes picks the class with s*)
    orig = renewal.closed_classes
    renewal.closed_classes = lambda p_, q_: [cs]
    _, res = excursion(p, q)
    renewal.closed_classes = orig
    key = (tuple((str(pr), str(g)) for pr, _, g, _ in res), tuple(tuple(st(s) for s in cl) for cl in other))
    summary.setdefault(key, []).append(c)
for k, v in summary.items():
    print(len(v), 'per-branch (prob, co-player net):', k[0], ' other class:', k[1])
print('total', sum(len(v) for v in summary.values()))
c = summary[next(iter(summary))][0]
print('example', c, 'W at x=1/4:', limit(p, code2q(c)))
