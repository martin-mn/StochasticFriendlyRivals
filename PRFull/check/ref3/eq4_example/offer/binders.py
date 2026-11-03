"""Which co-players bind at 5/13: margin 0 at x = 5/13 and negative above. Their excursion accounting,
their Moebius margin functions m(x) = (a + b x)/(c + d x), and prime patterns of the set."""
import itertools
from fractions import Fraction as F
from lo import code2q, genome, pstar, margin, M
from renewal import excursion, closed_classes, SSTAR, st

def load(name):
    return {int(c): F(m) for c, m in (l.split() for l in open(name))}

m513, m14, m25, m1100 = (load('margins_%s.txt' % s) for s in ('5_13', '1_4', '2_5', '1_100'))
newzero = sorted(c for c in m513 if m513[c] == 0 and m14[c] != 0)
neg25 = sorted(c for c in m25 if m25[c] < 0)
print('zero at 5/13 but not at 1/4:', len(newzero), '; negative at 2/5:', len(neg25), '; same set:', newzero == neg25)
zero_both = [c for c in m513 if m513[c] == 0 and m14[c] == 0]
print('zero at both 1/4 and 5/13:', len(zero_both), '; of these negative at 2/5:', sum(m25[c] < 0 for c in zero_both))
EIGHT = [36228, 36236, 40324, 40332, 52612, 52620, 56708, 56716]

def moebius(c):
    """reconstruct m(x) = (a+bx)/(g+dx) from three exact values (margin numerator and normaliser are affine in x)"""
    xs = [F(1, 4), F(2, 5), F(1, 100)]
    ms = [m14[c], m25[c], m1100[c]]
    # unknowns a, b, g, d with a + b x - m g - m d x = 0; fix g or d via null space
    rows = [[F(1), x, -m, -m * x] for x, m in zip(xs, ms)]
    import sympy as sp
    N = sp.Matrix(rows).nullspace()
    assert len(N) == 1
    v = N[0] / (N[0][2] if N[0][2] != 0 else N[0][3])
    a, b, g, d = [F(str(t)) for t in v]
    # check against 5/13 and a fourth point
    x4 = F(7, 10)
    assert (a + b * F(5, 13)) / (g + d * F(5, 13)) == m513[c]
    assert (a + b * x4) / (g + d * x4) == margin(pstar(x4), code2q(c))
    return a, b, g, d

groups = {}
for c in neg25:
    a, b, g, d = moebius(c)
    root = -a / b
    assert root == F(5, 13), (c, root)
    p = pstar(F(1, 4)); q = code2q(c)
    cls = closed_classes(p, q)
    key = (tuple(sorted(tuple(cl) for cl in cls)),)
    info = None
    if len(cls) == 1 and SSTAR in cls[0]:
        _, res = excursion(p, q)
        info = tuple((str(pr), str(gn), str(ln)) for pr, _, gn, ln in res)
    groups.setdefault((len(cls), info, q[M[SSTAR]]), []).append(c)

print('\nall 400 have numerator root exactly 5/13')
for k, v in sorted(groups.items(), key=lambda t: -len(t[1])):
    print('n=%d  closed classes=%d  q at offer=%d  excursion (prob, co-player net, rounds) at x=1/4: %s'
          % (len(v), k[0], k[2], k[1]))
    print('     e.g.', v[:6], ' margin at 2/5:', sorted({m25[c] for c in v})[:6])

# rank by margin above 5/13: which co-player is the most negative at several x in (5/13, 1)?
mob = {c: moebius(c) for c in neg25}
for x in [F(2, 5), F(1, 2), F(3, 4), F(99, 100)]:
    vals = {c: (a + b * x) / (g + d * x) for c, (a, b, g, d) in mob.items()}
    mn = min(vals.values())
    arg = sorted(c for c in vals if vals[c] == mn)
    print('x = %s: most negative margin %s by %s' % (x, mn, arg))

# prime patterns of the 400 (Quine-McCluskey merge)
S = {genome(code2q(c)) for c in neg25}
cur = set(S); primes = set()
while cur:
    nxt = set(); used = set()
    cl = list(cur)
    byk = {}
    for g in cl:
        for i in range(16):
            if g[i] != '*':
                k = g[:i] + '_' + g[i + 1:]
                byk.setdefault(k, []).append(g)
    for k, gs in byk.items():
        if len(gs) == 2:
            i = k.index('_')
            nxt.add(k.replace('_', '*')); used.update(gs)
    primes |= cur - used
    cur = nxt
print('\n%d prime patterns of the 400:' % len(primes))
for pat in sorted(primes, key=lambda s: (-s.count('*'), s)):
    print('  ', pat, pat.count('*'))
