"""IDs, pattern, exact margins of the eight co-players as functions of x = p_<DC,CD>."""
import time, itertools
from fractions import Fraction as F
from lo import *

IDS = [36228, 36236, 40324, 40332, 52612, 52620, 56708, 56716]
PAT = '001*00011011*0*1'

# 1. the coding convention: ALLD=0, ALLC=65535, TFT and WSLS genomes of PR2026, 19079 in pattern 111*000101*1*010
assert q2code([1] * 16) == 65535 and q2code([0] * 16) == 0
tft = [1 if (s // 4) in (0, 2) else 0 for s in range(16)]   # cooperate iff co-player cooperated last round
wsls = [1 if (s // 4) in (0, 3) else 0 for s in range(16)]
print('TFT genome', genome(tft), 'code', q2code(tft), '| WSLS genome', genome(wsls), 'code', q2code(wsls))

def expand(pat):
    stars = [i for i, ch in enumerate(pat) if ch == '*']
    out = []
    for bits in itertools.product('01', repeat=len(stars)):
        g = list(pat)
        for i, b in zip(stars, bits):
            g[i] = b
        out.append(q2code([int(ch) for ch in g]))
    return sorted(out)

print('pattern', PAT, '->', expand(PAT))
assert expand(PAT) == sorted(IDS)
assert 19079 in expand('111*000101*1*010')
for c in IDS:
    print(c, genome(code2q(c)))
print('co-player state seen when the focal is at <DC,CD> (s=9):', M[9], '; q at it for all eight:',
      {code2q(c)[M[9]] for c in IDS})
print('focal states whose co-player view is a wildcard position:', [s for s in range(16) if M[s] in (3, 12, 14)])

# 2. self-play of p(1/4) and of p(x) at a few x
for x in [F(1, 4), F(1, 100), F(5, 13), F(9, 10)]:
    print('self-play x =', x, limit(pstar(x), pstar(x)))

# 3. margins of the eight at many x; Moebius identity with (5-13x)/(2(25-13x))
xs = [F(k, 97) for k in range(1, 97, 6)] + [F(1, 4), F(5, 13), F(2, 5), F(1, 1000), F(999, 1000)]
t = time.time()
for c in IDS:
    q = code2q(c)
    for x in xs:
        m = margin(pstar(x), q)
        assert m == (5 - 13 * x) / (2 * (25 - 13 * x)), (c, x, m)
    w = limit(pstar(F(1, 4)), q)
    print(c, 'margin at 1/4 =', w[2] - w[1], ' W =', w, ' at 5/13:', margin(pstar(F(5, 13)), q))
print('all %d x-values match (5-13x)/(2(25-13x)) for all eight; %.2fs' % (len(xs), time.time() - t))
print('7/174 =', F(7, 174), ' formula at 1/4:', (5 - 13 * F(1, 4)) / (2 * (25 - 13 * F(1, 4))))
# corners x = 0 and x = 1 (other strata)
for xv in [0, 1]:
    print('corner x =', xv, ' self-play', limit(pstar(xv), pstar(xv)),
          ' margins of the eight', sorted({margin(pstar(xv), code2q(c)) for c in IDS}))
