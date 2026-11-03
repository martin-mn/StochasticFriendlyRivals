"""Independent cross-checks of prove.py, with the integer/Fraction code of PartnersRivalsM1f
(exact_partners.limitW, m1lib.stationary_exact), which shares no code with cyclelib.

A. Grid: p = (0,a,b,1), a,b in the 13 values k/n, n <= 6; games (u,v) = (i/12, j/12), |i|,|j| <= 24
   (2401 games, including every line a = -2v, b = 1-2u through grid points, and both corner games).
   Exact limit-Nash test against all 16 binary deviators vs the CORRECTED predicate and vs the
   ORIGINAL sentence.
B. Near the corners: a, b at distance 10^-k (k = 3, 6, 9) from (0,1) and (1,0): exact limits obey
   the open-square formula (the discontinuity sits at the corner point only).
C. Finite eps: at the two exceptional corners, exact stationary distributions at eps = 10^-4,
   10^-6, 10^-8 give gains that converge to the limit gains at rate O(eps).
"""
import sys
sys.path.insert(0, '/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
from fractions import Fraction as F
from itertools import product
import exact_partners as XP
import m1lib

QS = list(product([0, 1], repeat=4))


def limits(p, one):
    pa = [one * x for x in p]
    if one != 1:
        pa = [int(x) for x in pa]
    V = XP.limitW(pa, pa, one)[0]
    Ws = [XP.limitW([one * x for x in q], pa, one)[0] for q in QS]
    return V, Ws


def gains(V, Ws, u, v):
    pv = V[0] + u * V[1] + (1 + v) * V[2]
    return [W[0] + u * W[1] + (1 + v) * W[2] - pv for W in Ws]


def corrected(a, b, u, v):
    if (a, b) == (0, 1):
        return (u, v) == (0, 0)
    if (a, b) == (1, 0):
        return (u, v) == (F(1, 2), F(-1, 2))
    return a <= -2 * v and b <= 1 - 2 * u


def original(a, b, u, v):
    return a <= -2 * v and b <= 1 - 2 * u


GRID = sorted({F(k, n) for n in range(1, 7) for k in range(n + 1)})
GAMES = [(F(i, 12), F(j, 12)) for i in range(-24, 25) for j in range(-24, 25)]
nt = bad_corr = bad_orig = 0
orig_fail = []
for a in GRID:
    for b in GRID:
        p = (F(0), a, b, F(1))
        assert all((60 * x).denominator == 1 for x in p)
        V, Ws = limits(p, 60)
        for (u, v) in GAMES:
            st = max(gains(V, Ws, u, v)) <= 0
            nt += 1
            if st != corrected(a, b, u, v):
                bad_corr += 1
            if st != original(a, b, u, v):
                bad_orig += 1
                orig_fail.append((a, b, u, v))
print('A. %d exact (strategy, game) tests: corrected predicate disagrees %d times; original sentence %d times'
      % (nt, bad_corr, bad_orig))
print('   original fails only at (a,b) in', sorted({(str(x[0]), str(x[1])) for x in orig_fail}))
for ab in sorted({(x[0], x[1]) for x in orig_fail}):
    fl = [(x[2], x[3]) for x in orig_fail if (x[0], x[1]) == ab]
    print('   (a,b)=(%s,%s): %d games where the sentence says stable but p is not, e.g. (u,v)=%s' % (ab[0], ab[1], len(fl), [(str(x), str(y)) for x, y in fl[:3]]))

# B. near the corners, Fraction coordinates (one = 1)
nb = 0
for k in (3, 6, 9):
    d = F(1, 10 ** k)
    for (a, b) in [(d, 1 - d), (1 - d, d), (F(0), 1 - d), (d, F(1)), (F(1), d), (1 - d, F(0))]:
        p = (F(0), a, b, F(1))
        V, Ws = limits(p, 1)
        assert V == (F(1, 2), 0, 0, F(1, 2)), V
        u, v = F(-3, 7), F(-2, 9)
        g = gains(V, Ws, u, v)
        # ALLD (index 0) and ALLC (index 15) gains: (a+2v)/(2(2-a)) and (b+2u-1)/(2(1+b))
        assert g[0] == (a + 2 * v) / (2 * (2 - a)), (a, b, g[0])
        assert g[15] == (b + 2 * u - 1) / (2 * (1 + b)), (a, b, g[15])
        assert (max(g) <= 0) == original(a, b, u, v)
        nb += 1
print('B. %d near-corner strategies: self-play (1/2,0,0,1/2), ALLD/ALLC gains as on the open square' % nb)

# C. finite eps at the exceptional corners
for p, (u, v) in [((0, 0, 1, 1), (F(-1, 4), F(-1, 2))), ((0, 1, 0, 1), (F(-1, 4), F(-1, 2))),
                  ((0, 1, 0, 1), (F(0), F(-1)))]:
    V, Ws = limits(p, 1)
    gl = gains(V, Ws, u, v)
    best = max(range(16), key=lambda i: gl[i])
    line = []
    for kk in (4, 6, 8):
        eps = F(1, 10 ** kk)
        ws = m1lib.stationary_exact([F(x) for x in p], [F(x) for x in p], eps)
        wq = m1lib.stationary_exact([F(x) for x in QS[best]], [F(x) for x in p], eps)
        pv = ws[0] + u * ws[1] + (1 + v) * ws[2]
        pq = wq[0] + u * wq[1] + (1 + v) * wq[2]
        line.append('eps=1e-%d: gain-limit = %.3e' % (kk, float(pq - pv - gl[best])))
    print('C. p=%s (u,v)=(%s,%s): best deviator %s, limit gain %s;' % (
        ''.join(map(str, p)), u, v, ''.join(map(str, QS[best])), gl[best]), '; '.join(line))
