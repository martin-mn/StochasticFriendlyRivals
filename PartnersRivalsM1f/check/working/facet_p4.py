"""Facet {p4=0}, p1,p2,p3 interior.  Certify that all 16 limit-Nash rows are implied, at EVERY
game (u,v), by the two rows of ALLC (1111) and 0001, and that those two are positive multiples of
   (I*)  (1-p1) u + p3 <= 0        (II*) (1-p2) u + p3 (1+v) <= 0
which is the C/D mirror of (I),(II).  Master form a(u-1)+b v <= c."""
from fractions import Fraction as F
import random
from limlib import nash_rows, DET

random.seed(23)
G = {''.join(str(int(x)) for x in s): k for k, s in enumerate(DET)}
rnd = lambda: F(random.randint(1, 59), 60)

bad = ties = n = 0
for trial in range(40):
    p1, p2, p3 = rnd(), rnd(), rnd()
    p = (p1, p2, p3, F(0))
    rows, Wp = nash_rows(p)
    assert Wp == [0, 0, 0, 1]
    rI, rII = rows[G['1111']], rows[G['0001']]
    I = (1 - p1, F(0), -(p3 + 1 - p1))
    II = (1 - p2, p3, -(1 - p2 + p3))
    for r, g in ((rI, I), (rII, II)):        # positive multiple?
        assert r[0] > 0 and r[0] * g[1] == r[1] * g[0] and r[0] * g[2] == r[2] * g[0], (p, r, g)
    for k, (a, b, c) in enumerate(rows):
        n += 1
        if (a, b, c) == (0, 0, 0):
            ties += 1
            continue
        mu = b / rII[1]                       # rI has b = 0
        lam = (a - mu * rII[0]) / rI[0]
        nu = lam * rI[2] + mu * rII[2] - c
        if not (mu >= 0 and lam >= 0 and nu >= 0):
            bad += 1
            print('NO CERT', p, ''.join(str(int(x)) for x in DET[k]), (a, b, c), (lam, mu, nu))
print(f'facet p4=0: 40 random interior p, {n} rows, {ties} identically-zero rows, '
      f'{bad} rows without a Farkas certificate')
