"""Facet {p1=1}, p2,p3,p4 interior.  Certify that all 16 limit-Nash rows are implied, at EVERY
game (u,v) of the plane, by the two rows of ALLD (0000) and 0111, and that those two rows are
positive multiples of
   (I)  p4 v <= 1 - p2        (II) (1-p2)(u-1) + p3 v <= 0.
Exact Farkas: L_sigma = lam*L_I + mu*L_II - nu with lam,mu,nu >= 0."""
from fractions import Fraction as F
import random
from limlib import nash_rows, DET

random.seed(11)
G = {''.join(str(int(x)) for x in s): k for k, s in enumerate(DET)}
rnd = lambda: F(random.randint(1, 59), 60)

bad = ties = n = 0
for trial in range(40):
    p2, p3, p4 = rnd(), rnd(), rnd()
    p = (F(1), p2, p3, p4)
    rows, Wp = nash_rows(p)
    assert Wp == [1, 0, 0, 0]
    rI, rII = rows[G['0000']], rows[G['0111']]
    # the two rows are positive multiples of (I) and (II)
    assert rI[0] == 0 and rI[1] > 0 and rI[1] * (1 - p2) == rI[2] * p4, (p, rI)
    assert rII[2] == 0 and rII[0] > 0 and rII[0] * p3 == rII[1] * (1 - p2), (p, rII)
    for k, (a, b, c) in enumerate(rows):
        n += 1
        if (a, b, c) == (0, 0, 0):
            ties += 1
            continue
        mu = a / rII[0]
        lam = (b - mu * rII[1]) / rI[1]
        nu = lam * rI[2] + mu * rII[2] - c
        if not (mu >= 0 and lam >= 0 and nu >= 0):
            bad += 1
            print('NO CERT', p, ''.join(str(int(x)) for x in DET[k]), (a, b, c), (lam, mu, nu))
print(f'facet p1=1: 40 random interior p, {n} rows, {ties} identically-zero rows, '
      f'{bad} rows without a Farkas certificate')
