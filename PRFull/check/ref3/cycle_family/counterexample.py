"""(1) The evaluation's counterexample, exactly: p=(0,0,1,1), (u,v)=(-1/4,-1/2).
(a) eps as a free symbol: exact stationary distributions of self-play and of 0101 vs p.
(b) eps->0 limit of all 16 gains (sympy cofactors) and, independently, with the integer-arithmetic
    limitW of PartnersRivalsM1f/check/exact_partners.py."""
import sys
sys.path.insert(0, '.')
sys.path.insert(0, '/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
from fractions import Fraction as F
from cyclelib import *
import exact_partners as XP

p = (0, 0, 1, 1)
U, V = sp.Rational(-1, 4), sp.Rational(-1, 2)
for sig, nm in [(p, 'self-play'), ((0, 1, 0, 1), '0101 vs p')]:
    ds = tree_weights(sig, p)
    w = [sp.factor(d / sum(ds)) for d in ds]
    print(nm, 'exact w(eps) =', w)
# evaluation's formulas for 0101
ds = tree_weights((0, 1, 0, 1), p); w = [sp.cancel(d / sum(ds)) for d in ds]
print('0101: w_CC - (1/2 - e + e^2) =', sp.simplify(w[0] - (sp.Rational(1, 2) - e + e**2)),
      '; w_CD - (e - e^2) =', sp.simplify(w[1] - (e - e**2)))
ds = tree_weights(p, p); w = [sp.cancel(d / sum(ds)) for d in ds]
print('self-play uniform for every eps:', all(sp.simplify(x - sp.Rational(1, 4)) == 0 for x in w))

Wp, _, _ = limit_on_stratum(tuple(map(sp.Integer, p)), tuple(map(sp.Integer, p)), [])
pip = payoff(Wp).subs({u: U, v: V})
print('W(p,p) =', Wp, ' pi(p,p) =', pip)
Vx, _, _ = XP.limitW([60 * x for x in p], [60 * x for x in p])
for sig in BIN:
    W, _, _ = limit_on_stratum(tuple(map(sp.Integer, sig)), tuple(map(sp.Integer, p)), [])
    g = payoff(W).subs({u: U, v: V}) - pip
    Wx, _, _ = XP.limitW([60 * x for x in sig], [60 * x for x in p])
    gx = (Wx[0] + F(-1, 4) * Wx[1] + F(1, 2) * Wx[2]) - (Vx[0] + F(-1, 4) * Vx[1] + F(1, 2) * Vx[2])
    assert [sp.Rational(x.numerator, x.denominator) for x in Wx] == W, (sig, W, Wx)
    print('  ', label(sig), 'W =', W, ' gain =', g, ' (exact_partners: %s)' % gx)
