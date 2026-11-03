"""Claim D1 of certify.py, which its automatic Farkas search leaves open: on the switch line u+v=1 the strategies
(1,1,0,y), 0<y<=1, are stable only at (u,v)=(1,0).  Exact eps->0 limits (m1rows.py), the gains as functions of y and u
with v=1-u: every gain is a multiple of (u-1); ALLD gains 1-u and ALLC gains (u-1)/(y+2), so for u<1 ALLD and for
u>1 ALLC gains, and at u=1 every gain is 0.  (The search fails only because the coefficient (2y-1)/(y+4) of 0111
changes sign at y=1/2, which does not matter where all gains vanish.)"""
import itertools
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'farkas_display'))
from m1rows import stationary_limit  # noqa: E402

y, u = sp.symbols('y u', positive=True)
v = 1 - u
p = (sp.Integer(1), sp.Integer(1), sp.Integer(0), y)


def pay(W):
    return W[0] + u * W[1] + (1 + v) * W[2]


Wpp = stationary_limit(p, p)
ok = True
for bits in itertools.product([0, 1], repeat=4):
    W = stationary_limit(tuple(sp.Integer(x) for x in bits), p)
    g = sp.factor(sp.simplify(pay(W) - pay(Wpp)))
    ok &= sp.simplify(g.subs(u, 1)) == 0
    print(''.join(map(str, bits)), g)
gd = sp.simplify(pay(stationary_limit((0, 0, 0, 0), p)) - pay(Wpp))
gc = sp.simplify(pay(stationary_limit((1, 1, 1, 1), p)) - pay(Wpp))
ok &= sp.simplify(gd - (1 - u)) == 0 and sp.simplify(gc - (u - 1) / (y + 2)) == 0
print('D1 proved: every gain vanishes at u=1; ALLD gains 1-u, ALLC (u-1)/(y+2):', ok)
