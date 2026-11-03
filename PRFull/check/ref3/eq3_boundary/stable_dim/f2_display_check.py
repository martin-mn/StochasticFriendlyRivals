"""The proposed SI display of the 16 gains on the open facet p_DD = 0, checked against exact_partners.limitW
(Fractions, a separate code base) at random rational points and games: gain = pi(sigma,p) - pi(p,p), pi(p,p) = 0."""
import os
import random
import sys
from fractions import Fraction as F
sys.dont_write_bytecode = True
_M1_CAND = [os.path.join(os.path.dirname(os.path.abspath(__file__)), *(['..'] * 5), 'PartnersRivalsM1f', 'check'),
            '/Users/martin/Documents/CL1/PartnersRivalsM1f/check']
sys.path.insert(0, next(d for d in _M1_CAND if os.path.isfile(os.path.join(d, 'exact_partners.py'))))
import exact_partners as ep

def disp(s, cc, cd, dc, u, v):
    I = (1 - cc) * u + dc
    J = (1 - cd) * u + dc * (1 + v)
    D = {
        '0000': 0, '1000': 0, '1010': 0, '0010': 0, '0100': 0, '0110': 0,
        '1111': I / (1 - cc + dc), '1101': I / (1 - cc + dc),
        '1100': I / (2 * (1 - cc + dc)), '1110': 2 * I / (3 - 3 * cc + 2 * dc),
        '1001': J / (2 - 2 * cd + dc), '0001': J / (2 - 2 * cd + dc),
        '0011': (cd * I + J) / (2 * (1 - cc * cd + cd * dc)),
        '0101': ((1 - cd) * I + cc * J) / (1 - cd + dc * (2 + cc - 2 * cd)),
        '0111': (I + cc * J) / (1 + 2 * dc - cc * cd),
        '1011': (cd * I + (1 - cc) * J) / (2 - 2 * cc + cd * dc),
    }
    return D[s]

rng = random.Random(11)
n = 0
for _ in range(60):
    cc, cd, dc = (F(rng.randint(1, 96), 97) for _ in range(3))
    u, v = F(rng.randint(-300, 300), 97), F(rng.randint(-300, 300), 97)
    p = (cc, cd, dc, F(0))
    V, _, ok = ep.limitW(p, p, 1)
    assert ok and V == (0, 0, 0, 1)
    for q in ep.QS:
        s = ''.join(map(str, q))
        W, _, ok = ep.limitW(tuple(F(x) for x in q), p, 1)
        g = W[0] + W[1] * u + W[2] * (1 + v)
        assert g == disp(s, cc, cd, dc, u, v), (s, p, u, v)
        n += 1
print('proposed p_DD = 0 display agrees with exact_partners.limitW at %d (point, game, deviator) triples' % n)
