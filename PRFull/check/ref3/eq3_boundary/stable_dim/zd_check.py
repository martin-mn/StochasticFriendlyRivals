"""The equaliser plane E(u,v): p_CD = p_CC - v(1-p_CC+p_DD), p_DC = p_DD - u(1-p_CC+p_DD).
Check with exact_partners.limitW that every point of E inside the cube makes all 16 binary deviators tie, and that
points just off E do not.  Also check the facet points of E (one coordinate 0 or 1)."""
import os
import sys, random
from fractions import Fraction as F
sys.dont_write_bytecode = True
# exact_partners.py of the memory-one project: PartnersRivalsM1f/check, found from PRFull/check/ref3/eq3_boundary/<dir>/
_M1_CAND = [os.path.join(os.path.dirname(os.path.abspath(__file__)), *(['..'] * 5), 'PartnersRivalsM1f', 'check'),
            '/Users/martin/Documents/CL1/PartnersRivalsM1f/check']
sys.path.insert(0, next(d for d in _M1_CAND if os.path.isfile(os.path.join(d, 'exact_partners.py'))))
import exact_partners as ep
rng = random.Random(7)

def pays(p, u, v):
    out = []
    for q in ep.QS:
        W, _, ok = ep.limitW(tuple(F(x) for x in q), p, 1)
        assert ok
        out.append(W[0] + W[1] * u + W[2] * (1 + v))
    return out

n_in = n_off = 0
tries = 0
while n_in < 200:
    tries += 1
    u = F(rng.randint(-200, 200), 100); v = F(rng.randint(-200, 200), 100)
    pcc = F(rng.randint(1, 99), 100); pdd = F(rng.randint(1, 99), 100)
    k = 1 - pcc + pdd
    pcd = pcc - v * k; pdc = pdd - u * k
    if not (0 < pcd < 1 and 0 < pdc < 1):
        continue
    p = (pcc, pcd, pdc, pdd)
    P = pays(p, u, v)
    assert len(set(P)) == 1, (p, u, v, P)
    n_in += 1
    p2 = (pcc, pcd + F(1, 1000), pdc, pdd)
    if 0 < p2[1] < 1:
        assert len(set(pays(p2, u, v))) > 1
        n_off += 1
print('interior points of E: %d, all 16 deviators tie exactly; %d nearby points off E: not all tie' % (n_in, n_off))
# facet points of E: e.g. p_CC = 0 (F3) and p_DC = 1 (F8) etc.
n_f = 0
for _ in range(4000):
    u = F(rng.randint(-200, 200), 100); v = F(rng.randint(-200, 200), 100)
    pcc = F(rng.choice([0, 1, rng.randint(1, 99)]), 1) if False else rng.choice([F(0), F(1), F(rng.randint(1, 99), 100)])
    pdd = rng.choice([F(0), F(1), F(rng.randint(1, 99), 100)])
    k = 1 - pcc + pdd
    if k == 0:
        continue
    pcd = pcc - v * k; pdc = pdd - u * k
    if not (0 <= pcd <= 1 and 0 <= pdc <= 1):
        continue
    p = (pcc, pcd, pdc, pdd)
    if sum(1 for x in p if x in (0, 1)) != 1:
        continue
    P = pays(p, u, v)
    n_f += 1
    if len(set(P)) != 1:
        print('facet point of E that is NOT an equaliser in the limit:', p, u, v, sorted(set(P))[:3])
print('facet points of E tested: %d' % n_f)
