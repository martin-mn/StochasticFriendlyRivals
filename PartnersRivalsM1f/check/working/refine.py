"""(1) the degenerate game (1,0);  (2) the boundary of the dim-3 region, checked with a fine
sweep near u=1 and v=0;  (3) the 16 vertices' Nash count on the same games."""
from fractions import Fraction as F
import random
from limlib import nash_rows, DET, limit_W

rnd = lambda r: F(r.randint(1, 199), 200)

# ---------- (1) the whole facet p1=1 is (weakly) limit-Nash at the degenerate game (1,0)
r = random.Random(7)
bad = 0
for _ in range(200):
    p = (F(1), rnd(r), rnd(r), rnd(r))
    rows, _ = nash_rows(p)
    if not all(a * (F(1) - 1) + b * F(0) <= c for a, b, c in rows):
        bad += 1
print(f'(1,0): 200 random points of the open facet p1=1, {bad} not limit-Nash  '
      f'-> the Nash set is 3-dimensional at (1,0)')
# and the open cube is not Nash there
bad = sum(0 if all(a*(F(1)-1)+b*F(0) <= c for a, b, c in nash_rows(
              (rnd(r), rnd(r), rnd(r), rnd(r)))[0]) else 1 for _ in range(60))
print(f'      60 random interior points, {bad} not limit-Nash -> no 4-dimensional piece at (1,0)')

# ---------- (2) boundary of the dim-3 region
def full3(spec, u, v, n=400, seed=3):
    rr = random.Random(seed)
    for _ in range(n):
        p = tuple(rnd(rr) if c == 'i' else F(c) for c in spec)
        rows, _ = nash_rows(p)
        nz = [r for r in rows if r != (0, 0, 0)]
        if all(a*(u-1) + b*v < c for a, b, c in nz):
            return True
    return False

print('\nboundary of the dim-3 region (strictly full-dimensional inside the facet):')
tests = [(F(1), F(-1,1000)), (F(1), F(0)), (F(1), F(1,1000)), (F(1001,1000), F(0)),
         (F(999,1000), F(0)), (F(3), F(-1,1000)), (F(3), F(0)), (F(1), F(-1,10**6))]
for u, v in tests:
    a = full3(('i','i','i',0), u, v, 200)
    b = full3((1,'i','i','i'), u, v, 200)
    print(f'  u={str(u):>12s} v={str(v):>12s}   p4=0 facet full: {str(a):5s}  p1=1 facet full: {b}')

# ---------- (3) the sixteen
print('\nthe sixteen: rows once, then the count at any game')
VROWS = {}
for s in DET:
    g = ''.join(str(int(x)) for x in s)
    VROWS[g] = nash_rows(s)[0]
import pickle
pickle.dump(VROWS, open('vrows.pkl', 'wb'))
def count16(u, v, strict=False):
    out = []
    for g, rows in VROWS.items():
        nz = [r for r in rows if r != (0, 0, 0)]
        ok = all(a*(u-1)+b*v <= c for a, b, c in nz)
        if strict:
            ok = all(a*(u-1)+b*v < c for a, b, c in nz)
        if ok:
            out.append(g)
    return out
for u, v in [(-3,-3),(-1,-1),(0,-1),(F(1,2),-1),(2,-1),(-1,F(1,2)),(F(1,2),F(1,2)),
             (2,F(1,2)),(-1,2),(F(1,2),2),(2,2),(1,0),(3,1),(1,2),(4,1),(0,0),(1,1)]:
    u, v = F(u), F(v)
    w = count16(u, v)
    print(f'  ({str(u):>4s},{str(v):>4s})  N={len(w):2d}  {sorted(w)}')
