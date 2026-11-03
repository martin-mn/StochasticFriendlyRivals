"""The dimension field, stated, and checked against every receipt in this directory."""
from fractions import Fraction as F
import random
from limlib import nash_rows

def dim_claim(u, v):
    if u < 1 or v < 0: return 3
    if (u, v) == (1, 0): return 3
    return 2

def carrier_p1(u, v):  return u < 1 or v < 0 or (u, v) == (1, 0)
def carrier_p4(u, v):  return u < 0

# targeted witnesses on and around the two boundary rays
random.seed(0)
tests = []
for u, v in [(F(1),F(0)),(F(1),F(-1,10**6)),(F(1),F(1,10**6)),(F(1+F(1,10**6)),F(0)),
             (F(999,1000),F(0)),(F(5),F(-1,10**4)),(F(5),F(0)),(F(2),F(3)),(F(1),F(5))]:
    # best shot on facet p1=1: push p2 -> 1
    ok = False
    for e in (F(1,10),F(1,10**3),F(1,10**6),F(1,10**9)):
        p = (F(1), 1-e, F(1,2), F(1,2))
        rows,_ = nash_rows(p)
        nz = [r for r in rows if r != (0,0,0)]
        if all(a*(u-1)+b*v < c for a,b,c in nz): ok = True; break
    wk = all(a*(u-1)+b*v <= c for a,b,c in nash_rows((F(1),F(1,3),F(1,4),F(1,5)))[0])
    # alternating sheet
    rs,_ = nash_rows((F(1,3),F(0),F(1),F(2,3)))
    alt = all(a*(u-1)+b*v < c for a,b,c in rs if (a,b,c)!=(0,0,0))
    print(f'u={str(u):>16s} v={str(v):>16s}  claim dim {dim_claim(u,v)}   '
          f'p1=1 facet strictly full {str(ok):5s}  whole facet weakly Nash {str(wk):5s}  '
          f'(0,1) sheet strictly full {alt}')
