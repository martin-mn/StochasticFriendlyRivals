"""V12: an explicit PARTNER at every game of the plane, verified in exact Fractions."""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib

def pi(sig,res,eps,u,v):
    w=m1lib.stationary_exact(sig,res,eps)
    return w[0]+w[1]*u+w[2]*(1+v)
def gap(p,u,v,eps):
    return pi(p,p,eps,u,v)-max(pi(s,p,eps,u,v) for s in m1lib.DET_EXACT)
def effdef(p,u,v,eps):
    return pi(p,p,eps,u,v)-max(F(1),(1+u+v)/2)

def altfeas(u,v,a,d):
    return ((d+1)*(1-u)-(1-d)*v<=2 and a*(u-1)+(2-a)*v>=0 and (d+1)*(u-1)+(2-a)*(1+v)>=0
            and (a*d+d+1)*(1-u)+(a*d-d-1)*v<=2 and (a+d+1)*(1-u)+(a+d-3)*v<=2
            and (1+d-a)*(u+v-1)+2*(1-a)>=0 and (a*d-d-1)*(u-1)-(a*d-2*a-d+3)*v<=2*(1-a))

ALTGRID=[(F(0),F(1,2**k)) for k in range(1,14)]+[(1-F(1,2**k),F(1,2**k)) for k in range(1,14)] \
        +[(1-F(1,2**k),F(1,2**j)) for k in range(1,14) for j in range(1,14)]

def partner(u,v):
    """returns (p, tag) -- an efficient limit-Nash strategy for the game (u,v)."""
    if u+v<=1:
        if v<=0: return (F(1),F(1),F(1),F(1)),'ALLC'
        b=F(1,2)/max(F(1),v,1-u)
        return (F(1),1-b*v,b*(1-u),b),'ZD kappa=1'
    for (a,d) in ALTGRID:
        if altfeas(u,v,a,d): return (a,F(0),F(1),d),'alternator'
    return None,'NONE'

import itertools
N=17
grid=[F(i,2) for i in range(-2*4,2*4+1)]   # -4 .. 4 in steps of 1/2
bad=[];worst_ratio=0;rows=0
for u in grid:
    for v in grid:
        p,tag=partner(u,v)
        if p is None: bad.append((u,v,'no candidate')); continue
        g1=gap(p,u,v,F(1,10**6)); g2=gap(p,u,v,F(1,10**8))
        e1=effdef(p,u,v,F(1,10**6)); e2=effdef(p,u,v,F(1,10**8))
        rows+=1
        # both must be <= 0 and O(eps): the ratio at 1e-6 over 1e-8 should be ~100
        ok = (g1<=0 and g2<=0 and e1<=0 and e2<=0
              and abs(g2)<=abs(g1) and abs(e2)<=abs(e1)
              and abs(float(g2))<1e-6 and abs(float(e2))<1e-6)
        if not ok: bad.append((u,v,tag,float(g1),float(g2),float(e1),float(e2)))
print(f"{rows} games on the 1/2-grid over [-4,4]^2; failures: {len(bad)}")
for b in bad[:20]: print("  ",b)

print()
print("worked examples (exact Fractions, eps = 1e-6 and 1e-8):")
for (u,v) in [(F(-2),F(3)),(F(1,2),F(-3)),(F(1,2),F(2,5)),(F(0),F(9,10)),(F(1,2),F(5)),
              (F(3),F(-3,2)),(F(1),F(0)),(F(1,2),F(1,2)),(F(4),F(-2)),(F(-3),F(5))]:
    p,tag=partner(u,v)
    g1,g2=gap(p,u,v,F(1,10**6)),gap(p,u,v,F(1,10**8))
    e1,e2=effdef(p,u,v,F(1,10**6)),effdef(p,u,v,F(1,10**8))
    side='below' if u+v<1 else ('on' if u+v==1 else 'above')
    print(f"  (u,v)=({u},{v}) {side:5s} {tag:11s} p={tuple(str(x) for x in p)}")
    print(f"        gap/eps: {float(g1*10**6):+.6f} -> {float(g2*10**8):+.6f} ;"
          f"  effdef/eps: {float(e1*10**6):+.6f} -> {float(e2*10**8):+.6f}")
