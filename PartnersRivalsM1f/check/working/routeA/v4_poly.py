"""V4: (i) the affine/polygon structure checked on the families that actually have a Nash region;
       (ii) D2, the alternators above the switch line; (iii) the exact polygon N(p) by LP."""
import sys, numpy as np
from fractions import Fraction as F
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib

GEN=[f"{i:04b}"[::-1] for i in range(16)]

def Wlim(sig,p,eps=F('1e-9')):
    """limit W(sigma,p) approximated exactly at tiny eps (Fractions, no float)."""
    return m1lib.stationary_exact(sig,p,eps)

def coeffs(p, eps=F(1,10**9)):
    """the 16 triples (B,C,D) with the limit-Nash condition  B(u-1)+C v <= D."""
    Ws=m1lib.stationary_exact(p,p,eps)
    out=[]
    for sig in m1lib.DET_EXACT:
        W=m1lib.stationary_exact(sig,p,eps)
        out.append((W[1]-Ws[1], W[2]-Ws[2], W[3]-Ws[3]))
    return out

def gap_from_coeffs(C,u,v):
    return -max(float(B)*(u-1)+float(Cc)*v-float(D) for (B,Cc,D) in C)

print("=== (i) affine form vs direct gap, exact chain at eps=1e-9, on structured p ===")
fams = {
 'ZD kappa=1 (u<=1,v>=0)': (F(1),F(1)-F(1,5)*F(2),F(1,5)*(1-F(1,2)),F(1,5)),   # u=1/2,v=2,b=1/5
 'alternator (p1,0,1,p4)': (F(3,10),F(0),F(1),F(4,5)),
 'alternator (p1,1,0,p4)': (F(3,10),F(1),F(0),F(4,5)),
 'constant x=1/3':         (F(1,3),)*4,
 'WSLS':                   (F(1),F(0),F(0),F(1)),
 'TFT':                    (F(1),F(0),F(1),F(0)),
 'interior generic':       (F(2,7),F(5,9),F(1,3),F(4,11)),
}
for name,p in fams.items():
    C=coeffs(p)
    worst=0
    for (u,v) in [(0.5,2),(-1,0.3),(3,-2),(0.2,0.9),(4,4),(-3,-3),(1.5,-0.5)]:
        g1=m1lib.nash_gap(np.array([float(x) for x in p]),u,v,1e-9)
        g2=gap_from_coeffs(C,u,v)
        worst=max(worst,abs(g1-g2))
    print(f"  {name:24s} max|direct - affine| = {worst:.3e}")

print()
print("=== (ii) D2: alternators above the switch line ===")
for p in [(F(3,10),F(0),F(1),F(4,5)),(F(1,2),F(0),F(1),F(1,2)),(F(0),F(0),F(1),F(1,2)),
          (F(9,10),F(0),F(1),F(1,100)),(F(1,100),F(1),F(0),F(99,100)),(F(2,5),F(1),F(0),F(3,5))]:
    pf=np.array([float(x) for x in p])
    rows=[]
    for (u,v) in [(0.5,5),(2,2),(0.1,1.5),(-3,6),(6,-3),(1,1),(0.9,0.2),(-1,-1)]:
        e=1e-9
        gap=m1lib.nash_gap(pf,u,v,e); eff=m1lib.efficiency_deficit(pf,u,v,e)
        rows.append((u,v,gap,eff))
    print(f"  p={tuple(str(x) for x in p)}")
    for (u,v,g,ef) in rows:
        side = 'above' if u+v>1 else ('on' if u+v==1 else 'below')
        print(f"     (u,v)=({u:5.1f},{v:5.1f}) {side:5s}  nash_gap={g:+.3e}  eff_deficit={ef:+.3e}"
              f"   {'NASH' if g>-1e-6 else '    '} {'EFF' if ef>-1e-6 else '   '}")
