import sys
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib as L
from fractions import Fraction as F

P1=[0.0,0.1,0.3,0.5,0.7,0.9,0.99]; P4=[0.01,0.1,0.3,0.5,0.7,0.9,1.0]
for (uu,vv) in [(0.0,1.5),(0.5,1.0),(1.25,0.25)]:
    for corner in [(0,1),(1,0)]:
        bad=[]
        for p1 in P1:
            for p4 in P4:
                if p1==0.0 and p4==1.0: continue
                p=(p1,corner[0],corner[1],p4)
                g=L.nash_gap(p,uu,vv,1e-7)
                if g<=-1e-5: bad.append((p1,p4,g))
        print(f"(u,v)=({uu},{vv}) corner={corner}: {len(bad)} failures ->", bad[:6])
print()
print("EXACT check of the two suspects, corner (0,1), (u,v)=(0,1.5):")
for (p1,p4) in [(0.0,0.9),(0.0,0.7),(0.9,1.0),(0.99,1.0),(0.0,0.5),(0.5,1.0)]:
    p=(F(p1).limit_denominator(100),F(0),F(1),F(p4).limit_denominator(100))
    print(f"  p={tuple(str(x) for x in p)}")
    for k in [3,4,5,6,8,10]:
        e=F(1,10**k)
        own=L.stationary_exact(p,p,e); pi_own=own[0]+own[1]*F(0)+own[2]*F(1+1.5).limit_denominator(10)
        best=max(sum(w[i]*c for i,c in zip(range(3),[F(1),F(0),F(5,2)]))
                 for w in (L.stationary_exact(s,p,e) for s in L.DET_EXACT))
        print(f"    eps=1e-{k}: gap = {float(pi_own-best):+.10f}")
