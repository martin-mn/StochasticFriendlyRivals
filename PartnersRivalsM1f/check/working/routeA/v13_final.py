import sys, numpy as np, itertools
from fractions import Fraction as F
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib
from v1_facets import f_sig, f_self, GEN, facets

print("=== (A) how many DISTINCT facet rows does an interior p give? ===")
rng=np.random.default_rng(1)
cnt=set()
for _ in range(200):
    pv=rng.uniform(0.05,0.95,4)
    F16,_=facets(pv)
    rows=[tuple(np.round(x,12)) for x in F16.values()]
    groups={}
    for g,r in zip(GEN,rows): groups.setdefault(r,[]).append(g)
    cnt.add(tuple(sorted(tuple(sorted(v)) for v in groups.values())))
print("  distinct coincidence patterns over 200 interior p:",len(cnt))
for c in cnt: print("   ",[x for x in c if len(x)>1], " -> #distinct rows =",len(c))

print()
print("=== (B) worked polygon: p = (1,1/2,1/2,1/2) ===")
p=np.array([1.,.5,.5,.5])
pred=lambda u,v: (u+v<=1+1e-12) and (v<=1+1e-12)
bad=0
for u in np.linspace(-5,5,101):
    for v in np.linspace(-5,5,101):
        g=m1lib.nash_gap(p,u,v,1e-9)
        if (g>-1e-6)!=pred(u,v): bad+=1
print("  claim: limit-Nash exactly on {u+v<=1, v<=1}.  mismatches on a 101x101 grid:",bad)
for (u,v) in [(F(0),F(0)),(F(-3),F(1)),(F(0),F(1)),(F(1,2),F(1,2)),(F(0),F(11,10)),(F(3,2),F(0))]:
    g=m1lib.nash_gap(np.array([1.,.5,.5,.5]),float(u),float(v),1e-9)
    print(f"    (u,v)=({u},{v}) in region={pred(float(u),float(v))}  nash_gap={g:+.3e}")

print()
print("=== (C) symmetry receipt: the self-play chain is invariant under p2 <-> p3 ===")
worst=0
for _ in range(300):
    pv=rng.uniform(0,1,4); ps=np.array([pv[0],pv[2],pv[1],pv[3]])
    w1=m1lib.self_play(pv,1e-7); w2=m1lib.self_play(ps,1e-7)
    worst=max(worst,abs(w1[0]-w2[0]),abs(w1[3]-w2[3]),abs(w1[1]-w2[2]),abs(w1[2]-w2[1]))
print("  max |W(p,p) - swap W(p',p')| over 300 random p, eps=1e-7:",worst)
w=m1lib.stationary_exact((F(1,3),F(1,5),F(2,7),F(3,11)),(F(1,3),F(1,5),F(2,7),F(3,11)),F(1,10**7))
print("  exact self-play W_CD - W_DC (must be exactly 0):", w[1]-w[2])

print()
print("=== (D) ZD equalisers: feasibility region, checked against a direct search ===")
def zd_exists(u,v):
    return (u<=1) and (v>=-1) and (v>=u-1)
bad=[]
for u in np.linspace(-3,3,61):
    for v in np.linspace(-3,3,61):
        lo,hi=max(0.,u),min(1.,1+v)
        if zd_exists(u,v) != (lo<=hi+1e-15): bad.append((u,v))
print("  {max(0,u) <= kappa <= min(1,1+v) nonempty}  ==  {u<=1, v>=-1, v>=u-1}:",len(bad),"mismatches")
print("  games NOT admitting a ZD equaliser but admitting a limit-Nash strategy, examples:")
for (u,v) in [(2,0.5),(0.5,-2),(3,-3),(1.5,0.2),(2,-0.5)]:
    x = u/(u+v) if (u>0 and v>0) else (1.0 if v<=0 else 0.0)
    p=np.array([x]*4); g=m1lib.nash_gap(p,u,v,1e-9)
    print(f"    (u,v)=({u},{v}): ZD? {zd_exists(u,v)}  ALLD? {u<=0}  constant x={x:.4f} nash_gap={g:+.3e}")
