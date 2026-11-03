import sys, numpy as np
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib
rng=np.random.default_rng(42)

print("=== (A) is the limit-Nash set ever full-dimensional (4-dim) in the cube? ===")
for (u,v) in [(-1,-1),(-2,-3),(0.5,0.5),(-1,2),(2,-1),(0.2,0.2),(-4,-4),(3,3)]:
    hits=0; N=4000
    for _ in range(N):
        p=rng.random(4)
        if m1lib.nash_gap(p,u,v,1e-9) > -1e-6: hits+=1
    print(f"  (u,v)=({u},{v}): {hits}/{N} uniformly random interior p are limit-Nash")

print()
print("=== (B) curvature probe: 2-parameter sub-families, is the boundary of the union a line? ===")
def nash_p1one(u,v,m,q3,q4):
    t=u-1.0; s=q3+q4-q3*q4
    return ((q4*v<=m)&(m*t+q3*v<=0)&(q4*(m*t+v)<=m*(1-q3))&(m*t+s*v<=m*(1-q3))
            &(m*(1-q4)*t+s*v<=m)&(m*(1-q4)*t+q3*v<=q3*m))
n=600
M,Q3=np.meshgrid(np.linspace(1e-6,1,n),np.linspace(0,1,n),indexing='ij')
for q4 in [0.25,0.5,0.75,1.0]:
    print(f"  slice p1=1, q4={q4}: for each u the largest v still Nash for some (q2,q3)")
    out=[]
    for u in [-3,-2,-1,0,0.5,0.9,1.0,1.5,2,3]:
        lo,hi=-20.,60.
        for _ in range(60):
            mid=(lo+hi)/2
            if nash_p1one(u,mid,M,Q3,q4).any(): lo=mid
            else: hi=mid
        out.append((u,lo))
    print("    ",[(u,round(v,5)) for u,v in out])
    # second differences of vmax(u) -- zero on a straight segment
    xs=np.array([o[0] for o in out]); ys=np.array([o[1] for o in out])
    print("     vmax(u) values:",np.round(ys,5))
