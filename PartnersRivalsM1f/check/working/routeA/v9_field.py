"""V9: the partner field below the switch line.  p1 = 1 is forced (Eff_CC), so the partner set
lives in the (q2,q3,q4) cube; measure its volume as a function of the game."""
import numpy as np, sys
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib

def nash_p1one(u,v,m,q3,q4):
    """m = 1-q2.  The six non-trivial facets of the slice p1=1 (derived symbolically in v8)."""
    t=u-1.0; s=q3+q4-q3*q4
    G1 = q4*v <= m
    G2 = m*t + q3*v <= 0
    G3 = q4*(m*t+v) <= m*(1-q3)
    G4 = m*t + s*v <= m*(1-q3)
    G5 = m*(1-q4)*t + s*v <= m
    G6 = m*(1-q4)*t + q3*v <= q3*m
    return G1&G2&G3&G4&G5&G6

# --- validate the closed form against the chain
rng=np.random.default_rng(11)
bad=0
for _ in range(600):
    q2,q3,q4=rng.uniform(0.02,0.98,3); u,v=rng.uniform(-5,5,2)
    g=m1lib.nash_gap(np.array([1.,q2,q3,q4]),u,v,1e-9)
    f=bool(nash_p1one(u,v,1-q2,q3,q4))
    if (g>-1e-6)!=f: bad+=1
print("p1=1 closed form vs chain, 600 random (q,u,v):",bad,"disagreements")

# --- the volume field
NS=20000
Q=rng.random((3,NS))          # q2,q3,q4 uniform on the cube, common random numbers
M=1-Q[0]; Q3=Q[1]; Q4=Q[2]
G=161; us=np.linspace(-4,4,G); vs=np.linspace(-4,4,G)
vol=np.zeros((G,G))
for i,u in enumerate(us):
    for j,v in enumerate(vs):
        vol[i,j]=nash_p1one(u,v,M,Q3,Q4).mean()
np.savez('/private/tmp/claude-501/-Users-martin-Documents-CL1-PartnersRivalsM1f/c1127d6a-179d-4b6f-be64-0786413868ce/scratchpad/A/vol.npz',vol=vol,us=us,vs=vs)
below=np.array([[us[i]+vs[j]<1 for j in range(G)] for i in range(G)])
print("below the line: min volume =",vol[below].min()," fraction of below-line games with volume 0:",
      (vol[below]==0).mean())
print("on/above the line: max volume =",vol[~below].max())
# where is the volume exactly 1 (whole cube)?
print("volume==1 set: fraction",(vol==1).mean(), "  matches {u<=1,v<=0}? ",
      np.array_equal(vol==1, np.array([[ (us[i]<=1 and vs[j]<=0) for j in range(G)] for i in range(G)])))
# print a coarse map of the volume, below the line
print()
print("volume of the partner set (fraction of the (q2,q3,q4) cube), rows u, cols v:")
ui=[np.argmin(abs(us-x)) for x in [-3,-2,-1,0,0.5,0.9,1.0,1.1,1.5,2,3]]
vj=[np.argmin(abs(vs-x)) for x in [-3,-2,-1,-0.5,0,0.25,0.5,1,2,3]]
print("        "+"".join(f"v={vs[j]:5.2f} " for j in vj))
for i in ui:
    print(f"u={us[i]:5.2f} "+"".join(f"{vol[i,j]:7.4f} " for j in vj))
