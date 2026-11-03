"""V3 (items 4,5): existence of a limit-Nash strategy at EVERY game, via constant strategies."""
import sys, numpy as np
from fractions import Fraction as F
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib

def const(x): return (F(x),)*4

def gap_exact(p,u,v,eps):
    e=F(eps); u=F(u); v=F(v)
    own = m1lib.stationary_exact(p,p,e)
    own = own[0]+own[1]*u+own[2]*(1+v)
    best = max((lambda w: w[0]+w[1]*u+w[2]*(1+v))(m1lib.stationary_exact(s,p,e))
               for s in m1lib.DET_EXACT)
    return own-best

print("=== (1) mixed stage-NE constant strategy x* = u/(u+v), anti-coordination games u,v>0 ===")
for (u,v) in [(F(1,2),F(1,2)),(F(1,4),3),(2,F(1,10)),(3,7),(F(1,100),F(1,3)),(5,5),(F(1,2),5),(8,F(1,5))]:
    x = F(u)/(F(u)+F(v)); p = const(x)
    row=[]
    for eps in ['1/1000','1/100000','1/10000000']:
        g = gap_exact(p,u,v,eps); row.append(float(g)/float(F(eps)))
    A = (F(u)+F(v))*(1-2*x)     # = v-u
    pred = -float(abs(A))*float(x if A>0 else 1-x)
    print(f"  (u,v)=({u},{v})  x*={x}  gap/eps at eps=1e-3,1e-5,1e-7: "
          f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f}   predicted leading {pred:.6f}   (v-u={float(A)})")

print()
print("=== (2) the pure ones ===")
for (u,v,name,p) in [(-1,2,'ALLD',const(0)),(F(-1,10),5,'ALLD',const(0)),
                     (3,-1,'ALLC',const(1)),(2,F(-1,100),'ALLC',const(1)),
                     (0,3,'ALLD',const(0)),(3,0,'ALLC',const(1))]:
    g=[float(gap_exact(p,u,v,e)) for e in ['1/1000','1/100000']]
    print(f"  (u,v)=({u},{v}) {name}: gap = {g[0]:.3e}, {g[1]:.3e}")

print()
print("=== (3) EXISTENCE SWEEP: is some limit-Nash memory-one strategy present at every game? ===")
rng=np.random.default_rng(7)
bad=[]
for u in np.linspace(-6,6,49):
    for v in np.linspace(-6,6,49):
        if u>0 and v>0: x=u/(u+v)
        elif v<=0:      x=1.0
        else:           x=0.0
        p=np.array([x]*4)
        g=m1lib.nash_gap(p,u,v,1e-8)
        if g < -1e-6: bad.append((u,v,x,g))
print(f"   49x49 grid over [-6,6]^2 with the stage-NE constant strategy: {len(bad)} games with gap < -1e-6")
if bad: print("   worst:", sorted(bad,key=lambda t:t[3])[:5])

print()
print("=== (4) the 16 vertices: where do THEY have no limit-Nash equilibrium? ===")
N=241; us=np.linspace(-4,4,N); vs=np.linspace(-4,4,N)
have=np.zeros((N,N),bool)
for i,u in enumerate(us):
    for j,v in enumerate(vs):
        have[i,j]=any(m1lib.nash_gap(s,u,v,1e-8) > -1e-6 for s in m1lib.DET)
frac=1-have.mean()
print(f"   fraction of the box [-4,4]^2 with NO limit-Nash vertex: {frac*100:.3f} %")
# is that set exactly {u>0,v>0} minus the unit square?
pred=np.zeros((N,N),bool)
for i,u in enumerate(us):
    for j,v in enumerate(vs):
        pred[i,j]= (u>0 and v>0 and not (u<=1 and v<=1))
print(f"   predicted {{u>0,v>0}} \\ [0,1]^2 : {pred.mean()*100:.3f} % ; mismatching cells: {(pred!=(~have)).sum()}")
mm=np.argwhere(pred!=(~have))
print("   sample mismatches (u,v,hasNash):", [(round(us[i],3),round(vs[j],3),bool(have[i,j])) for i,j in mm[:12]])
