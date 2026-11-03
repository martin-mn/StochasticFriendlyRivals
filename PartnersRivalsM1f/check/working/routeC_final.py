"""Route C final checks: criterion vs closed form on a dense grid; closed form vs the
float pair chain by Richardson extrapolation on a random sample; witnesses; D2 revisited."""
import sys, random, numpy as np
from fractions import Fraction as F
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib as L
from routeC_verify import piA, piB, GEN, pred_gap

random.seed(7); np.random.seed(7)

# ---- 1. reduced criterion vs full 16-inequality predicted gap, dense grid, exact rationals
bad=0; n=0
for corner in [(0,1),(1,0)]:
    for ai in range(0,20):
        a=F(ai,20)
        for di in range(1,21):
            d=F(di,20)
            if a==0 and d==1: continue
            for ui in range(-12,13):
                u=F(ui,4)
                for vi in range(-12,13):
                    v=F(vi,4)
                    if u+v<=1: continue
                    n+=1
                    X=u-v-1; Y=u+v-1
                    crit = ((d*X+Y+2>=0) and ((a-1)*X+Y>=0)) if corner==(0,1) else (X==0)
                    if (pred_gap(corner,a,d,u,v)>=0)!=crit:
                        bad+=1
print(f"[1] criterion vs 16 inequalities: {n} cases, {bad} mismatches")

# ---- 2. closed-form limit vs float pair chain (Richardson in eps), random sample
def float_gap(corner,a,d,u,v,eps):
    p=(a,corner[0],corner[1],d)
    return L.nash_gap(p,u,v,eps)
worst=0.0; N=3000
for _ in range(N):
    corner=random.choice([(0,1),(1,0)])
    a=random.random()*0.999; d=0.001+random.random()*0.999
    while True:
        u=random.uniform(-4,6); v=random.uniform(-4,6)
        if u+v>1.05: break
    pg=float(pred_gap(corner,F(a).limit_denominator(10**9),F(d).limit_denominator(10**9),
                      F(u).limit_denominator(10**9),F(v).limit_denominator(10**9)))
    g7=float_gap(corner,a,d,u,v,1e-7); g8=float_gap(corner,a,d,u,v,1e-8)
    rich=(10*g8-g7)/9
    worst=max(worst,abs(rich-pg))
print(f"[2] float Richardson vs closed form, {N} random points: max |diff| = {worst:.3e}")

# ---- 3. corner B is beaten by ALLD (T>S) or ALLC (T<S): explicit
print("[3] corner B witnesses")
for (u,v) in [(0.5,5.0),(3.0,-0.5),(0.0,2.0),(5.0,-3.0)]:
    p=(0.4,1,0,0.6); Emax=max(1.0,(1+u+v)/2)
    wD=L.stationary(p,p,1e-9)
    piD=L.payoff(L.stationary((0,0,0,0),p,1e-9),u,v)
    piC=L.payoff(L.stationary((1,1,1,1),p,1e-9),u,v)
    print(f"  (u,v)=({u},{v}) Emax={Emax:.4f}  pi(p,p)={L.payoff(wD,u,v):.6f}"
          f"  pi(ALLD,p)={piD:.6f} (=T={1+v})  pi(ALLC,p)={piC:.6f} (=S={u})")

# ---- 4. D2 revisited: its four points, and the counterexamples to the general claim
print("[4] D2 revisited at (u,v)=(0.5,5): corner A all Nash (u>=0,v>=0)")
for (p1,p4) in [(0.42436,0.42436),(0.3,0.8),(0.7,0.2),(0.5,0.5),(0.0,0.001),(0.99,1.0)]:
    print(f"   A p=({p1},0,1,{p4}) gap(1e-8)={L.nash_gap((p1,0,1,p4),0.5,5.0,1e-8):+.3e}")
for (p1,p4) in [(0.42436,0.42436),(0.3,0.8)]:
    print(f"   B p=({p1},1,0,{p4}) gap(1e-8)={L.nash_gap((p1,1,0,p4),0.5,5.0,1e-8):+.3e}")
print("   counterexamples inside corner A:")
for (u,v,p1,p4,why) in [(-3.0,4.5,0.5,0.9,'u<0, d>d*=5/17'),(-3.0,4.5,0.5,0.2,'u<0, d<d*'),
                        (5.0,-3.0,0.5,0.5,'v<0, a<a*=6/7'),(5.0,-3.0,0.95,0.5,'v<0, a>a*')]:
    print(f"   A p=({p1},0,1,{p4}) at (u,v)=({u},{v}): gap(1e-8)={L.nash_gap((p1,0,1,p4),u,v,1e-8):+.6f}   [{why}]")

# ---- 5. area of the partner slice in (p1,p4)
def area(u,v):
    X=u-v-1; Y=u+v-1
    assert Y>0
    ad=1.0 if X<=0 else max(0.0,min(1.0,Y/X))       # a-constraint
    dd=1.0 if X>=0 else max(0.0,min(1.0,(Y+2)/(-X)))# d-constraint
    return ad*dd
print("[5] area of the corner-A partner slice, Monte-Carlo vs formula")
for (u,v) in [(0.5,5.0),(-3.0,4.5),(5.0,-3.0),(-1.0,3.0),(4.0,-1.0),(2.0,1.0),(0.0,4.0)]:
    m=0; M=40000
    for _ in range(M):
        a=random.random(); d=random.random()
        X=u-v-1; Y=u+v-1
        if d*X+Y+2>=0 and (a-1)*X+Y>=0: m+=1
    print(f"  (u,v)=({u},{v}): MC={m/M:.4f}  formula={area(u,v):.4f}")
