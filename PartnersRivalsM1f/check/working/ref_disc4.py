import sys, numpy as np, itertools
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib as L
SWAP=(0,2,1,3)
def occf(p0,p,q0,q,d):
    p=np.asarray(p,float); q=np.asarray(q,float)
    M=np.empty((4,4))
    for s in range(4):
        a=p[s]; b=q[SWAP[s]]
        M[s]=[a*b,a*(1-b),(1-a)*b,(1-a)*(1-b)]
    v0=np.array([p0*q0,p0*(1-q0),(1-p0)*q0,(1-p0)*(1-q0)])
    return (1-d)*np.linalg.solve((np.eye(4)-d*M).T, v0)
u,v=1.5,0.5   # T = S  (u-v=1),  u+v=2 > 1
rng=np.random.default_rng(1)
worst=0
for _ in range(4000):
    p0=rng.random(); p=rng.random(4); q0=rng.random(); q=rng.random(4); d=rng.random()*0.999
    w=occf(p0,p,q0,q,d)
    m=(w[2]-w[1])*(1+v-u)
    worst=max(worst,abs(m))
print("T=S ray: max |pi(p,q)-pi(q,p)| over 4000 random (p,q,delta) =",worst)
# the V*-optimal strategy there is therefore unbeatable trivially
E=(1+u+v)/2; 
for d in [0.5,0.9,0.99]:
    a=(2*E-(1-d)-np.sqrt((2*E-(1-d))**2-4*d*E))/(2*d)
    p=[a,0.0,1.0,a]
    w=occf(a,p,a,p,d); print(f" d={d} a*={a:.4f} selfpay={L.payoff(w,u,v):.6f} Emax={E}")
# project's undiscounted predicates on the ray
for p in [(0.5,0,1,0.5),(0.3,0,1,0.7)]:
    print(" undiscounted:",p,"in_eff_alt",L.in_eff_alt(np.array(p,float),u,v) if hasattr(L,'in_eff_alt') else '?')
