import numpy as np, itertools
from scipy.optimize import minimize
SWAP=(0,2,1,3)
def occf(p0,p,q0,q,d):
    p=np.asarray(p,float); q=np.asarray(q,float)
    M=np.empty((4,4))
    for s in range(4):
        a=p[s]; b=q[SWAP[s]]
        M[s]=[a*b,a*(1-b),(1-a)*b,(1-a)*(1-b)]
    v0=np.array([p0*q0,p0*(1-q0),(1-p0)*q0,(1-p0)*(1-q0)])
    return (1-d)*np.linalg.solve((np.eye(4)-d*M).T, v0)
def pay(w,u,v): return w[0]+w[1]*u+w[2]*(1+v)

def selfpay(x,u,v,d):
    p0=x[0]; p=x[1:]
    return pay(occf(p0,p,p0,p,d),u,v)

def Vstar_closed(u,v,d):
    sig=1+u+v; E=sig/2
    # min over a of h(a)/(1-d+2 d a(1-a)) * (1-d)
    a=np.linspace(0,1,2000001)
    h=sig/2-sig*a+(sig-1)*a**2
    D=(1-d)*h/(1-d+2*d*a*(1-a))
    i=np.argmin(D)
    return E-D[i], a[i]

for (u,v) in [(2.0,0.5),(0.5,1.5),(1.0,1.0),(1.5,0.5)]:
    sig=1+u+v; E=sig/2
    for d in [0.0,0.5,0.9,0.99,0.999]:
        Vc,astar=Vstar_closed(u,v,d)
        # unrestricted numeric optimum over (p0,p1,p2,p3,p4)
        best=-9
        rng=np.random.default_rng(0)
        for t in range(60):
            x0=rng.random(5)
            r=minimize(lambda x:-selfpay(np.clip(x,0,1),u,v,d),x0,bounds=[(0,1)]*5,method="L-BFGS-B")
            best=max(best,-r.fun)
        print(f"u={u} v={v} d={d}: Emax={E:.6f} closed V*={Vc:.8f} a*={astar:.4f} numeric max={best:.8f} deficit/(1-d)={(E-best)/(1-d) if d<1 else 0:.6f}  sqrt(E(E-1))={np.sqrt(E*(E-1)):.6f}")
    print()
