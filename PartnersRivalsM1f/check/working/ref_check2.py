import numpy as np, math
def selfpay(p0,p,d,u,v):
    sw={0:0,1:2,2:1,3:3}
    M=np.zeros((4,4))
    for s in range(4):
        a=p[s]; b=p[sw[s]]
        M[s]=[a*b,a*(1-b),(1-a)*b,(1-a)*(1-b)]
    v0=np.array([p0*p0,p0*(1-p0),(1-p0)*p0,(1-p0)**2])
    w=(1-d)*v0@np.linalg.inv(np.eye(4)-d*M)
    return w@np.array([1.0,u,1+v,0.0])

for (u,v) in [(2.0,0.5),(0.5,1.5),(1.0,1.0)]:
    E=(1+u+v)/2.0
    print(f"\n(u,v)=({u},{v})  E={E}  sqrt(E(E-1))={math.sqrt(E*(E-1)):.7f}")
    for d in [0.99,0.999,0.9999,0.99999]:
        B=2*E-(1-d)
        a=(B-math.sqrt(B*B-4*d*E))/(2*d)
        h=(1-d)*(a*a*(E-1)+(1-a)**2*E)/(1-d+2*d*a*(1-a))
        # brute force over the alternating family + fine grid on (p0,p1,p4)
        best=-1e9
        for x in np.linspace(0,1,2001):
            for (p2,p3) in [(0.0,1.0),(1.0,0.0)]:
                val=selfpay(x,[x,p2,p3,x],d,u,v); best=max(best,val)
        # unrestricted random search
        rng=np.random.default_rng(0)
        for _ in range(20000):
            z=rng.random(5); best=max(best,selfpay(z[0],list(z[1:]),d,u,v))
        print(f"  d={d}: closed-form deficit={h:.8f}  h/(1-d)={h/(1-d):.6f} | "
              f"numeric best deficit={E-best:.8f} ratio={(E-best)/(1-d):.6f}  a_delta={a:.6f}")
