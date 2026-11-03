"""High-precision finite-eps deficits (convention-B / finite-error analysis)."""
import mpmath as mp, itertools
mp.mp.dps = 60
SWAP=(0,2,1,3)
VERT=[tuple(b) for b in itertools.product([0,1],repeat=4)]

def cl(x,eps): return eps+(1-2*eps)*mp.mpf(x) if not isinstance(x,mp.mpf) else eps+(1-2*eps)*x

def stat(p,q,eps):
    pc=[eps+(1-2*eps)*mp.mpf(x) for x in p]; qc=[eps+(1-2*eps)*mp.mpf(x) for x in q]
    M=mp.zeros(4,4)
    for s in range(4):
        a,b=pc[s],qc[SWAP[s]]
        M[s,0]=a*b; M[s,1]=a*(1-b); M[s,2]=(1-a)*b; M[s,3]=(1-a)*(1-b)
    A=mp.zeros(4,4); rhs=mp.matrix([0,0,0,1])
    for i in range(3):
        for j in range(4): A[i,j]=M[j,i]-(1 if i==j else 0)
    for j in range(4): A[3,j]=1
    return mp.lu_solve(A,rhs)

def selfXYZ(p,eps):
    a,b,c,d=[eps+(1-2*eps)*mp.mpf(x) for x in p]
    beta=(1-b)*(1-c); gam=b*c
    Xn=d*(d*beta+gam*(2-d)); Yn=2*d*(1-a)*(1+a-d); Zn=(1-a)*((1+a)*beta+gam*(1-a))
    S=Xn+Yn+Zn
    return Xn/S,Yn/S,Zn/S

def edef(p,u,v,eps):
    X,Y,Z=selfXYZ(p,eps)
    return X+mp.mpf(1)/2*Y*(1+u+v)-max(mp.mpf(1),(1+u+v)/2)

def rm(p,sign,eps):
    return min(sign*(w[2]-w[1]) for w in (stat(p,q,eps) for q in VERT))

def d_eps(p,u,v,eps):
    sign = 1 if (u-v)<1 else -1
    return min(edef(p,u,v,eps), rm(p,sign,eps))
