import itertools, numpy as np
from fractions import Fraction as F
SWAP=(0,2,1,3)

def occ(p0,p,q0,q,d):
    # exact discounted occupancy, Fractions
    M=[[0]*4 for _ in range(4)]
    for s in range(4):
        a=p[s]; b=q[SWAP[s]]
        M[s]=[a*b,a*(1-b),(1-a)*b,(1-a)*(1-b)]
    v0=[p0*q0,p0*(1-q0),(1-p0)*q0,(1-p0)*(1-q0)]
    # solve w (I - d M) = (1-d) v0   ->  A^T x = b with A = I - d M
    A=[[ (1 if i==j else 0) - d*M[i][j] for j in range(4)] for i in range(4)]
    # w A = (1-d) v0  => A^T w^T = (1-d) v0^T
    B=[[A[j][i] for j in range(4)] for i in range(4)]
    b=[(1-d)*x for x in v0]
    # gaussian elim exact
    n=4
    Mat=[row[:]+[b[i]] for i,row in enumerate(B)]
    for c in range(n):
        piv=next(r for r in range(c,n) if Mat[r][c]!=0)
        Mat[c],Mat[piv]=Mat[piv],Mat[c]
        pv=Mat[c][c]
        Mat[c]=[x/pv for x in Mat[c]]
        for r in range(n):
            if r!=c and Mat[r][c]!=0:
                f=Mat[r][c]
                Mat[r]=[Mat[r][k]-f*Mat[c][k] for k in range(n+1)]
    return [Mat[i][n] for i in range(n)]

def pay(w,u,v): return w[0]+w[1]*u+w[2]*(1+v)
def payco(w,u,v): return w[0]+w[1]*(1+v)+w[2]*u

# ---- 1. Akin identity check, random rationals
import random
random.seed(7)
def rr(): return F(random.randint(0,9),9)
bad=0
for _ in range(200):
    p0=rr(); p=[rr() for _ in range(4)]; q0=rr(); q=[rr() for _ in range(4)]
    d=F(random.randint(1,19),20)
    w=occ(p0,p,q0,q,d)
    assert sum(w)==1
    c=[1,1,0,0]
    lhs=sum(w[i]*(c[i]-d*p[i]) for i in range(4))
    if lhs != (1-d)*p0: bad+=1
print("Akin identity failures:",bad,"/200")

# ---- 2. wedge-S rival threshold: brute force over ALL deterministic co-players (32) incl opening
def min_margin(p0,p,d,u,v,extra_random=0):
    m=None
    cops=[(q0,q) for q0 in (F(0),F(1)) for q in itertools.product([F(0),F(1)],repeat=4)]
    for q0,q in cops:
        w=occ(p0,list(p),q0,list(q),d)
        val=pay(w,u,v)-payco(w,u,v)
        m=val if m is None else min(m,val)
    return m

u,v=F(7,5),F(-9,10)   # wedge S: u-v=2.3>1, u+v=0.5<1
print("game u,v",u,v,"T-S",1+v-u,"2R-(T+S)",1-u-v)
rows=[]
for d in [F(1,4),F(1,2),F(3,5),F(9,10),F(99,100)]:
    thr=2-1/d
    err=0; tot=0
    for _ in range(300):
        p0=F(1); p=[F(1)]+[rr() for _ in range(3)]
        m=min_margin(p0,p,d,u,v)
        pred = (p[1]+p[2] >= thr)
        got = (m>=0)
        tot+=1
        if pred!=got: err+=1; print("  MISMATCH",d,p,m,float(p[1]+p[2]),float(thr))
    print("delta",d,"thr",float(thr),"mismatches",err,"/",tot)
