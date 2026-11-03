from fractions import Fraction as F
from itertools import product
import random

def occ(p0,p,q0,q,d):
    """exact discounted occupancy w = (1-d) v0 (I - d M)^-1 ; states CC,CD,DC,DD (focal first)"""
    # focal prob of C in state s: p[s]; co-player prob of C in state s (its own view):
    # co-player's state for focal-state CC->CC, CD->DC, DC->CD, DD->DD
    sw={0:0,1:2,2:1,3:3}
    M=[[F(0)]*4 for _ in range(4)]
    for s in range(4):
        a=p[s]; b=q[sw[s]]
        M[s][0]=a*b; M[s][1]=a*(1-b); M[s][2]=(1-a)*b; M[s][3]=(1-a)*(1-b)
    v0=[p0*q0, p0*(1-q0), (1-p0)*q0, (1-p0)*(1-q0)]
    # solve w (I - d M) = (1-d) v0  ->  A^T w^T = ...
    A=[[ (F(1) if i==j else F(0)) - d*M[i][j] for j in range(4)] for i in range(4)]
    # unknown row vector w: sum_i w_i A[i][j] = (1-d) v0[j]
    B=[[A[i][j] for i in range(4)] for j in range(4)]   # B[j][i]
    rhs=[(1-d)*v0[j] for j in range(4)]
    n=4
    for c in range(n):
        piv=next(r for r in range(c,n) if B[r][c]!=0)
        B[c],B[piv]=B[piv],B[c]; rhs[c],rhs[piv]=rhs[piv],rhs[c]
        pv=B[c][c]
        B[c]=[x/pv for x in B[c]]; rhs[c]=rhs[c]/pv
        for r in range(n):
            if r!=c and B[r][c]!=0:
                f=B[r][c]
                B[r]=[B[r][k]-f*B[c][k] for k in range(n)]; rhs[r]=rhs[r]-f*rhs[c]
    return rhs

def payoffdiff(p0,p,q0,q,d,u,v):
    w=occ(p0,p,q0,q,d)
    T=1+v; S=u
    return (T-S)*(w[2]-w[1])

def selfpay(p0,p,d,u,v):
    w=occ(p0,p,p0,p,d)
    return w[0]*F(1)+w[1]*u+w[2]*(1+v)+w[3]*F(0)

def Emax(u,v):
    return max(F(1),(u+1+v)/2)

DET=[tuple(F(b) for b in bits) for bits in product([0,1],repeat=4)]
COP=[(F(q0),q) for q in DET for q0 in (0,1)]

def is_rival(p0,p,d,u,v,extra=()):
    for q0,q in list(COP)+list(extra):
        if payoffdiff(p0,p,q0,q,d,u,v)<0: return False
    return True

def is_eff(p0,p,d,u,v):
    return selfpay(p0,p,d,u,v)==Emax(u,v)

# ---------- TEST 1: wedge S characterisation, exact, several games & deltas ----------
print("TEST 1  wedge-S FR set  vs  {p0=p1=1, p2+p3 >= 2-1/d}")
games=[(F(2,5),F(-1)),(F(7,5),F(-9,10)),(F(3,2),F(-1,2)),(F(1),F(0)),(F(1,2),F(-1,2)),(F(0),F(-2))]
deltas=[F(1,10),F(1,2),F(3,5),F(3,4),F(9,10),F(99,100)]
grid=[F(k,4) for k in range(5)]
bad=0; tot=0
for (u,v) in games:
    for d in deltas:
        for p0 in [F(0),F(1,2),F(1)]:
            for p in product(grid,repeat=4):
                tot+=1
                fr = is_eff(p0,list(p),d,u,v) and is_rival(p0,list(p),d,u,v)
                pred = (u+v<=1) and p0==1 and p[0]==1 and ( (u-v>1 and p[1]+p[2]>=2-F(1)/d) or (u-v==1) )
                if fr!=pred:
                    bad+=1
                    if bad<8: print("  MISMATCH",(u,v),d,p0,p,"meas",fr,"pred",pred)
print("  cases",tot,"mismatches",bad)

# ---------- TEST 2: emptiness off wedge S ----------
print("\nTEST 2  emptiness on T>S and on u+v>1")
off=[(F(1,5),F(1,5),'W'),(F(-1,2),F(1,2),'W'),(F(3,10),F(7,10),'W/T=S-below? u-v=-0.4'),
     (F(2),F(1,2),'E'),(F(1,2),F(3,2),'N'),(F(1),F(1),'N/E bdry u-v=0'),
     (F(3,2),F(1,2),'T=S upper ray'),(F(3,2),F(1,2),'dup')]
grid=[F(k,4) for k in range(5)]
for (u,v,lab) in off:
    for d in [F(0),F(1,10),F(1,2),F(9,10),F(99,100)]:
        hits=0; effs=0
        for p0 in [F(0),F(1,4),F(1,2),F(3,4),F(1)]:
            for p in product(grid,repeat=4):
                e=is_eff(p0,list(p),d,u,v)
                if e:
                    effs+=1
                    if is_rival(p0,list(p),d,u,v): hits+=1
        print(f"  (u,v)=({u},{v}) {lab:12s} d={d}  efficient={effs:5d}  friendly_rivals={hits}")

# ---------- TEST 3: delta = 0 edge case ----------
print("\nTEST 3  delta = 0 : is efficiency really {p0=1,p1=1}?")
for (u,v) in [(F(2,5),F(-1)),(F(1,5),F(1,5))]:
    for p in [(F(0),F(0),F(0),F(0)),(F(1),F(1),F(1),F(1)),(F(1,2),F(0),F(1),F(0))]:
        print("  (u,v)",(u,v),"p0=1 p=",p,"selfpay",selfpay(F(1),list(p),F(0),u,v),"Emax",Emax(u,v),
              "eff",is_eff(F(1),list(p),F(0),u,v),"rival",is_rival(F(1),list(p),F(0),u,v))

# ---------- TEST 4: ALLC on the whole closed wedge S, many delta ----------
print("\nTEST 4  ALLC friendly rival across closed wedge S")
ok=True
for iu in range(0,41):
    for iv in range(-40,41):
        u=F(iu,10); v=F(iv,10)
        if not (u+v<=1 and u-v>=1): continue
        for d in [F(0),F(1,3),F(9,10),F(999,1000)]:
            if not (is_eff(F(1),[F(1)]*4,d,u,v) and is_rival(F(1),[F(1)]*4,d,u,v)):
                print("  ALLC FAILS at",u,v,d); ok=False
print("  ALLC ok everywhere on closed wedge S grid:",ok)
