"""V6: the union over the alternator family of the Nash polygons -- is its boundary a line
arrangement or an algebraic curve?"""
import numpy as np, sys
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib

def feasible(u,v,A,D):
    """vectorised over arrays A,D of (a,d): all 7 constraints of the alternator (a,0,1,d)."""
    F1 = (D+1)*(1-u) - (1-D)*v <= 2
    F2 = A*(u-1) + (2-A)*v >= 0
    F3 = (D+1)*(u-1) + (2-A)*(1+v) >= 0
    F4 = (A*D+D+1)*(1-u) + (A*D-D-1)*v <= 2
    F5 = (A+D+1)*(1-u) + (A+D-3)*v <= 2
    F6 = (1+D-A)*(u+v-1) + 2*(1-A) >= 0
    F7 = (A*D-D-1)*(u-1) - (A*D-2*A-D+3)*v <= 2*(1-A)
    return F1&F2&F3&F4&F5&F6&F7

n=200
aa=np.linspace(0,1,n,endpoint=False); dd=np.linspace(1.0/n,1,n)
A,D=np.meshgrid(aa,dd,indexing='ij')

# sanity: the closed form vs the actual chain, for a handful of (a,d,u,v)
rng=np.random.default_rng(3)
worst=0
for _ in range(400):
    a=rng.uniform(0.02,0.98); d=rng.uniform(0.02,1.0); u=rng.uniform(-5,5); v=rng.uniform(-5,5)
    g=m1lib.nash_gap(np.array([a,0.,1.,d]),u,v,1e-9)
    f=feasible(u,v,np.array([a]),np.array([d]))[0]
    if (g>-1e-6) != bool(f): worst+=1
print("closed-form vs chain disagreements over 400 random (a,d,u,v):",worst)

M=241; us=np.linspace(-4,4,M); vs=np.linspace(-4,4,M)
ok=np.zeros((M,M),bool)
for i,u in enumerate(us):
    for j,v in enumerate(vs):
        ok[i,j]=feasible(u,v,A,D).any()
print("fraction of [-4,4]^2 covered by SOME alternator being limit-Nash:",ok.mean())

# compare with the candidate closed region {u>=0, u+v>=1}
cand=np.zeros((M,M),bool)
for i,u in enumerate(us):
    for j,v in enumerate(vs):
        cand[i,j]= (u>=0) and (u+v>=1)
print("candidate {u>=0,u+v>=1}:",cand.mean(),"mismatch cells:",(ok!=cand).sum())
mm=np.argwhere(ok!=cand)
print("sample mismatches (u,v,ok,cand):",[(round(us[i],3),round(vs[j],3),bool(ok[i,j]),bool(cand[i,j])) for i,j in mm[:20]])

# where exactly does the alternator region live relative to the switch line?
above=np.array([[us[i]+vs[j]>1 for j in range(M)] for i in range(M)])
print("ok and below-or-on line:",(ok&~above).sum(), " ok and above:",(ok&above).sum())
