"""V10: above the switch line the only efficient strategies are the alternators
Eff_alt = {(p1,0,1,p4)} u {(p1,1,0,p4)}; measure the area of the (a,d) square on which each is
limit-Nash, as a function of the game."""
import numpy as np, sys
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib

def feas01(u,v,A,D):      # p = (a,0,1,d)
    return (((D+1)*(1-u)-(1-D)*v<=2) & (A*(u-1)+(2-A)*v>=0) & ((D+1)*(u-1)+(2-A)*(1+v)>=0) &
            ((A*D+D+1)*(1-u)+(A*D-D-1)*v<=2) & ((A+D+1)*(1-u)+(A+D-3)*v<=2) &
            ((1+D-A)*(u+v-1)+2*(1-A)>=0) & ((A*D-D-1)*(u-1)-(A*D-2*A-D+3)*v<=2*(1-A)))

rng=np.random.default_rng(5)
NS=40000
A=rng.random(NS); D=rng.random(NS)          # (a,d) uniform on the unit square
G=161; us=np.linspace(-4,4,G); vs=np.linspace(-4,4,G)
ar=np.zeros((G,G))
for i,u in enumerate(us):
    for j,v in enumerate(vs):
        ar[i,j]=feas01(u,v,A,D).mean()
above=np.array([[us[i]+vs[j]>1 for j in range(G)] for i in range(G)])
print("strictly above the line: min area =",ar[above].min(),"  #games with area 0:",(ar[above]==0).sum())
print("below/on: max area =",ar[~above].max())
ui=[np.argmin(abs(us-x)) for x in [-3,-1,0,0.5,1,1.5,2,3]]
vj=[np.argmin(abs(vs-x)) for x in [-3,-1,0,0.5,1,2,3,4]]
print("area of the alternator (a,d) square that is limit-Nash:")
print("        "+"".join(f"v={vs[j]:5.2f} " for j in vj))
for i in ui:
    print(f"u={us[i]:5.2f} "+"".join(f"{ar[i,j]:7.4f} " for j in vj))

print()
print("the other alternating corner p=(a,1,0,d): any limit-Nash anywhere?")
bad=[]
for _ in range(4000):
    a,d=rng.random(2); u,v=rng.uniform(-6,6,2)
    if m1lib.nash_gap(np.array([a,1.,0.,d]),u,v,1e-9) > -1e-6: bad.append((a,d,u,v))
print("  hits over 4000 random (a,d,u,v):",len(bad), bad[:5])
