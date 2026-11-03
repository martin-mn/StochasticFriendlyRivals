import numpy as np, sys
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib
def feas01(u,v,A,D):
    return (((D+1)*(1-u)-(1-D)*v<=2) & (A*(u-1)+(2-A)*v>=0) & ((D+1)*(u-1)+(2-A)*(1+v)>=0) &
            ((A*D+D+1)*(1-u)+(A*D-D-1)*v<=2) & ((A+D+1)*(1-u)+(A+D-3)*v<=2) &
            ((1+D-A)*(u+v-1)+2*(1-A)>=0) & ((A*D-D-1)*(u-1)-(A*D-2*A-D+3)*v<=2*(1-A)))
rng=np.random.default_rng(5); NS=40000
A=rng.random(NS); D=rng.random(NS)
G=161; us=np.linspace(-4,4,G); vs=np.linspace(-4,4,G)
zeros=[]
for i,u in enumerate(us):
    for j,v in enumerate(vs):
        if u+v>1 and not feas01(u,v,A,D).any(): zeros.append((u,v))
print(len(zeros),"games strictly above the line with MC area 0; u+v-1 values:",
      sorted(set(round(u+v-1,4) for u,v in zeros)))
print("sample:",zeros[:10])
# fine deterministic search on a log grid clustered at the (a,d)=(1,0) corner
ee=np.array([1-10.0**-k for k in np.arange(0.2,9,0.2)]); dd=np.array([10.0**-k for k in np.arange(0.2,9,0.2)])
AA,DD=np.meshgrid(np.concatenate([ee,np.linspace(0,1,200)]),np.concatenate([dd,np.linspace(1e-9,1,200)]),indexing='ij')
still=[]
for (u,v) in zeros:
    F=feas01(u,v,AA,DD)
    if F.any():
        k=np.argwhere(F)[0]; a0,d0=AA[k[0],k[1]],DD[k[0],k[1]]
        g=m1lib.nash_gap(np.array([a0,0.,1.,d0]),u,v,1e-11)
        ef=m1lib.efficiency_deficit(np.array([a0,0.,1.,d0]),u,v,1e-11)
        still.append(None)
        print(f"  ({u:.3f},{v:.3f}) u+v-1={u+v-1:.4f}: witness a={a0:.10f} d={d0:.3e}"
              f"  chain nash_gap={g:+.2e} eff_def={ef:+.2e}")
    else:
        still.append((u,v))
print("still no witness:",[z for z in still if z])
