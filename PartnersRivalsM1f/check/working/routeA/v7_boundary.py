import numpy as np
def feasible(u,v,A,D):
    F1 = (D+1)*(1-u) - (1-D)*v <= 2
    F2 = A*(u-1) + (2-A)*v >= 0
    F3 = (D+1)*(u-1) + (2-A)*(1+v) >= 0
    F4 = (A*D+D+1)*(1-u) + (A*D-D-1)*v <= 2
    F5 = (A+D+1)*(1-u) + (A+D-3)*v <= 2
    F6 = (1+D-A)*(u+v-1) + 2*(1-A) >= 0
    F7 = (A*D-D-1)*(u-1) - (A*D-2*A-D+3)*v <= 2*(1-A)
    return F1&F2&F3&F4&F5&F6&F7
names=['F1(ALLD)','F2(ALLC)','F3(1100)','F4(0001)','F5(0101)','F6(0110/WSLS)','F7(0111)']
def which(u,v,a,d):
    A=np.array([a]);D=np.array([d])
    F=[ (D+1)*(1-u)-(1-D)*v-2, -(A*(u-1)+(2-A)*v), -((D+1)*(u-1)+(2-A)*(1+v)),
        (A*D+D+1)*(1-u)+(A*D-D-1)*v-2, (A+D+1)*(1-u)+(A+D-3)*v-2,
        -((1+D-A)*(u+v-1)+2*(1-A)), (A*D-D-1)*(u-1)-(A*D-2*A-D+3)*v-2*(1-A)]
    return [(n,float(x[0])) for n,x in zip(names,F)]

n=400
aa=np.linspace(0,1,n,endpoint=False); dd=np.linspace(1.0/n,1,n)
A,D=np.meshgrid(aa,dd,indexing='ij')

print("For each u, the smallest v (grid 1/200) with an alternator Nash AND u+v>1:")
for u in [-4,-3,-2,-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,4,6]:
    got=None
    for v in np.arange(max(-6, 1-u+1e-9), 12, 0.005):
        if u+v<=1: continue
        F=feasible(u,v,A,D)
        if F.any():
            got=v
            idx=np.argwhere(F); a0,d0=aa[idx[0][0]],dd[idx[0][1]]
            break
    print(f"  u={u:6.2f}:  vmin={got}   (u+v-1 = {None if got is None else round(u+got-1,4)})  witness (a,d)=({a0:.3f},{d0:.3f})" if got is not None else f"  u={u:6.2f}:  none up to v=12")
print()
print("same but WITHOUT the above-the-line restriction (pure Nash region of the alternator family):")
for u in [-4,-3,-2,-1,0,1,2,4]:
    got=None
    for v in np.arange(-8,12,0.005):
        if feasible(u,v,A,D).any(): got=v; break
    print(f"  u={u:6.2f}: vmin={got}")
