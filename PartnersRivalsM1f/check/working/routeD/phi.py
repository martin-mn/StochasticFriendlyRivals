import numpy as np, fast, itertools
def dW(P): return np.maximum.reduce([1-P[:,0], P[:,3], np.maximum(P[:,1]+P[:,2]-1,0)/2])
def dN(P): return np.maximum.reduce([P[:,1], 1-P[:,2], P[:,3]])
def dE(P): return np.maximum.reduce([1-P[:,0], P[:,1], 1-P[:,2]])
W={'W':(-0.5,0.5,1,dW),'N':(0.2,1.5,1,dN),'E':(1.5,-0.2,-1,dE)}
eps=1e-9
def run(P,tag):
    for name,(u,v,sgn,df) in W.items():
        e=fast.eff_deficit(P,u,v,eps); r=fast.rival_margin(P,sgn,eps)
        d=np.minimum(e,r); D=df(P)
        print(f'--- {tag} {name}:  max d over all p = {d.max():+.10f}  (argmax {P[d.argmax()]})')
        for delta in (0.01,0.02,0.05,0.10,0.20,0.30,0.50):
            m=D>=delta-1e-12
            if m.sum()==0: continue
            print(f'      Phi({delta:.2f}) = {d[m].max():+.10f}   -delta/(1+delta) = {-delta/(1+delta):+.10f}'
                  f'   {"OK" if d[m].max()<=-delta/(1+delta)+1e-7 else "EXCEEDS"}')
h=0.02; g=np.linspace(0,1,int(round(1/h))+1)
# grid too big for one array; chunk it
P=np.array(list(itertools.product(g,repeat=4)))
run(P,'grid h=0.02')
rng=np.random.default_rng(7)
run(rng.random((400000,4)),'random interior 4e5')
# random on the closure of the relevant faces
Q=rng.random((400000,4)); Q[:,0]=1.0; Q[:,3]=0.0
run(Q,'random p1=1,p4=0')
