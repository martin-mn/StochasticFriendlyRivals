import numpy as np, fast, itertools

def dist_W(P):   # Lambda_W = {p1=1, p2+p3<=1, p4=0}
    return np.maximum.reduce([1-P[:,0], P[:,3], np.maximum(P[:,1]+P[:,2]-1,0)/2])
def dist_N(P):   # Lambda_N = {(p1,0,1,0)}
    return np.maximum.reduce([P[:,1], 1-P[:,2], P[:,3]])
def dist_E(P):   # Lambda_E = {(1,0,1,p4)}
    return np.maximum.reduce([1-P[:,0], P[:,1], 1-P[:,2]])

WEDGE={'W':(-0.5,0.5,1,dist_W),'N':(0.2,1.5,1,dist_N),'E':(1.5,-0.2,-1,dist_E)}
h=0.05; g=np.linspace(0,1,int(round(1/h))+1)
P=np.array(list(itertools.product(g,repeat=4)))
eps=1e-9
for name,(u,v,sgn,dfun) in WEDGE.items():
    e=fast.eff_deficit(P,u,v,eps); r=fast.rival_margin(P,sgn,eps)
    d=np.minimum(e,r); D=dfun(P)
    bound=-D/(1+D)
    viol=d>bound+1e-7
    print(f'{name}: violations of  d <= -dist/(1+dist) : {viol.sum()} / {len(P)}')
    if viol.sum():
        idx=np.argsort(-(d-bound))[:5]
        for i in idx: print('   ',P[i],'d=',d[i],'bound=',bound[i],'dist=',D[i])
    # tightness: how close is the best d at each distance level to the bound
    lev=np.unique(np.round(D,10))[:6]
    for L in lev:
        m=np.isclose(D,L)
        print(f'    dist={L:.4f}: max d = {d[m].max():+.8f}   bound = {-L/(1+L):+.8f}')
