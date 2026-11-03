import numpy as np, fast, itertools, sys

def grid(h):
    n = int(round(1/h))+1
    g = np.linspace(0,1,n)
    return np.array(list(itertools.product(g,repeat=4)))

def topd(P,u,v,eps,k=8):
    sign = 1 if (u-v)<1 else -1
    e = fast.eff_deficit(P,u,v,eps)
    r = fast.rival_margin(P,sign,eps)
    d = np.minimum(e,r)
    idx = np.argsort(-d)[:k]
    return [(P[i], d[i], e[i], r[i]) for i in idx]

WEDGE = {'W':(-0.5,0.5), 'S':(0.5,-0.5), 'E':(1.5,-0.2), 'N':(0.2,1.5)}
# check wedge membership: W u+v<1,u-v<1 ; S u+v<1,u-v>1 ; E u+v>1,u-v>1 ; N u+v>1,u-v<1
for k,(u,v) in WEDGE.items():
    print(k,'u,v=',u,v,' u+v=',u+v,' u-v=',u-v)
print()
for h in (0.05, 0.02):
    P = grid(h)
    print('=== h =',h,' grid',len(P))
    for wname in ('W','N','E','S'):
        u,v = WEDGE[wname]
        for p,d,e,r in topd(P,u,v,1e-7,5):
            print(f'  {wname}  p={np.round(p,4)}  d={d:.8f}  edef={e:.8f}  rm={r:.8f}   d/h={d/h:.4f}')
        print()
