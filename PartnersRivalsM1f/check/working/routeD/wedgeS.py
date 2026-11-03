import numpy as np, fast, itertools
def frs_rule(P,tol=1e-12):
    p1,p2,p3,p4=P.T
    return ((p1>=1-tol)&(p2+p3>=1-tol)&(p2>tol)&(p3>tol)
            &((p4>tol)|(np.maximum(p2,p3)>=1-tol)))
h=0.05; g=np.linspace(0,1,int(round(1/h))+1)
P=np.array(list(itertools.product(g,repeat=4)))
eps=1e-9
# --- S wedge: FR_S closed form vs numerics
u,v=0.5,-0.5
e=fast.eff_deficit(P,u,v,eps); r=fast.rival_margin(P,-1,eps)
d=np.minimum(e,r)
pred=frs_rule(P); obs=d>-1e-6
print('S wedge, grid h=0.05:  predicted FR', pred.sum(), ' observed FR', obs.sum(),
      ' disagreements', (pred!=obs).sum())
if (pred!=obs).sum():
    for i in np.where(pred!=obs)[0][:8]: print('   ',P[i],d[i],pred[i],obs[i])
# --- switch line with T>S: u+v=1, u-v<1.   efficiency deficit = -W_DD there
for (u,v) in [(0.3,0.7),(-0.5,1.5),(0.0,1.0)]:
    assert abs(u+v-1)<1e-12 and u-v<1
    e=fast.eff_deficit(P,u,v,eps); r=fast.rival_margin(P,+1,eps)
    d=np.minimum(e,r)
    i=int(np.argmax(d))
    print(f'switch line u={u},v={v} (T>S): max d = {d.max():+.8f} at p={P[i]}  '
          f'(edef={e[i]:+.8f}, rm={r[i]:+.8f})')
# --- switch line with T<S
for (u,v) in [(1.5,-0.5)]:
    e=fast.eff_deficit(P,u,v,eps); r=fast.rival_margin(P,-1,eps)
    d=np.minimum(e,r)
    print(f'switch line u={u},v={v} (T<S): max d = {d.max():+.8f}, #FR = {(d>-1e-6).sum()}')
