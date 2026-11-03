import numpy as np, fast, itertools
eps=1e-9
g=np.linspace(0,1,41)
def best_near(mask_fn,u,v,sgn,tag,n=400000,seed=11):
    rng=np.random.default_rng(seed)
    P=mask_fn(rng)
    e=fast.eff_deficit(P,u,v,eps); r=fast.rival_margin(P,sgn,eps); d=np.minimum(e,r)
    i=int(np.argmax(d)); print(f'  {tag}: max d = {d[i]:+.8f} at {np.round(P[i],8)}')
uN,vN=0.2,1.5; uE,vE=1.5,-0.2; uW,vW=-0.5,0.5
print('N wedge (T>S, above line):')
best_near(lambda r:(lambda P:(P.__setitem__((slice(None),1),1-r.random(400000)**6),
                              P.__setitem__((slice(None),2),r.random(400000)**6),P)[2])
          (r.random((400000,4))), uN,vN,1,'near the (p2,p3)=(1,0) corner face')
best_near(lambda r:(lambda P:(P.__setitem__((slice(None),1),r.random(400000)**6),
                              P.__setitem__((slice(None),2),1-r.random(400000)**6),P)[2])
          (r.random((400000,4))), uN,vN,1,'near the (p2,p3)=(0,1) corner face')
print('E wedge (T<S, above line):')
best_near(lambda r:(lambda P:(P.__setitem__((slice(None),1),1-r.random(400000)**6),
                              P.__setitem__((slice(None),2),r.random(400000)**6),P)[2])
          (r.random((400000,4))), uE,vE,-1,'near the (1,0) corner face')
best_near(lambda r:(lambda P:(P.__setitem__((slice(None),1),r.random(400000)**6),
                              P.__setitem__((slice(None),2),1-r.random(400000)**6),P)[2])
          (r.random((400000,4))), uE,vE,-1,'near the (0,1) corner face')
print('W wedge: near the deleted edge E={(p1,1,0,0)}')
best_near(lambda r:(lambda P:(P.__setitem__((slice(None),1),1-r.random(400000)**6),
                              P.__setitem__((slice(None),2),r.random(400000)**6),
                              P.__setitem__((slice(None),3),r.random(400000)**6),P)[3])
          (r.random((400000,4))), uW,vW,1,'p2->1,p3->0,p4->0')
