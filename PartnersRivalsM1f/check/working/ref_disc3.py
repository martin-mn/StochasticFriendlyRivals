import sys, numpy as np, itertools
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib as L
u,v=1.5,0.5      # T = S = 1.5 ; u+v=2 > 1 (above switch line); on the ray u-v=1
Emax=max(1.0,(u+(1+v))/2)
print("u,v",u,v,"T-S",1+v-u,"Emax",Emax)
print("undiscounted eps->0 self-play payoffs (convention A):")
for p in [(0,0,1,1),(0,0,1,0),(0.5,0,1,0.5),(1,0,1,0),(0,1,0,1),(0.3,0,1,0.7),(0,0,1,0.5)]:
    for eps in [1e-3,1e-5,1e-7]:
        w=L.stationary(np.array(p,float),np.array(p,float),eps)
        pi=L.payoff(w,u,v)
    print("  p=",p," pi(eps=1e-7)=",round(pi,8)," efficient?",abs(pi-Emax)<1e-5)
# does the project's own predicate agree?
names=[n for n in dir(L) if 'eff' in n.lower() or 'riv' in n.lower()]
print("m1lib predicates:",names)
