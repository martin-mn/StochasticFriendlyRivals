import lib, sympy as sp, itertools
from fractions import Fraction as F

def eff_cc_rule(p):
    p1,p2,p3,p4=p
    return (p1==1 and (p2,p3) not in ((0,1),(1,0))
            and (p4>0 or (1-p2)*(1-p3)==0))

def alt_rule(p):
    p1,p2,p3,p4=p
    return ((p2,p3) in ((0,1),(1,0)) and p1<1 and p4>0 and not (p1==0 and p4==1))

n=7
g=[F(i,n-1) for i in range(n)]
badX=badY=0; nX=nY=0
for p in itertools.product(g,repeat=4):
    X,Y,Z=lib.selfplay_limit(list(p))
    cx = (X==1); cy=(Y==1)
    nX+=cx; nY+=cy
    if cx!=eff_cc_rule(p): badX+=1; print('EFF_CC disagree',p,X)
    if cy!=alt_rule(p):    badY+=1; print('ALT disagree',p,Y)
print(f'grid {n}^4={n**4}: W_CC=1 -> {nX} points, rule disagreements {badX}')
print(f'grid {n}^4={n**4}: alternators -> {nY} points, rule disagreements {badY}')
# p2<->p3 symmetry of self-play
import random
random.seed(1); bad=0
for _ in range(200):
    p=[F(random.randrange(0,13),12) for _ in range(4)]
    q=[p[0],p[2],p[1],p[3]]
    if lib.selfplay_limit(p)!=lib.selfplay_limit(q): bad+=1
print('p2<->p3 self-play symmetry violations:',bad,'/200')
