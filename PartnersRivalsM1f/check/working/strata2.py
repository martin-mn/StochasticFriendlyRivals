from fractions import Fraction as F
from itertools import product
import pickle, numpy as np
res = pickle.load(open('census.pkl','rb'))
STRATA = list(product([0,'i',1], repeat=4))
GU = [F(-3),F(-1),F(-1,4),F(0),F(1,4),F(1),F(3,2),F(3)]
GAMES=[(u,v) for u in GU for v in GU]
name=lambda s:''.join('*' if c=='i' else str(c) for c in s)
IDX={name(s):i for i,s in enumerate(STRATA)}
info={name(STRATA[si]):(k,pts) for si,k,pts in res}

def fullmap(nm):
    k,pts=info[nm]
    Z=set(range(16))
    for p,rows,W in pts: Z &= {j for j,r in enumerate(rows) if r==(0,0,0)}
    out={}
    for u,v in GAMES:
        out[(u,v)]=any(all((r[0]*(u-1)+r[1]*v < r[2]) if k else (r[0]*(u-1)+r[1]*v <= r[2])
                       for j,r in enumerate(rows) if j not in Z) for p,rows,W in pts)
    return out

print('self-play limit-recurrent support W(p,p) on the candidate carrier strata:')
for nm in ['1***','***0','*01*','*10*','1**0','0**1','****']:
    k,pts=info[nm]; Ws={tuple(W) for p,rows,W in pts}
    print(f'  {nm}  k={k}  W(p,p) in {[tuple(str(x) for x in w) for w in list(Ws)[:2]]}')

print('\nwhere each candidate carrier is FULL-dimensional (grid of 64 games):')
Q=lambda u,v: (u>=1 and v>=0)
for nm in ['1***','***0','*01*','*10*','1**0','0**1','****','0***','*0**','**0*','***1','*1**','**1*']:
    fm=fullmap(nm); gs=[g for g in GAMES if fm[g]]
    inQ=[g for g in gs if Q(*g)]
    pred=''
    if nm=='1***': pred='  [predicted: u<1 or v<0, plus (1,0)]'
    if nm=='***0': pred='  [predicted: u<0]'
    print(f'  {nm:5s} k={info[nm][0]}  full at {len(gs):2d}/64, of which {len(inQ):2d} in the quadrant u>=1,v>=0{pred}')
    if 0 < len(gs) <= 20: print(f'         games: {[(str(a),str(b)) for a,b in gs]}')
