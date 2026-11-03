"""Global check of the structural theorem across all 81 strata.

THEOREM (binding states).  Let R = {s : W_s(p,p) > 0} be the limit-recurrent support of
self-play and I = {i : 0 < p_i < 1}.  A state s is BINDING if s in R and p_s is interior.
The advantage of C over D at state s in the deviator's average-reward MDP is
   Delta_s = beta + q_s (alpha - beta),   q_CC=p1, q_CD=p3, q_DC=p2, q_DD=p4,
so two binding states with different q force alpha = beta = 0, i.e. p is a limit equaliser.
Consequently a stratum can carry a Nash set of its own full dimension only if it has NO
binding state.  Prediction: full-dimensional  ==>  I cap R = empty."""
from itertools import product
import pickle
res = pickle.load(open('census.pkl', 'rb'))
STRATA = list(product([0,'i',1], repeat=4))
from fractions import Fraction as F
GU = [F(-3),F(-1),F(-1,4),F(0),F(1,4),F(1),F(3,2),F(3)]
GAMES = [(u,v) for u in GU for v in GU]
name = lambda s: ''.join('*' if c=='i' else str(c) for c in s)

viol = 0; nofull_with_free = []
tab = []
for si, k, pts in res:
    spec = STRATA[si]
    I = {i for i,c in enumerate(spec) if c=='i'}
    Rs = {frozenset(i for i,x in enumerate(W) if x != 0) for p,rows,W in pts}
    nobind = all(not (I & R) for R in Rs)
    Z = set(range(16))
    for p,rows,W in pts: Z &= {j for j,r in enumerate(rows) if r==(0,0,0)}
    everfull = any(any(all((r[0]*(u-1)+r[1]*v < r[2]) if k else (r[0]*(u-1)+r[1]*v <= r[2])
                       for j,r in enumerate(rows) if j not in Z) for p,rows,W in pts)
                   for u,v in GAMES)
    if everfull and not nobind:
        viol += 1
        print('VIOLATION', name(spec), 'I=',sorted(I), 'R=',[sorted(R) for R in Rs])
    tab.append((name(spec), k, nobind, everfull, sorted(Rs, key=len)[0]))
print(f'81 strata: {viol} violations of  (full-dimensional ==> no binding state)')
nb = [t for t in tab if t[2]]
print(f'{len(nb)} strata have no binding state; {sum(1 for t in tab if t[3])} are ever full-dimensional')
print('\nstrata with NO binding state, by dimension (these are the only possible carriers):')
for kk in (4,3,2,1,0):
    row = [t[0] for t in nb if t[1]==kk]
    print(f'  k={kk}: {len(row):2d}  {row}')
print('\nno-binding strata that are NEVER full-dimensional at any of the 64 games:')
print(' ', [t[0] for t in nb if not t[3]])
print('\nstates recurrent in self-play, for the maximal carriers:')
for nm in ['1***','***0','*01*','*10*','1**0','0**1']:
    t = [x for x in tab if x[0]==nm][0]
    print(f'  {nm}  R = {sorted(t[4])}   (0=CC 1=CD 2=DC 3=DD)')
