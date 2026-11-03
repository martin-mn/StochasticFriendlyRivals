import numpy as np
from fractions import Fraction as F
import itertools
sw=[0,2,1,3]
def p0(x9):
    return [1,F(1,2),1,F(1,2), 0,0,1,F(1,2), 1,x9,1,F(1,4), F(1,4),0,F(1,4),0]
# --- self-play after one error, eps=0: exact absorbing chain with exploitation-difference tracking
def selfplay(p):
    # state: (x,y) for A, plus d = (#A exploited B) - (#B exploited A); count rounds of CD and DC for A
    from collections import defaultdict
    start=((2,0),1)  # A defected by error: A's last outcome DC, before CC; d=1
    # value iteration on expected counts (exact fractions, absorbing chain, finite since probabilities)
    # build reachable states
    states=set(); frontier=[start]
    trans={}
    while frontier:
        s=frontier.pop()
        if s in states: continue
        states.add(s)
        (x,y),d=s
        if (x,y)==(0,0): continue
        pa=p[4*x+y]; pb=p[4*sw[x]+sw[y]]
        out=[]
        for ca,pra in ((1,pa),(0,1-pa)):
            for cb,prb in ((1,pb),(0,1-pb)):
                pr=pra*prb
                if pr==0: continue
                o=2*(1-ca)+(1-cb)
                nd=d+(1 if o==2 else -1 if o==1 else 0)
                ns=((o,x),nd); out.append((ns,pr,o))
                frontier.append(ns)
        trans[s]=out
    S=sorted(states,key=str); idx={s:i for i,s in enumerate(S)}
    import sympy as sp
    n=len(S)
    # expected number of rounds with o==1 (A exploited) from each state, and absorption distribution of d
    M=sp.zeros(n,n); bA=sp.zeros(n,1); bB=sp.zeros(n,1)
    for s in S:
        i=idx[s]; M[i,i]=1
        if s not in trans: continue
        for ns,pr,o in trans[s]:
            M[i,idx[ns]]-=sp.Rational(pr.numerator,pr.denominator) if isinstance(pr,F) else pr
            if o==1: bA[i]+=pr
            if o==2: bB[i]+=pr
    eA=M.LUsolve(bA); eB=M.LUsolve(bB)
    absorb={s for s in S if s not in trans}
    print('reachable states',n,'absorbing d values',sorted({s[1] for s in absorb}))
    i=idx[start]
    print('expected times A (the erring player) is exploited after the error:',eA[i],' B exploited (incl. the error):',eB[i])
selfplay(p0(F(1,4)))
# --- margins against all binary co-players, small eps
def margins(p,eps=1e-8):
    p=np.array([float(v) for v in p]); pt=eps+(1-2*eps)*p
    Q=np.array(list(itertools.product([0,1],repeat=16)))[:, ::-1].astype(float)  # all binary co-players
    Q=eps+(1-2*eps)*Q
    N=len(Q); P=np.zeros((N,16,16))
    for x in range(4):
        for y in range(4):
            s=4*x+y; pa=pt[s]; pb=Q[:,4*sw[x]+sw[y]]
            for ca,pra in ((1,pa),(0,1-pa)):
                for cb in (1,0):
                    prb=pb if cb else 1-pb
                    o=2*(1-ca)+(1-cb)
                    P[:,s,4*o+x]+=pra*prb
    A=np.transpose(P,(0,2,1))-np.eye(16); A[:,-1,:]=1
    b=np.zeros((N,16)); b[:,-1]=1
    pi=np.linalg.solve(A,b)
    w=pi.reshape(N,4,4).sum(axis=2)   # last-round outcome
    return (w[:,2]-w[:,1]), Q
for x9 in [F(1,4),F(1,3),F(5,13)-F(1,100),F(5,13)+F(1,100),F(2,5),F(1,2)]:
    m,Q=margins(p0(x9))
    k=np.argmin(m)
    print('p_<DC,CD>=%-8s min margin %.3e  argmin co-player %s'%(x9,m[k],''.join(str(int(round(v))) for v in Q[k])))
lo,hi=0.37,0.40
for _ in range(22):
    mid=(lo+hi)/2
    m,_=margins(p0(F(mid)),eps=1e-10)
    if m.min()>-1e-7: lo=mid
    else: hi=mid
print('threshold in [%.7f, %.7f]; 5/13 = %.7f'%(lo,hi,5/13))
