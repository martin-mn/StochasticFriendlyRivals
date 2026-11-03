"""V5: the alternator family p = (a, 0, 1, d), 0<=a<1, 0<d<=1, (a,d)!=(0,1):
     exact limit W(sigma,p) for the 16 deviators, the 16 facets, and the Nash polygon."""
import sympy as sp
from itertools import product
SWAP=(0,2,1,3)
e=sp.symbols('epsilon',positive=True)
a,d,u,v=sp.symbols('a d u v')

def clamp(x): return e+(1-2*e)*x
def stat(sig,res):
    A=[clamp(sig[s]) for s in range(4)]; B=[clamp(res[SWAP[s]]) for s in range(4)]
    M=sp.Matrix(4,4,lambda s,t:[A[s]*B[s],A[s]*(1-B[s]),(1-A[s])*B[s],(1-A[s])*(1-B[s])][t])
    X=(sp.eye(4)-M)
    num=[sp.expand(X[[i for i in range(4) if i!=j],[i for i in range(4) if i!=j]].det()) for j in range(4)]
    tot=sp.expand(sum(num))
    return [sp.cancel(n/tot) for n in num]

res=[a,sp.Integer(0),sp.Integer(1),d]
Ws=[sp.simplify(sp.limit(x,e,0,'+')) for x in stat(res,res)]
print("self-play limit W(p,p) =",[sp.simplify(x) for x in Ws])

rows={}
for bits in product([0,1],repeat=4):
    g=''.join(map(str,bits))
    W=[sp.simplify(sp.limit(x,e,0,'+')) for x in stat(list(bits),res)]
    rows[g]=W
    B=sp.simplify(W[1]-Ws[1]); C=sp.simplify(W[2]-Ws[2]); D=sp.simplify(W[3]-Ws[3])
    print(f"{g}  W={[sp.nsimplify(sp.factor(x)) for x in W]}")
    print(f"      facet: ({sp.factor(B)})(u-1) + ({sp.factor(C)}) v <= ({sp.factor(D)})")
import pickle
pickle.dump({k:[sp.srepr(x) for x in w] for k,w in rows.items()},
            open('/private/tmp/claude-501/-Users-martin-Documents-CL1-PartnersRivalsM1f/c1127d6a-179d-4b6f-be64-0786413868ce/scratchpad/A/alt.pkl','wb'))
