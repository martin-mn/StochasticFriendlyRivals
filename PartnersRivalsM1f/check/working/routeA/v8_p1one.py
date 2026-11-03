"""V8: the slice p1 = 1 (where every partner below the switch line lives).
Exact limit W(sigma,p) for p = (1,q2,q3,q4), q interior, and the 16 facets."""
import sympy as sp
from itertools import product
SWAP=(0,2,1,3)
e=sp.symbols('epsilon',positive=True)
q2,q3,q4=sp.symbols('q2 q3 q4',positive=True)
res=[sp.Integer(1),q2,q3,q4]
def clamp(x): return e+(1-2*e)*x
def statlim(sig,res):
    A=[clamp(sig[s]) for s in range(4)]; B=[clamp(res[SWAP[s]]) for s in range(4)]
    M=sp.Matrix(4,4,lambda s,t:[A[s]*B[s],A[s]*(1-B[s]),(1-A[s])*B[s],(1-A[s])*(1-B[s])][t])
    X=sp.eye(4)-M
    num=[sp.expand(X[[i for i in range(4) if i!=j],[i for i in range(4) if i!=j]].det()) for j in range(4)]
    tot=sp.expand(sum(num))
    return [sp.simplify(sp.limit(sp.cancel(n/tot),e,0,'+')) for n in num]
Ws=statlim(res,res)
print("self-play limit:",Ws)
print()
for bits in product([0,1],repeat=4):
    g=''.join(map(str,bits))
    W=statlim(list(bits),res)
    B=sp.factor(sp.simplify(W[1]-Ws[1])); C=sp.factor(sp.simplify(W[2]-Ws[2])); D=sp.factor(sp.simplify(W[3]-Ws[3]))
    print(f"{g}: W = {[sp.factor(x) for x in W]}")
    print(f"     ({B})(u-1) + ({C})v <= ({D})")
