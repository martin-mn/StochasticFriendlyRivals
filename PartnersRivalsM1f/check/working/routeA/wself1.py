import sympy as sp
SWAP=(0,2,1,3)
p=sp.symbols('p1 p2 p3 p4')
def chain(sig,res):
    a=[sig[s] for s in range(4)]; b=[res[SWAP[s]] for s in range(4)]
    return sp.Matrix(4,4,lambda s,t:[a[s]*b[s],a[s]*(1-b[s]),(1-a[s])*b[s],(1-a[s])*(1-b[s])][t])
def stat_tree(M):
    """w_j proportional to cofactor C_jj of (I-M) (matrix-tree / Markov chain tree theorem)."""
    A=sp.eye(4)-M
    num=[]
    for j in range(4):
        idx=[i for i in range(4) if i!=j]
        num.append(sp.expand(A[idx,idx].det()))
    return num
M=chain(list(p),list(p))
num=stat_tree(M)
tot=sp.expand(sum(num))
print("numerators (unnormalised weights), common denominator = their sum:")
for i,x in enumerate(num):
    print(" ",['CC','CD','DC','DD'][i],"=",sp.factor(x))
print(" sum =",sp.factor(tot))
print(" CD-DC diff:",sp.simplify(num[1]-num[2]))
