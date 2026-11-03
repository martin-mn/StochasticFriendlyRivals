"""Independent adversarial re-check of the discounted friendly-rival dossier.
Exact Fraction arithmetic.  Conventions: p=(p1..p4)=P(C|CC,CD,DC,DD), focal first,
plus explicit opening p0.  States CC=0 CD=1 DC=2 DD=3, SWAP=(0,2,1,3)."""
from fractions import Fraction as F
from itertools import product

SWAP = (0, 2, 1, 3)

def M_of(p, q):
    M = []
    for s in range(4):
        a, b = p[s], q[SWAP[s]]
        M.append([a*b, a*(1-b), (1-a)*b, (1-a)*(1-b)])
    return M

def occupancy(p0, p, q0, q, d):
    """w = (1-d) v0 (I - d M)^{-1}, exact."""
    v0 = [p0*q0, p0*(1-q0), (1-p0)*q0, (1-p0)*(1-q0)]
    M = M_of(p, q)
    # solve w (I - dM) = (1-d) v0  ->  (I - dM)^T w^T = (1-d) v0^T
    A = [[ (F(1) if i == j else F(0)) - d*M[j][i] for j in range(4)] for i in range(4)]
    b = [(1-d)*v0[i] for i in range(4)]
    # gaussian elimination
    n = 4
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]; b[c], b[piv] = b[piv], b[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]/A[c][c]
                for k in range(c, n): A[r][k] -= f*A[c][k]
                b[r] -= f*b[c]
    return [b[i]/A[i][i] for i in range(n)]

def payoffs(w, u, v):
    R, S, T, P = F(1), u, 1+v, F(0)
    return (w[0]*R + w[1]*S + w[2]*T + w[3]*P,
            w[0]*R + w[1]*T + w[2]*S + w[3]*P)

DETQ = [(tuple(F(x) for x in b), F(q0)) for b in product([0,1], repeat=4) for q0 in (0,1)]

def margin(p0, p, d, u, v):
    """min over the 32 deterministic co-players of pi(p,q)-pi(q,p)."""
    best = None
    for q, q0 in DETQ:
        w = occupancy(p0, p, q0, q, d)
        a, b = payoffs(w, u, v)
        m = a - b
        if best is None or m < best: best = m
    return best

def selfplay(p0, p, d, u, v):
    w = occupancy(p0, p, p0, p, d)
    return payoffs(w, u, v)[0]

def Emax(u, v):
    return max(F(1), (1+u+v)/2)
