"""V1: the 16 limit-W formulas and the affine (polygon) structure of the limit-Nash region."""
import sys, pickle, numpy as np, sympy as sp
sys.path.insert(0, '/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib
from itertools import product

p1,p2,p3,p4 = sp.symbols('p1 p2 p3 p4')
P=(p1,p2,p3,p4)

# --- the 16 limit distributions W(sigma,p), interior p, from wsig.py (eps->0 of the clamped chain)
Wsig = {
 '0000': (0, 0, -p4/(p2-p4-1), (p2-1)/(p2-p4-1)),
 '0001': (-p4*(p2-1)/(p1*p4-2*p2-p3*p4+p3+2), (p2-1)*(p4-1)/(p1*p4-2*p2-p3*p4+p3+2),
          (p1*p4-p3*p4+p3)/(p1*p4-2*p2-p3*p4+p3+2), -(p2-1)/(p1*p4-2*p2-p3*p4+p3+2)),
 '0010': (-p2*p4/(p1*p2-p2*p3+p3-2*p4-1), p4*(p2-1)/(p1*p2-p2*p3+p3-2*p4-1),
          -p4/(p1*p2-p2*p3+p3-2*p4-1), (p1*p2-p2*p3+p3-1)/(p1*p2-p2*p3+p3-2*p4-1)),
 '0011': (-(p2*p3-p3*p4+p4)/(2*(p1*p2-p1*p4-p2*p3+p3*p4-1)),
          (p1*p2-p1*p4+p4-1)/(2*(p1*p2-p1*p4-p2*p3+p3*p4-1)),
          -(p1*p4-p3*p4+p3)/(2*(p1*p2-p1*p4-p2*p3+p3*p4-1)),
          (p1*p2-p2*p3+p3-1)/(2*(p1*p2-p1*p4-p2*p3+p3*p4-1))),
 '0100': (0, 0, -p4/(p2-p4-1), (p2-1)/(p2-p4-1)),
 '0101': (-p3*(p2-1)/(p1*p3-2*p2*p3+p2*p4-p2+2*p3-p4+1), (p2-1)*(p4-1)/(p1*p3-2*p2*p3+p2*p4-p2+2*p3-p4+1),
          p1*p3/(p1*p3-2*p2*p3+p2*p4-p2+2*p3-p4+1), -p3*(p2-1)/(p1*p3-2*p2*p3+p2*p4-p2+2*p3-p4+1)),
 '0110': (-p3*p4/(p1*p3+p2*p4-2*p3*p4-p3-p4), p4*(p2-1)/(p1*p3+p2*p4-2*p3*p4-p3-p4),
          -p3*p4/(p1*p3+p2*p4-2*p3*p4-p3-p4), p3*(p1-1)/(p1*p3+p2*p4-2*p3*p4-p3-p4)),
 '0111': (-p3/(p1*p2-p1*p4-2*p3+p4-1), (p1*p2-p1*p4+p4-1)/(p1*p2-p1*p4-2*p3+p4-1),
          -p1*p3/(p1*p2-p1*p4-2*p3+p4-1), p3*(p1-1)/(p1*p2-p1*p4-2*p3+p4-1)),
 '1000': (0, 0, -p4/(p2-p4-1), (p2-1)/(p2-p4-1)),
 '1001': (-p4*(p2-1)/(2*p1*p2-p1*p3-2*p1-p2*p4-2*p2+p3+p4+2), (p1-1)*(p2-1)/(2*p1*p2-p1*p3-2*p1-p2*p4-2*p2+p3+p4+2),
          -p3*(p1-1)/(2*p1*p2-p1*p3-2*p1-p2*p4-2*p2+p3+p4+2), (p1-1)*(p2-1)/(2*p1*p2-p1*p3-2*p1-p2*p4-2*p2+p3+p4+2)),
 '1010': (p2*p4/(p1*p3-2*p1*p4-p1+p2*p4-p3+2*p4+1), -p4*(p1-1)/(p1*p3-2*p1*p4-p1+p2*p4-p3+2*p4+1),
          -p4*(p1-1)/(p1*p3-2*p1*p4-p1+p2*p4-p3+2*p4+1), (p1-1)*(p3-1)/(p1*p3-2*p1*p4-p1+p2*p4-p3+2*p4+1)),
 '1011': (-(p2*p3-p3*p4+p4)/(2*p1-p2*p3+p3*p4-p4-2), (p1-1)/(2*p1-p2*p3+p3*p4-p4-2),
          p3*(p1-1)/(2*p1-p2*p3+p3*p4-p4-2), -(p1-1)*(p3-1)/(2*p1-p2*p3+p3*p4-p4-2)),
 '1100': (-p3/(2*(p1-p3-1)), (p1-1)/(2*(p1-p3-1)), -p4/(2*(p2-p4-1)), (p2-1)/(2*(p2-p4-1))),
 '1101': (-p3/(p1-p3-1), (p1-1)/(p1-p3-1), 0, 0),
 '1110': (-p3/(p1-p3-1), (p1-1)/(p1-p3-1), 0, 0),
 '1111': (-p3/(p1-p3-1), (p1-1)/(p1-p3-1), 0, 0),
}

# --- self-play W(p,p), interior p (matrix-tree cofactors, common factor (p2-p3-1) divided out)
_num = (-p4*(2*p2*p3-p2*p4-p3*p4+p4),
         p4*(p1-1)*(p1-p4+1),
         p4*(p1-1)*(p1-p4+1),
        -(p1-1)*(p1*p2+p1*p3-p1-2*p2*p3+p2+p3-1))
_tot = sum(_num)
Wself = tuple(sp.cancel(x/_tot) for x in _num)

GEN = [''.join(str(b) for b in bits) for bits in product([0,1],repeat=4)]
f_sig = {g: sp.lambdify(P, Wsig[g], 'numpy') for g in GEN}
f_self = sp.lambdify(P, Wself, 'numpy')

def facets(pv):
    """16 triples (B,C,D): the limit-Nash region of p is  { B(u-1) + C v <= D }."""
    Ws = np.array(f_self(*pv), float)
    out = {}
    for g in GEN:
        W = np.array(f_sig[g](*pv), float)
        d = W - Ws
        out[g] = (d[1], d[2], d[3])
    return out, Ws

def nash_pred(pv, u, v):
    F,_ = facets(pv)
    return max(B*(u-1)+C*v-D for (B,C,D) in F.values())

if __name__ == '__main__':
    rng = np.random.default_rng(20260909)
    # (a) do the symbolic limit W(sigma,p) agree with the numerical chain at small eps?
    worst = 0.0
    for _ in range(300):
        pv = rng.uniform(0.05,0.95,4)
        for g,sig in zip(GEN, m1lib.DET):
            Wn = m1lib.stationary(sig, pv, 1e-9)
            Wa = np.array(f_sig[g](*pv), float)
            worst = max(worst, np.abs(Wn-Wa).max())
        Wn = m1lib.self_play(pv, 1e-9); Wa = np.array(f_self(*pv), float)
        worst = max(worst, np.abs(Wn-Wa).max())
    print("(a) max |W_symbolic - W_numeric(eps=1e-9)| over 300 interior p x 17 chains:", worst)

    # (b) does the affine form reproduce nash_gap?  gap = -max_sigma[B(u-1)+Cv-D]
    worst2 = 0.0
    for _ in range(200):
        pv = rng.uniform(0.05,0.95,4)
        for _ in range(10):
            u,v = rng.uniform(-6,6,2)
            g1 = m1lib.nash_gap(pv,u,v,1e-9)
            g2 = -nash_pred(pv,u,v)
            worst2 = max(worst2, abs(g1-g2))
    print("(b) max |nash_gap(eps=1e-9) - affine prediction| over 2000 (p,u,v):", worst2)

    # (c) convexity: the region is an intersection of 16 half planes -> check it is convex
    #     by sampling: if two points are limit-Nash so is their midpoint.
    bad = 0; tested = 0
    for _ in range(400):
        pv = rng.uniform(0.05,0.95,4)
        F,_ = facets(pv)
        A = np.array([[B,C] for (B,C,D) in F.values()]); b = np.array([B+D for (B,C,D) in F.values()])
        # region: A @ (u,v) <= b   (B u + C v <= D + B)
        pts=[]
        for _ in range(400):
            z = rng.uniform(-8,8,2)
            if np.all(A@z <= b + 1e-12): pts.append(z)
        if len(pts) >= 2:
            for _ in range(50):
                i,j = rng.integers(0,len(pts),2)
                m = 0.5*(pts[i]+pts[j]); tested += 1
                if not np.all(A@m <= b + 1e-9): bad += 1
    print("(c) convexity midpoint test:", tested, "midpoints,", bad, "violations")
