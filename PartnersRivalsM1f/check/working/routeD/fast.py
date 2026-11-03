"""Vectorised near-limit evaluation.  Self-play uses the cancellation-free closed form;
the rival margin uses batched 4x4 solves at a moderately small eps."""
import numpy as np, itertools
SWAP = (0, 2, 1, 3)
VERT = np.array(list(itertools.product([0, 1], repeat=4)), float)

def selfplay_XYZ(P, eps):
    """P: (n,4).  Returns X=W_CC, Y=W_CD+W_DC, Z=W_DD, each (n,).  No cancellation."""
    a, b, c, d = [eps + (1-2*eps)*P[:, i] for i in range(4)]
    beta = (1-b)*(1-c); gam = b*c
    Xn = d*(d*beta + gam*(2-d))
    Yn = 2*d*(1-a)*(1+a-d)
    Zn = (1-a)*((1+a)*beta + gam*(1-a))
    S = Xn+Yn+Zn
    return Xn/S, Yn/S, Zn/S

def margins(P, eps):
    """(n,16) array of w_DC - w_CD for focal P against each of the 16 vertices."""
    n = len(P)
    pc = eps + (1-2*eps)*P                       # (n,4)
    out = np.empty((n, 16))
    for k, q in enumerate(VERT):
        qc = eps + (1-2*eps)*q                   # (4,)
        M = np.empty((n, 4, 4))
        for s in range(4):
            A = pc[:, s]; B = qc[SWAP[s]]
            M[:, s, 0] = A*B; M[:, s, 1] = A*(1-B)
            M[:, s, 2] = (1-A)*B; M[:, s, 3] = (1-A)*(1-B)
        L = np.transpose(M, (0, 2, 1)) - np.eye(4)
        L[:, 3, :] = 1.0
        rhs = np.zeros((n, 4)); rhs[:, 3] = 1.0
        w = np.linalg.solve(L, rhs)
        out[:, k] = w[:, 2] - w[:, 1]
    return out

def rival_margin(P, sign, eps):
    return (sign*margins(P, eps)).min(axis=1)

def eff_deficit(P, u, v, eps):
    X, Y, Z = selfplay_XYZ(P, eps)
    pi = X + 0.5*Y*(1+u+v)
    return pi - max(1.0, (1+u+v)/2)
