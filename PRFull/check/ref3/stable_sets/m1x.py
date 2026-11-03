"""Exact eps->0 limits for memory-one pairs, in Fraction arithmetic (no floats, no tolerance).

Conventions (ms.tex): p = (p_CC, p_CD, p_DC, p_DD), focal action first; states CC=0, CD=1, DC=2, DD=3;
clamp p~ = eps + (1-2eps) p; the co-player sees CD and DC exchanged; (R,S,T,P) = (1,u,1+v,0).

limitW(sig, p): the eps->0 stationary distribution of the pair chain, from sig's (focal) view, by the
Markov-chain-tree theorem: d_X(eps) = det(I-M) with row/col X removed (polynomials in eps with Fraction
coefficients), W_X = lc(d_X)/lc(Dtot) if ord(d_X) = ord(Dtot), else 0.

gains(p): for each of the 16 binary deviators sigma, the limit gain pi(sigma,p) - pi(p,p) as an affine
form A*u + B*v + C (A = d_CD, B = d_DC, C = -d_CD - d_DD, d = W(sigma,p) - W(p,p)).
"""
from fractions import Fraction as F
from itertools import product

SW = (0, 2, 1, 3)
BIN = [tuple(F(b) for b in q) for q in product([0, 1], repeat=4)]
BINNAME = ["".join(str(int(x)) for x in q) for q in BIN]


def pmul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return r


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def psub(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)]


def det3(M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    return padd(psub(pmul(a, psub(pmul(e, i), pmul(f, h))),
                     pmul(b, psub(pmul(d, i), pmul(f, g)))),
                pmul(c, psub(pmul(d, h), pmul(e, g))))


def tree_weights(sig, p):
    L = []
    for s in range(4):
        x = [F(sig[s]), 1 - 2 * F(sig[s])]          # eps + (1-2eps) x  = x + eps(1-2x)
        nx = [1 - F(sig[s]), -(1 - 2 * F(sig[s]))]
        y = [F(p[SW[s]]), 1 - 2 * F(p[SW[s]])]
        ny = [1 - F(p[SW[s]]), -(1 - 2 * F(p[SW[s]]))]
        row = (pmul(x, y), pmul(x, ny), pmul(nx, y), pmul(nx, ny))
        L.append([psub([F(1)] if s == j else [F(0)], row[j]) for j in range(4)])
    out = []
    for k in range(4):
        idx = [i for i in range(4) if i != k]
        out.append(det3([[L[i][j] for j in idx] for i in idx]))
    return out


def order_lc(P):
    for i, c in enumerate(P):
        if c:
            return i, c
    return None, F(0)


def limitW(sig, p):
    ds = tree_weights(sig, p)
    Dt = padd(padd(ds[0], ds[1]), padd(ds[2], ds[3]))
    oD, cD = order_lc(Dt)
    W = []
    for X in range(4):
        oX, cX = order_lc(ds[X])
        W.append(F(0) if oX is None or oX > oD else cX / cD)
        assert oX is None or oX >= oD
    assert sum(W) == 1
    return tuple(W)


def gains(p):
    """[(name, A, B, C)] with gain = A u + B v + C, plus the self-play limit V."""
    V = limitW(p, p)
    out = []
    for name, q in zip(BINNAME, BIN):
        W = limitW(q, p)
        d = [W[k] - V[k] for k in range(4)]
        out.append((name, d[1], d[2], -d[1] - d[3], W))
    return V, out


def maxgain(G, u, v):
    return max(A * u + B * v + C for (_, A, B, C, _) in G)


def stable(G, u, v):
    return all(A * u + B * v + C <= 0 for (_, A, B, C, _) in G)


def gainers(G, u, v):
    return [(n, A * u + B * v + C) for (n, A, B, C, _) in G if A * u + B * v + C > 0]


def efficient(V, u, v):
    pi = V[0] + V[1] * u + V[2] * (1 + v)
    return pi == max(F(1), (1 + u + v) / 2)


def selfpay(V, u, v):
    return V[0] + V[1] * u + V[2] * (1 + v)
