"""Exact rational version of the eps->0 leading-order GTH engine of ../cannon/m2f.c.

Coefficients are Fractions and exponents are integers, so every limit frequency produced here is
an exact rational number: no floating point anywhere.  Used to certify the headline claims.
"""
from fractions import Fraction as F
NS = 16
EINF = 10**9
def swapo(o): return ((o & 1) << 1) | (o >> 1)
def mirror(j): return 4*swapo(j >> 2) + swapo(j & 3)

def lead(x):
    """leading order of (Pr(C), Pr(D)) for a probability x in [0,1], under clamp eps+(1-2eps)x."""
    cC = (x, 0) if x > 0 else (F(1), 1)
    cD = (1 - x, 0) if x < 1 else (F(1), 1)
    return cC, cD

def build(p, q):
    c = [[F(0)]*NS for _ in range(NS)]
    e = [[EINF]*NS for _ in range(NS)]
    for j in range(NS):
        jt = mirror(j)
        (sC, sD) = lead(p[j]); (tC, tD) = lead(q[jt])
        cs = [sC, sD]; ct = [tC, tD]
        for dF in (0, 1):
            for dA in (0, 1):
                nxt = 4*(2*dF + dA) + (j >> 2)
                cc = cs[dF][0]*ct[dA][0]; ee = cs[dF][1] + ct[dA][1]
                if ee < e[j][nxt]: e[j][nxt] = ee; c[j][nxt] = cc
                elif ee == e[j][nxt]: c[j][nxt] += cc
    return c, e

def gth(c, e):
    for k in range(NS-1, 0, -1):
        es = EINF; cs = F(0)
        for j in range(k):
            if e[k][j] < es: es = e[k][j]; cs = c[k][j]
            elif e[k][j] == es and es < EINF: cs += c[k][j]
        for i in range(k):
            if e[i][k] < EINF: c[i][k] /= cs; e[i][k] -= es
        for i in range(k):
            if e[i][k] >= EINF: continue
            eik = e[i][k]; cik = c[i][k]
            for j in range(k):
                if e[k][j] >= EINF: continue
                en = eik + e[k][j]; cn = cik*c[k][j]
                if en < e[i][j]: e[i][j] = en; c[i][j] = cn
                elif en == e[i][j]: c[i][j] += cn
    pc = [F(0)]*NS; pe = [EINF]*NS; pc[0] = F(1); pe[0] = 0
    for k in range(1, NS):
        ek = EINF; ck = F(0)
        for i in range(k):
            if e[i][k] >= EINF or pe[i] >= EINF: continue
            en = pe[i] + e[i][k]; cn = pc[i]*c[i][k]
            if en < ek: ek = en; ck = cn
            elif en == ek: ck += cn
        pc[k] = ck; pe[k] = ek
    emin = min(pe); tot = sum(pc[k] for k in range(NS) if pe[k] == emin)
    return [pc[k]/tot if pe[k] == emin else F(0) for k in range(NS)]

def freqs(p, q):
    w = gth(*build([F(x) for x in p], [F(x) for x in q]))
    W = [F(0)]*4
    for j in range(NS): W[j >> 2] += w[j]
    return W

def g2p(g): return [F((g >> j) & 1) for j in range(NS)]

def rival_margin(p, sign, stop_at_negative=True):
    best = None; arg = None
    for t in range(65536):
        W = freqs(p, g2p(t))
        d = sign*(W[2] - W[1])
        if best is None or d < best:
            best = d; arg = t
            if stop_at_negative and best < 0: break
    return best, arg
