"""Exact rational eps->0 leading-order engine (no floating point anywhere).

A copy of PartnersRivalsM2f/check/exact.py (same build() and gth(), same conventions), extended to
return what the certificates need:

    lead_exact(p, q) -> dict
        w16  : the limit distribution over the 16 pair states (exact rationals)
        W    : (w_CC, w_CD, w_DC, w_DD) from the focal player's view
        E    : the minimum eps-exponent over ALL spanning trees of the pair chain (integer)
        Ek   : Ek[k] = minimum exponent of a spanning tree rooted at k (integer)
        D    : D* = the eps^E coefficient of D(eps) = sum_k A_k(eps), A_k the tree sum rooted at k
        supp : bitmask of the states with Ek == E (the support of the limit)
    so that  w16[k] = A*_k / D*  with A*_k = D* w16[k] the eps^E coefficient of A_k.

How E and D* are obtained: GTH eliminates states 15..1; the pivot s_k is the exit probability of
state k in the chain censored on {0..k}, and prod_k s_k = det((I-P) without row/col 0) = A_0(eps)
(matrix-tree theorem).  All intermediates are sums, products and quotients of positive quantities,
so leading terms are exact: E_0 = sum_k exp(s_k), A_0* = prod_k coef(s_k), E_k = E_0 + exp(pi_k/pi_0)
and D* = A_0* * sum_{k in supp} coef(pi_k/pi_0).

Arithmetic: gmpy2.mpq when available (about 4x faster than fractions.Fraction), else Fraction;
results are returned as fractions.Fraction either way (set as_fraction=False to keep mpq).
Inputs may be int, Fraction, mpq, float (converted exactly) or decimal strings like '1/3'.
"""
from fractions import Fraction
try:
    import gmpy2
    _mpq = gmpy2.mpq
except Exception:          # pragma: no cover
    gmpy2 = None; _mpq = None

NS = 16
EINF = 10**9
def swapo(o): return ((o & 1) << 1) | (o >> 1)
def mirror(j): return 4*swapo(j >> 2) + swapo(j & 3)
MIRROR = [mirror(j) for j in range(NS)]

def _num(kind):
    if kind == 'fraction' or (kind == 'auto' and _mpq is None): return Fraction
    if _mpq is None: raise RuntimeError("gmpy2 not available")
    return _mpq

_MPQT = type(_mpq(0)) if _mpq is not None else ()
def _conv(x, Q):
    if isinstance(x, str): x = Fraction(x)
    elif not isinstance(x, (int, float, Fraction, _MPQT)):                   # numpy scalars etc.
        x = Fraction(x); x = Fraction(int(x.numerator), int(x.denominator))
    if Q is Fraction: return Fraction(x)
    if isinstance(x, Fraction): return _mpq(x.numerator, x.denominator)
    return _mpq(x)

def _tofrac(x):
    if isinstance(x, Fraction): return x
    return Fraction(int(x.numerator), int(x.denominator))

def _lead(x, one):
    """leading order (coef, exp) of (Pr(C), Pr(D)) for x in [0,1] under x -> eps+(1-2eps)x"""
    cC = (x, 0) if x > 0 else (one, 1)
    cD = (one - x, 0) if x < 1 else (one, 1)
    return cC, cD

def _build(p, q, one, zero):
    c = [[zero]*NS for _ in range(NS)]
    e = [[EINF]*NS for _ in range(NS)]
    for j in range(NS):
        (sC, sD) = _lead(p[j], one); (tC, tD) = _lead(q[MIRROR[j]], one)
        cs = (sC, sD); ct = (tC, tD)
        for dF in (0, 1):
            for dA in (0, 1):
                nxt = 4*(2*dF + dA) + (j >> 2)
                cc = cs[dF][0]*ct[dA][0]; ee = cs[dF][1] + ct[dA][1]
                if ee < e[j][nxt]: e[j][nxt] = ee; c[j][nxt] = cc
                elif ee == e[j][nxt]: c[j][nxt] += cc
    return c, e

def _gth(c, e, one, zero):
    E0 = 0; A0 = one
    for k in range(NS-1, 0, -1):
        es = EINF; cs = zero
        ck = c[k]; ek = e[k]
        for j in range(k):
            if ek[j] < es: es = ek[j]; cs = ck[j]
            elif ek[j] == es and es < EINF: cs += ck[j]
        E0 += es; A0 *= cs
        for i in range(k):
            if e[i][k] < EINF: c[i][k] /= cs; e[i][k] -= es
        for i in range(k):
            eik = e[i][k]
            if eik >= EINF: continue
            cik = c[i][k]; ci = c[i]; ei = e[i]
            for j in range(k):
                if ek[j] >= EINF: continue
                en = eik + ek[j]; cn = cik*ck[j]
                if en < ei[j]: ei[j] = en; ci[j] = cn
                elif en == ei[j]: ci[j] += cn
    pc = [zero]*NS; pe = [EINF]*NS; pc[0] = one; pe[0] = 0
    for k in range(1, NS):
        ek = EINF; ckk = zero
        for i in range(k):
            if e[i][k] >= EINF or pe[i] >= EINF: continue
            en = pe[i] + e[i][k]; cn = pc[i]*c[i][k]
            if en < ek: ek = en; ckk = cn
            elif en == ek: ckk += cn
        pc[k] = ckk; pe[k] = ek
    emin = min(pe)
    tot = zero
    for k in range(NS):
        if pe[k] == emin: tot += pc[k]
    w = [pc[k]/tot if pe[k] == emin else zero for k in range(NS)]
    Ek = [E0 + pe[k] if pe[k] < EINF else EINF for k in range(NS)]
    supp = sum(1 << k for k in range(NS) if pe[k] == emin)
    return w, E0 + emin, Ek, A0*tot, supp

def lead_exact(p, q, num='auto', as_fraction=True):
    """Exact leading-order limit of the pair (p, q); see the module docstring for the fields."""
    Q = _num(num); one = Q(1); zero = Q(0)
    P = [_conv(x, Q) for x in p]; R = [_conv(x, Q) for x in q]
    for x in P + R:
        if x < 0 or x > 1: raise ValueError("coordinate outside [0,1]")
    w, E, Ek, D, supp = _gth(*_build(P, R, one, zero), one, zero)
    W = [zero]*4
    for j in range(NS): W[j >> 2] += w[j]
    if as_fraction:
        w = [_tofrac(x) for x in w]; W = [_tofrac(x) for x in W]; D = _tofrac(D)
    return dict(w16=w, W=W, E=E, Ek=Ek, D=D, supp=supp)

def freqs(p, q, num='auto'):
    """exact (w_CC, w_CD, w_DC, w_DD), focal view"""
    return lead_exact(p, q, num)['W']

def margin(p, q, sign, num='auto'):
    """exact competitive margin sign*(w_DC - w_CD) of focal p against co-player q"""
    W = lead_exact(p, q, num)['W']
    return sign*(W[2] - W[1])

def selfplay(p, num='auto'):
    return lead_exact(p, p, num)['W']

def g2p(g): return [Fraction((g >> j) & 1) for j in range(NS)]

# ---------------------------------------------------------------- finite-eps reference (tests)
def finite_eps(p, q, eps):
    """Exact stationary data of the clamped pair chain at a fixed rational eps>0, computed without
    any leading-order arithmetic: returns (pi (16 Fractions), A (16 tree sums A_k(eps) as
    determinants of (I-P) with row/col k deleted)).  Independent check of E, D* and w (tests only)."""
    eps = Fraction(eps); P = [Fraction(x) for x in p]; R = [Fraction(x) for x in q]
    T = [[Fraction(0)]*NS for _ in range(NS)]
    for j in range(NS):
        a = eps + (1 - 2*eps)*P[j]; b = eps + (1 - 2*eps)*R[MIRROR[j]]
        sF = (a, 1 - a); sA = (b, 1 - b)
        for dF in (0, 1):
            for dA in (0, 1):
                T[j][4*(2*dF + dA) + (j >> 2)] += sF[dF]*sA[dA]
    L = [[(Fraction(1) if i == j else Fraction(0)) - T[i][j] for j in range(NS)] for i in range(NS)]
    def det(M):
        M = [r[:] for r in M]; n = len(M); d = Fraction(1)
        for c in range(n):
            piv = next((r for r in range(c, n) if M[r][c] != 0), None)
            if piv is None: return Fraction(0)
            if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
            d *= M[c][c]; inv = 1/M[c][c]
            for r in range(c+1, n):
                f = M[r][c]*inv
                if f:
                    Mr = M[r]; Mc = M[c]
                    for cc in range(c, n): Mr[cc] -= f*Mc[cc]
        return d
    A = []
    for k in range(NS):
        idx = [i for i in range(NS) if i != k]
        A.append(det([[L[i][j] for j in idx] for i in idx]))
    tot = sum(A)
    return [a/tot for a in A], A
