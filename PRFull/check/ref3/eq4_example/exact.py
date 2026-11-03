"""Exact eps->0 limits by GTH over the field Q(x, eps) (no floating point)."""
import sympy as sp
from sympy.polys.domains import QQ
SW = [0, 2, 1, 3]
x, e = sp.symbols('x eps')
K = QQ.frac_field(x, e)
EX = [1, sp.Rational(1,2), 1, sp.Rational(1,2), 0, 0, 1, sp.Rational(1,2),
      1, sp.Rational(1,4), 1, sp.Rational(1,4), sp.Rational(1,4), 0, sp.Rational(1,4), 0]

def code2q(c):
    return [(c >> s) & 1 for s in range(16)]

def cst(s):
    a, b = divmod(s, 4); return 4*SW[a] + SW[b]

def trans(p, q):
    """p, q: lists of sympy numbers/expressions in x. Returns 16x16 list of K elements."""
    T = [[K.zero]*16 for _ in range(16)]
    E = K.from_sympy(e)
    for s in range(16):
        a = s // 4
        pc = E + (1 - 2*E) * K.from_sympy(sp.sympify(p[s]))
        qc = E + (1 - 2*E) * K.from_sympy(sp.sympify(q[cst(s)]))
        for fD in (0, 1):
            pf = pc if fD == 0 else 1 - pc
            for qD in (0, 1):
                pq = qc if qD == 0 else 1 - qc
                o = 2*fD + qD
                T[s][4*o + a] += pf * pq
    return T

def gth(T):
    A = [row[:] for row in T]
    n = len(A)
    for k in range(n-1, 0, -1):
        s = sum((A[k][j] for j in range(k)), K.zero)
        for i in range(k):
            A[i][k] = A[i][k] / s
        for i in range(k):
            if A[i][k] == K.zero: continue
            for j in range(k):
                if A[k][j] == K.zero: continue
                A[i][j] += A[i][k] * A[k][j]
    pi = [K.zero]*n; pi[0] = K.one
    for k in range(1, n):
        pi[k] = sum((pi[i]*A[i][k] for i in range(k)), K.zero)
    tot = sum(pi, K.zero)
    return [v/tot for v in pi]

def lim0(f):
    """limit eps->0 of an element of K, returned as sympy expression in x."""
    ex = K.to_sympy(f)
    num, den = sp.fraction(sp.together(ex))
    pn = sp.Poly(sp.expand(num), e); pd = sp.Poly(sp.expand(den), e)
    # lowest-order coefficients
    def low(P):
        terms = sorted(P.terms(), key=lambda t: t[0][0])
        for (k,), c in terms:
            if sp.simplify(c) != 0: return k, c
        return None, 0
    kn, cn = low(pn); kd, cd = low(pd)
    if kn is None: return sp.Integer(0)
    assert kn >= kd, (kn, kd)
    return sp.factor(cn/cd) if kn == kd else sp.Integer(0)

def outcome_freq(pi):
    return [sum((pi[4*o+b] for b in range(4)), K.zero) for o in range(4)]

def margin_limit(p, q):
    w = outcome_freq(gth(trans(p, q)))
    return lim0(w[2] - w[1]), [lim0(v) for v in w]
