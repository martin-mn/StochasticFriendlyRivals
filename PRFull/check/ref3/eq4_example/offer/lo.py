"""Independent exact eps->0 limits for memory-two pairs (Ref4 item 3).

Leading-order arithmetic: a quantity is (c, e) meaning c*eps^e + higher order, c > 0 a Fraction.
State reduction (GTH) forms only sums, products and quotients of positive quantities, so leading
terms never cancel and the limit is exact. Conventions (ms.tex SI 'Memory two', Methods):
state s = 4x+y, x = last outcome, y = the one before, CC=0 CD=1 DC=2 DD=3 focal first;
co-player cooperates with q[m(s)], m exchanges CD and DC in both remembered outcomes;
binary code = sum_s c_s 2^s, c_s = 1 for C (ALLD = 0, ALLC = 65535).
"""
from fractions import Fraction as F

SW = (0, 2, 1, 3)
M = tuple(4 * SW[s // 4] + SW[s % 4] for s in range(16))
ONE = F(1)


def code2q(c):
    return [(c >> s) & 1 for s in range(16)]


def q2code(q):
    return sum(int(v) << s for s, v in enumerate(q))


def genome(q):
    return ''.join(str(int(v)) for v in q)


def coop(v):
    """leading term of the clamped cooperation probability eps+(1-2eps)v"""
    v = F(v)
    return (v, 0) if v > 0 else (ONE, 1)


def dfct(v):
    v = F(v)
    return (1 - v, 0) if v < 1 else (ONE, 1)


def add(a, b):
    if a is None:
        return b
    if b is None:
        return a
    if a[1] < b[1]:
        return a
    if b[1] < a[1]:
        return b
    return (a[0] + b[0], a[1])


def mul(a, b):
    return (a[0] * b[0], a[1] + b[1])


def div(a, b):
    return (a[0] / b[0], a[1] - b[1])


def trans(p, q):
    T = [[None] * 16 for _ in range(16)]
    for s in range(16):
        a = s // 4
        pc, pd = coop(p[s]), dfct(p[s])
        qc, qd = coop(q[M[s]]), dfct(q[M[s]])
        for fD, pf in ((0, pc), (1, pd)):
            for qD, pq in ((0, qc), (1, qd)):
                o = 2 * fD + qD
                T[s][4 * o + a] = add(T[s][4 * o + a], mul(pf, pq))
    return T


def gth(T):
    A = [row[:] for row in T]
    n = 16
    for k in range(n - 1, 0, -1):
        S = None
        for j in range(k):
            S = add(S, A[k][j])
        assert S is not None
        colk = []
        for i in range(k):
            if A[i][k] is not None:
                A[i][k] = div(A[i][k], S)
                colk.append(i)
        rowk = [j for j in range(k) if A[k][j] is not None]
        for i in colk:
            aik = A[i][k]
            Ai = A[i]
            for j in rowk:
                Ai[j] = add(Ai[j], mul(aik, A[k][j]))
    pi = [None] * n
    pi[0] = (ONE, 0)
    for k in range(1, n):
        acc = None
        for i in range(k):
            if A[i][k] is not None:
                acc = add(acc, mul(pi[i], A[i][k]))
        pi[k] = acc
    emin = min(e for c, e in pi)
    tot = sum(c for c, e in pi if e == emin)
    return [c / tot if e == emin else F(0) for c, e in pi]


def limit(p, q):
    """W = (w_CC, w_CD, w_DC, w_DD) in the limit eps->0, from p's view (last outcome)."""
    pi = gth(trans(p, q))
    return tuple(sum(pi[4 * o + b] for b in range(4)) for o in range(4))


def margin(p, q):
    w = limit(p, q)
    return w[2] - w[1]


def pstar(x):
    """the point Eq. (4) of the main text with p_<DC,CD> = p_9 = x"""
    h, t = F(1, 2), F(1, 4)
    return [1, h, 1, h, 0, 0, 1, h, 1, F(x), 1, t, t, 0, t, 0]
