"""Renewal accounting at eps = 0 around the focal decision state s* = <DC,CD> (s = 9), exact Fractions.
Also the error-free recovery of self-play after one error (main-text worked example)."""
from fractions import Fraction as F
from lo import M, code2q, pstar, margin, limit

NAMES = ['CC', 'CD', 'DC', 'DD']
SSTAR = 9
IDS = [36228, 36236, 40324, 40332, 52612, 52620, 56708, 56716]


def st(s):
    return '<%s,%s>' % (NAMES[s // 4], NAMES[s % 4])


def step0(p, q, s):
    """eps = 0 transitions from pair state s (focal view): list of (prob, next, outcome)."""
    a = s // 4
    out = []
    for fD in (0, 1):
        pf = F(p[s]) if fD == 0 else 1 - F(p[s])
        for qD in (0, 1):
            pq = F(q[M[s]]) if qD == 0 else 1 - F(q[M[s]])
            if pf * pq:
                o = 2 * fD + qD
                out.append((pf * pq, 4 * o + a, o))
    return out


def closed_classes(p, q):
    succ = {s: {n for _, n, _ in step0(p, q, s)} for s in range(16)}
    reach = {}
    for s in range(16):
        seen, st_ = {s}, [s]
        while st_:
            u = st_.pop()
            for v in succ[u]:
                if v not in seen:
                    seen.add(v); st_.append(v)
        reach[s] = seen
    cls = {frozenset(reach[s]) for s in range(16) if all(s in reach[v] for v in reach[s])}
    return [sorted(c) for c in cls]


def solve(A, b):
    n = len(A)
    A = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][n] / A[i][i] for i in range(n)]


def first_passage(p, q, target, reward):
    """expected (reward sum, rounds) from each state until the chain enters `target` (eps = 0).
    Solved on the states that reach target with probability 1 (checked)."""
    S = list(range(16))
    # value equations V(s) = sum pr*(r(o) + [n != target] V(n))
    A = [[F(0)] * 16 for _ in S]; bR = [F(0)] * 16; bL = [F(0)] * 16; bH = [F(0)] * 16
    for s in S:
        A[s][s] += 1
        for pr, n, o in step0(p, q, s):
            bR[s] += pr * reward(o); bL[s] += pr; bH[s] += pr * (1 if n == target else 0)
            if n != target:
                A[s][n] -= pr
    return A, bR, bL, bH


def excursion(p, q):
    """from s* split by the focal's action: (prob, co-player net exploitations, rounds) until s* recurs"""
    cls = closed_classes(p, q)
    # restrict to the closed class containing s*
    C = next(c for c in cls if SSTAR in c)
    idx = {s: i for i, s in enumerate(C)}
    n = len(C)
    def sys(reward):
        A = [[F(0)] * n for _ in range(n)]; b = [F(0)] * n
        for s in C:
            i = idx[s]; A[i][i] += 1
            for pr, nx, o in step0(p, q, s):
                b[i] += pr * reward(o)
                if nx != SSTAR:
                    A[i][idx[nx]] -= pr
        return solve(A, b)
    co_net = lambda o: {1: 1, 2: -1}.get(o, 0)       # co-player exploits focal at CD (focal view)
    G = sys(co_net); L = sys(lambda o: 1)
    res = []
    for pr, nx, o in step0(p, q, SSTAR):
        g = co_net(o) + (0 if nx == SSTAR else G[idx[nx]])
        l = 1 + (0 if nx == SSTAR else L[idx[nx]])
        res.append((pr, 'focal ' + ('C' if o < 2 else 'D') + ', outcome ' + NAMES[o], g, l))
    return cls, res


def paths(p, q, start_state, first_outcome, maxlen=12):
    """enumerate the eps = 0 paths from s* after the given first outcome, until s* recurs"""
    out = []
    def rec(s, seq, pr):
        if len(seq) > maxlen:
            out.append((pr, seq + ['...'])); return
        for pp, n, o in step0(p, q, s):
            if n == SSTAR:
                out.append((pr * pp, seq + [NAMES[o]]))
            else:
                rec(n, seq + [NAMES[o]], pr * pp)
    rec(4 * first_outcome + SSTAR // 4, [NAMES[first_outcome]], F(1))
    return out


if __name__ == '__main__':
    for x in [F(1, 4), F(5, 13), F(1, 7), F(9, 10)]:
        p = pstar(x)
        print('=== x = p_<DC,CD> =', x)
        for c in IDS:
            q = code2q(c)
            cls, res = excursion(p, q)
            assert len(cls) == 1 and SSTAR in cls[0]
            line = []
            for pr, lab, g, l in res:
                line.append('%s: prob %s, co-player net %s over %s rounds' % (lab, pr, g, l))
            Eg = sum(pr * g for pr, _, g, _ in res); El = sum(pr * l for pr, _, _, l in res)
            ren = -Eg / El
            ex = margin(p, q)
            assert ren == ex == (5 - 13 * x) / (2 * (25 - 13 * x))
            assert Eg == (13 * x - 5) / 8
            if c == IDS[0] or x == F(1, 4):
                print(c, 'closed class (focal view):', [st(s) for s in cls[0]])
                print('   ', ' | '.join(line), '| mean gain per decision', Eg, '= (13x-5)/8;',
                      'renewal margin', ren, '= exact limit', ex)
    # paths of the no-offer branch for one co-player (x-independent)
    p = pstar(F(1, 4)); q = code2q(36228)
    print('offer branch paths:', paths(p, q, SSTAR, 1))
    nb = paths(p, q, SSTAR, 3)
    print('no-offer branch: %d paths' % len(nb))
    for pr, seq in sorted(nb, key=lambda t: -t[0])[:12]:
        print('   ', pr, ','.join(seq))
    # which states of the closed class have a free coordinate of p?
    C = closed_classes(p, q)[0]
    print('free focal coordinates visited:', [(st(s), p[s]) for s in C if 0 < p[s] < 1])
