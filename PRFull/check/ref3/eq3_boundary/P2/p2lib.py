"""Prover 2: eps->0 limit of the deviator/resident pair chain by REGULAR PERTURBATION of the
unperturbed (eps = 0) chain.  No cofactors, no Markov-chain-tree formula.

Conventions.  States 0=CC, 1=CD, 2=DC, 3=DD, deviator's action first.  The resident p sees the
state with CD <-> DC swapped.  Each intended action is flipped with prob eps, i.e. a strategy q
is played as eps + (1-2 eps) q.  M_eps = M0 + eps M1 + eps^2 M2 exactly (each entry is a product
of two affine functions of eps).

THE LIMIT LEMMA (proved in README / docstring of limit_W, used for every stratum and deviator):
For every eps in (0,1/2) all 16 entries of M_eps are positive, so pi_eps is unique.  Let pi* be
any accumulation point of pi_eps as eps->0.  (i) pi* M0 = pi*, so pi* = sum_i w_i mu_i, where
C_1..C_k are the closed communicating classes of M0 and mu_i is the unique stationary law of
M0 restricted to C_i (irreducible, periodic or not).  (ii) Let h_j be the M0-absorption
probability into C_j (h_j = 1 on C_j, 0 on other closed classes, the unique solution of
(I - M0_TT) h_T = M0_{T,C_j} 1 on the transient set T).  M0 h_j = h_j, so
0 = pi_eps (M_eps - I) h_j = eps pi_eps M1 h_j + eps^2 pi_eps M2 h_j  exactly; dividing by eps
and letting eps -> 0 gives  pi* M1 h_j = 0,  i.e.  w Q = 0  with  Q_ij = mu_i M1 h_j.
Q_ij >= 0 for i != j (M1 >= 0 off the class, h_j = 0 on C_i) and the rows of Q sum to 0
(sum_j h_j = 1, M1 1 = 0), so Q is a generator.  (iii) If Q has a unique stationary law w
(the directed graph of its positive off-diagonal entries has exactly one closed class), every
accumulation point equals sum_i w_i mu_i, and since the simplex is compact pi_eps converges to
it.  If Q has several closed classes the first-order equations do not fix w and a second-order
step is needed (limit_W raises; it never happens on the strata treated here, see the run).

What makes this valid at EVERY parameter value of an open stratum (not generically): the
support of M0 depends only on which of a, b, c are 0, 1 or in (0,1), so the closed classes and
the transient set are the same at every point of the stratum; the linear systems for mu_i and
h_j are then nonsingular at every point; every denominator that appears is certified positive
on the whole open stratum; and every off-diagonal Q_ij is either identically 0 or certified
positive, so the closed-class structure of Q is also constant on the stratum.
"""
import itertools
import sympy as sp

A, B, C = sp.symbols('a b c', positive=True)
EPS = sp.Symbol('eps', positive=True)
SW = (0, 2, 1, 3)
NAMES = ['CC', 'CD', 'DC', 'DD']
SIGMAS = [''.join(map(str, t)) for t in itertools.product([0, 1], repeat=4)]


def stratum_values(code):
    """code e.g. '0*1' for (a,b,c) -> (a,b,c) sympy values and list of free symbols."""
    vals, free = [], []
    for ch, sym in zip(code, (A, B, C)):
        if ch == '0':
            vals.append(sp.Integer(0))
        elif ch == '1':
            vals.append(sp.Integer(1))
        else:
            vals.append(sym); free.append(sym)
    return tuple(vals), free


def chain(sig, p, eps=EPS):
    """Exact transition matrix M_eps (sympy) of deviator sig (4 numbers) vs resident p."""
    M = sp.zeros(4, 4)
    for s in range(4):
        x = eps + (1 - 2 * eps) * sig[s]
        y = eps + (1 - 2 * eps) * p[SW[s]]
        M[s, 0] = x * y; M[s, 1] = x * (1 - y); M[s, 2] = (1 - x) * y; M[s, 3] = (1 - x) * (1 - y)
    return M.applyfunc(sp.expand)


def expansion(sig, p):
    M = chain(sig, p)
    M0 = M.subs(EPS, 0).applyfunc(sp.expand)
    M1 = M.diff(EPS).subs(EPS, 0).applyfunc(sp.expand)
    M2 = (M.diff(EPS, 2) / 2).applyfunc(sp.expand)
    assert (M - M0 - EPS * M1 - EPS ** 2 * M2).applyfunc(sp.expand) == sp.zeros(4, 4)
    return M0, M1, M2


# ------------------------------------------------------------------ sign certificates
def t_transform(f, free):
    """For a polynomial f in the free params x_i in (0,1), return the polynomial
    F(t) = f(t/(1+t)) * prod (1+t_i)^{deg_i f}  (a Poly in t_i).  Since prod(1+t_i)^d > 0 and
    x -> t/(1+t) maps (0,inf) onto (0,1), sign f = sign F on the open box."""
    ts = sp.symbols('t0:%d' % len(free), positive=True)
    P = sp.Poly(sp.expand(f), *free)
    degs = [P.degree(x) for x in free]
    g = sp.expand(f)
    g = g.subs({x: t / (1 + t) for x, t in zip(free, ts)}, simultaneous=True)
    for t, d in zip(ts, degs):
        g = g * (1 + t) ** d
    g = sp.cancel(sp.together(g))
    num, den = sp.fraction(g)
    assert den.is_number and den > 0, den
    return sp.Poly(sp.expand(num), *ts), ts


def poly_sign_cert(expr, free):
    """Certificate that a POLYNOMIAL expr in the free params is >0 ('+'), <0 ('-') or
    identically 0 ('0') on the open box (0,1)^n; (None, why) if no certificate found.
    Method: x = t/(1+t), clear prod(1+t)^deg; all coefficients of one sign (and not all zero)
    => sign-definite and never zero for t > 0."""
    expr = sp.expand(expr)
    if expr == 0:
        return '0', 'identically zero'
    if not (set(expr.free_symbols) - set()) or not free:
        assert expr.is_number
        return ('+' if expr > 0 else '-'), 'constant %s' % expr
    F, ts = t_transform(expr, free)
    coeffs = [sp.Rational(c) for c in F.coeffs()]
    if all(c >= 0 for c in coeffs):
        return '+', 't-subst: %d coeffs >= 0' % len(coeffs)
    if all(c <= 0 for c in coeffs):
        return '-', 't-subst: %d coeffs <= 0' % len(coeffs)
    return None, 'mixed coefficient signs'


def factor_sign_cert(expr, free):
    """Sign certificate for a RATIONAL function: factor numerator and denominator into
    irreducibles; each factor must carry a poly_sign_cert (manifest factors such as a, 1-a,
    b, 1-b, c, 1-c are covered: e.g. 1-a -> 1/(1+t) cleared -> coefficient 1).
    Returns (sign, list of (factor, multiplicity, sign, how)) or (None, ...)."""
    expr = sp.factor(sp.cancel(sp.together(expr)))
    if expr == 0:
        return '0', [('0', 1, '0', 'identically zero')]
    num, den = sp.fraction(expr)
    total = 1
    rep = []
    for part, ispart_den in ((num, False), (den, True)):
        cst, facs = sp.factor_list(sp.expand(part), *free) if free else (part, [])
        if cst < 0:
            total = -total
        elif cst == 0:
            return '0', []
        for f, m in facs:
            s, how = poly_sign_cert(f, free)
            if s is None:
                return None, rep + [(f, m, None, how)]
            if s == '0':
                if ispart_den:
                    return None, rep + [(f, m, '0', 'zero denominator')]
                return '0', []
            if s == '-' and m % 2 == 1:
                total = -total
            rep.append((f, m, s, how))
    return ('+' if total > 0 else '-'), rep


# ------------------------------------------------------------------ graph structure of M0
def support(M0, free):
    """0/1 support of M0 on the open stratum; asserts every entry is 0 or manifestly positive
    (a product of factors x, 1-x, positive constants)."""
    S = [[0] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            e = sp.factor(M0[i, j])
            if e == 0:
                continue
            s, rep = factor_sign_cert(e, free)
            assert s == '+', (M0[i, j], rep)
            S[i][j] = 1
    return S


def closed_classes(S):
    n = len(S)
    reach = [[bool(S[i][j]) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    seen, classes = set(), []
    for i in range(n):
        if i in seen:
            continue
        cls = [j for j in range(n) if reach[i][j] and reach[j][i]]
        seen.update(cls)
        closed = all((not reach[x][y]) or (y in cls) for x in cls for y in range(n))
        if closed:
            classes.append(cls)
    return classes


def stationary_on(P, idx):
    """Stationary law of the irreducible stochastic block P[idx, idx] (sympy)."""
    n = len(idx)
    if n == 1:
        return {idx[0]: sp.Integer(1)}
    xs = sp.symbols('m0:%d' % n)
    eqs = [sp.Eq(sum(xs[i] * P[idx[i], idx[j]] for i in range(n)), xs[j]) for j in range(n - 1)]
    eqs.append(sp.Eq(sum(xs), 1))
    sol = sp.solve(eqs, xs, dict=True)
    assert len(sol) == 1
    return {idx[i]: sp.factor(sp.cancel(sol[0][xs[i]])) for i in range(n)}


def absorption(M0, classes):
    """h[j][s]: probability, under M0, of absorption in classes[j] from state s."""
    rec = set(x for cl in classes for x in cl)
    T = [s for s in range(4) if s not in rec]
    h = []
    for j, cl in enumerate(classes):
        hj = {s: sp.Integer(1 if s in cl else 0) for s in rec}
        if T:
            ys = sp.symbols('y0:%d' % len(T))
            eqs = []
            for i, s in enumerate(T):
                rhs = sum(M0[s, t] for t in cl) + sum(M0[s, T[k]] * ys[k] for k in range(len(T)))
                eqs.append(sp.Eq(ys[i], rhs))
            sol = sp.solve(eqs, ys, dict=True)
            assert len(sol) == 1
            for i, s in enumerate(T):
                hj[s] = sp.factor(sp.cancel(sol[0][ys[i]]))
        h.append(hj)
    return h, T


def limit_W(sig, p, free, log=None):
    """Returns (W as list of 4 sympy expressions, info dict with all objects to certify)."""
    M0, M1, M2 = expansion(sig, p)
    S = support(M0, free)
    classes = closed_classes(S)
    mus = [stationary_on(M0, cl) for cl in classes]
    info = {'classes': classes, 'mus': mus, 'S': S}
    to_certify_pos = []  # expressions that must be > 0 (or, for '*-det'/'h-den', != 0)
    for cl, mu in zip(classes, mus):
        # identity check: mu is a law and stationary for the block, as rational functions
        assert sp.simplify(sum(mu.values()) - 1) == 0
        for j in cl:
            assert sp.simplify(sum(mu[i] * M0[i, j] for i in cl) - mu[j]) == 0
        for s, val in mu.items():
            to_certify_pos.append(('mu', val))
        if len(cl) > 1:
            # the 'normaliser' of this step: det of the stationary system of the block
            n = len(cl)
            Ms = sp.Matrix(n, n, lambda r, q: (1 if q == r else 0) - M0[cl[q], cl[r]] if r < n - 1 else 1)
            to_certify_pos.append(('stat-det', sp.factor(Ms.det())))
    if len(classes) == 1:
        W = [mus[0].get(s, sp.Integer(0)) for s in range(4)]
        info.update(k=1, Q=None, w=[1], h=None, T=None, pos=to_certify_pos, zero_or_pos=[])
        return W, info
    h, T = absorption(M0, classes)
    if T:
        to_certify_pos.append(('trans-det', sp.factor((sp.eye(len(T)) - M0.extract(T, T)).det())))
    for hj in h:
        # identity check: h_j is M0-harmonic
        for s in range(4):
            assert sp.simplify(sum(M0[s, t] * hj[t] for t in range(4)) - hj[s]) == 0
        for s in T:
            d = sp.fraction(sp.together(hj[s]))[1]
            to_certify_pos.append(('h-den', d))
    assert all(sp.simplify(sum(hj[s] for hj in h) - 1) == 0 for s in range(4))
    k = len(classes)
    Q = sp.zeros(k, k)
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            Q[i, j] = sp.factor(sp.cancel(sum(mus[i][s] * sum(M1[s, t] * h[j][t] for t in range(4))
                                              for s in classes[i])))
        Q[i, i] = -sum(Q[i, j] for j in range(k) if j != i)
    # rows of Q must sum to zero identically (consistency of the lemma)
    zero_or_pos = []
    SQ = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                s, rep = factor_sign_cert(Q[i, j], free)
                assert s in ('0', '+'), ('Q entry not certified', Q[i, j], rep)
                SQ[i][j] = 1 if s == '+' else 0
                zero_or_pos.append(('Q', Q[i, j]))
    qclasses = closed_classes(SQ)
    if len(qclasses) != 1:
        raise RuntimeError('first-order reduced chain has %d closed classes: second order needed'
                           % len(qclasses))
    D = qclasses[0]
    if len(D) == 1:
        w = [sp.Integer(1) if i == D[0] else sp.Integer(0) for i in range(k)]
    else:
        # stationary law of Q restricted to D: solve w Q_D = 0, sum w = 1
        xs = sp.symbols('w0:%d' % len(D))
        eqs = [sp.Eq(sum(xs[a_] * Q[D[a_], D[b_]] for a_ in range(len(D))), 0) for b_ in range(len(D) - 1)]
        eqs.append(sp.Eq(sum(xs), 1))
        sol = sp.solve(eqs, xs, dict=True)
        assert len(sol) == 1
        w = [sp.Integer(0)] * k
        for a_, i in enumerate(D):
            w[i] = sp.factor(sp.cancel(sol[0][xs[a_]]))
            to_certify_pos.append(('w', w[i]))
    W = [sp.factor(sp.cancel(sum(w[i] * mus[i].get(s, 0) for i in range(k)))) for s in range(4)]
    info.update(k=k, Q=Q, w=w, h=h, T=T, qclasses=qclasses, pos=to_certify_pos, zero_or_pos=zero_or_pos)
    return W, info
