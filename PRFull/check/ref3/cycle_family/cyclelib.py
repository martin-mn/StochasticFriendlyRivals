"""Exact eps->0 limits for the family p = (0, a, b, 1) against the 16 binary memory-one deviators.

Conventions (as in ms.tex): states CC=0, CD=1, DC=2, DD=3, focal (deviator) first; the co-player
sees CD and DC exchanged; errors act as x -> eps + (1-2 eps) x; payoffs (R,S,T,P) = (1,u,1+v,0).

The limit is taken by the Markov-chain-tree theorem: the unnormalised weights d_X are the
principal 3x3 cofactors of I - P, polynomials in eps (and the free coordinates).  On a stratum,
with k the lowest eps-order of D = sum d_X whose coefficient is not identically zero, the limit is
W_X = [eps^k] d_X / [eps^k] D at every point of the stratum PROVIDED (i) every coefficient of
eps^j, j < k, of every d_X vanishes identically on the stratum and (ii) [eps^k] D has no zero on
the open stratum.  Both are checked; (ii) by a positivity certificate (see certify_positive).
"""
import itertools
import sympy as sp

e, a, b, u, v = sp.symbols('e a b u v')
s_, t_ = sp.symbols('s t', positive=True)
M = [0, 2, 1, 3]  # CD <-> DC
BIN = [tuple(bits) for bits in itertools.product([0, 1], repeat=4)]
NAME = {(0, 0, 0, 0): 'ALLD', (1, 1, 1, 1): 'ALLC', (1, 0, 1, 0): 'TFT', (1, 0, 0, 1): 'WSLS',
        (0, 1, 0, 1): 'ATFT'}


def label(sig):
    return ''.join(str(x) for x in sig)


def clamp(x):
    return e + (1 - 2 * e) * x


def tree_weights(sig, p):
    P = sp.zeros(4, 4)
    for st in range(4):
        x = clamp(sig[st]); y = clamp(p[M[st]])
        P[st, 0] = x * y; P[st, 1] = x * (1 - y); P[st, 2] = (1 - x) * y; P[st, 3] = (1 - x) * (1 - y)
    A = sp.eye(4) - P
    out = []
    for i in range(4):
        idx = [k for k in range(4) if k != i]
        out.append(sp.expand(A.extract(idx, idx).det()))
    return out


def certify_positive(f, free):
    """Certify f > 0 on the open box (0,1)^len(free).  Returns a short string naming the
    certificate, or raises.  Method: x = s/(1+s) maps (0,1) onto (0,inf); after clearing the
    positive denominator (1+s)^deg, a polynomial with all coefficients >= 0 (not all zero) is
    > 0 for s > 0.  If that fails, factor f and certify each factor's sign the same way."""
    f = sp.expand(f)
    if not free:
        if f > 0:
            return 'constant %s' % f
        raise ValueError('nonpositive constant %s' % f)

    def onefactor(g):
        g = sp.expand(g)
        subs = {}
        mult = 1
        for x, y in zip(free, (s_, t_)):
            d = sp.degree(g, x)
            subs[x] = y / (1 + y)
            mult *= (1 + y) ** d
        h = sp.expand(sp.cancel(g.subs(subs) * mult))
        if h == 0:
            return 0
        cs = sp.Poly(h, *[y for _, y in zip(free, (s_, t_))]).coeffs()
        if all(c >= 0 for c in cs):
            return 1
        if all(c <= 0 for c in cs):
            return -1
        return None

    sgn = onefactor(f)
    if sgn == 1:
        return 'all coefficients >= 0 after x=s/(1+s)'
    c, facs = sp.factor_list(f)
    total = 1 if c > 0 else -1
    for g, m in facs:
        sg = onefactor(g)
        if sg is None or sg == 0:
            raise ValueError('cannot certify sign of factor %s of %s' % (g, f))
        total *= sg ** m
    if total != 1:
        raise ValueError('f is negative: %s' % f)
    return 'factorwise: ' + ' * '.join('(%s)^%d' % (g, m) for g, m in facs) + ' const %s' % c


def limit_on_stratum(sig, p, free):
    """Return (W, k, cert): exact limit valid on the whole open stratum."""
    ds = tree_weights(sig, p)
    D = sp.expand(sum(ds))
    PD = sp.Poly(D, e)
    coeffs = {mon[0]: sp.expand(c) for mon, c in zip(PD.monoms(), PD.coeffs())}
    k = min(j for j, c in coeffs.items() if c != 0)
    cert = certify_positive(coeffs[k], free)
    W = []
    for d in ds:
        Pd = sp.Poly(d, e)
        cd = {mon[0]: sp.expand(c) for mon, c in zip(Pd.monoms(), Pd.coeffs())}
        for j, c in cd.items():
            if j < k and c != 0:
                raise RuntimeError('lower-order numerator term on stratum')
        W.append(sp.factor(sp.cancel(cd.get(k, 0) / coeffs[k])))
    assert sp.simplify(sum(W) - 1) == 0
    return W, k, cert


def payoff(W):
    return W[0] + u * W[1] + (1 + v) * W[2]


STRATA = {
    'square 0<a,b<1': ((0, a, b, 1), [a, b]),
    'edge a=0 (0<b<1)': ((0, 0, b, 1), [b]),
    'edge a=1 (0<b<1)': ((0, 1, b, 1), [b]),
    'edge b=0 (0<a<1)': ((0, a, 0, 1), [a]),
    'edge b=1 (0<a<1)': ((0, a, 1, 1), [a]),
    'corner (0,0,1,1)': ((0, 0, 1, 1), []),
    'corner (0,1,0,1)': ((0, 1, 0, 1), []),
    'corner (0,0,0,1)': ((0, 0, 0, 1), []),
    'corner (0,1,1,1)': ((0, 1, 1, 1), []),
}
