"""Limit rows of the 16 binary deviators against a memory-one resident p, eps->0, symbolic.
Usage: python3 m1rows.py  (prints, for the open facet p=(1,a,b,c), each deviator's limit
frequencies W(sigma,p) and its stability row  d_CD*(u-1) + d_DC*v - d_DD  (<=0 means no gain)).
States CC=0, CD=1, DC=2, DD=3, focal (deviator) first."""
import itertools, sys
import sympy as sp

e = sp.symbols('e', positive=True)
M = [0, 2, 1, 3]  # swap CD <-> DC


def clamp(x):
    return e + (1 - 2 * e) * x


def stationary_limit(sig, p):
    """sig, p: 4-tuples (sympy). Returns limit W as list of 4 sympy exprs."""
    P = sp.zeros(4, 4)
    for s in range(4):
        x = clamp(sig[s]); y = clamp(p[M[s]])
        P[s, 0] = x * y; P[s, 1] = x * (1 - y); P[s, 2] = (1 - x) * y; P[s, 3] = (1 - x) * (1 - y)
    # Markov chain tree theorem via cofactors of (I - P): w_i proportional to det of (I-P) minus row/col i
    A = sp.eye(4) - P
    N = []
    for i in range(4):
        idx = [k for k in range(4) if k != i]
        N.append(sp.expand(A.extract(idx, idx).det()))
    tot = sp.expand(sum(N))
    polys = [sp.Poly(n, e) for n in N]
    ptot = sp.Poly(tot, e)
    # lowest order with non-zero coefficient in tot
    terms = sorted(ptot.terms(), key=lambda t: t[0][0])
    for (k,), coeff in terms:
        if sp.simplify(coeff) != 0:
            kmin = k; cden = sp.factor(coeff); break
    W = []
    for pl in polys:
        c = pl.coeff_monomial(e ** kmin) if kmin > 0 else pl.coeff_monomial(1)
        # lower orders in numerators must vanish
        for (k,), cc in pl.terms():
            if k < kmin and sp.simplify(cc) != 0:
                raise RuntimeError('numerator has lower order')
        W.append(sp.factor(sp.simplify(c / cden)))
    return W


if __name__ == '__main__':
    a, b, c = sp.symbols('a b c', positive=True)
    u, v = sp.symbols('u v')
    p = (sp.Integer(1), a, b, c)
    Wpp = stationary_limit(p, p)
    print('W(p,p) =', Wpp)
    for bits in itertools.product([0, 1], repeat=4):
        sig = tuple(sp.Integer(x) for x in bits)
        W = stationary_limit(sig, p)
        d = [sp.simplify(W[k] - Wpp[k]) for k in range(4)]
        row = sp.factor(sp.simplify(d[1] * (u - 1) + d[2] * v - d[3]))
        print(''.join(map(str, bits)), 'W =', [sp.factor(w) for w in W], ' row =', row)
        sys.stdout.flush()
