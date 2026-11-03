#!/usr/bin/env python3
"""facets_cert.py -- what the dimension of the memory-one stable set needs beyond the facet p_CC = 1 (SI 'The dimension
of the stable set'), computed with the exact routines of the boundary-strata prover (check/ref3/eq3_boundary/P1).

    python3 facets_cert.py            (about 1 minute; one core)

A three-dimensional stable set must meet the interior or one of the eight open facets in a three-dimensional set,
since the other strata have dimension at most 2.  The facet p_CC = 1 is settled by the display of SI 'Partners'.
This script does the other seven open facets:

 F2  p = (a,b,c,0), p_DD = 0.  Self-play limit (0,0,0,1).  The 16 gains as l1 I* + l2 II* - mu with
     I* = (1-a)u + c, II* = (1-b)u + c(1+v) (the SI's two conditions), l1, l2, mu >= 0 certified; ALLC gains
     I*/(1-a+c) and 0001 gains II*/(...), so stability <=> I* <= 0, II* <= 0, and no stable point if u >= 0.
 F3..F8  p_CC = 0, p_DD = 1, p_CD = 0, p_CD = 1, p_DC = 0, p_DC = 1.  (i) The self-play limit is positive on all
     four states at every point of the open facet (each W_X a rational function certified positive), so the three
     free states bind; (ii) at every game (u,v) some pair of deviators has a payoff difference that is not
     identically 0 on the facet: the numerators N_sigma(a,b,c;u,v) of pi(sigma,p) - pi(ALLD,p) are affine in (u,v),
     and the linear system 'every coefficient of every N_sigma vanishes' has no solution (u,v).  So the equalisers of
     the open facet lie, at every game, in the zero set of a non-zero polynomial: dimension <= 2.
 A   the alternators (x,0,1,y), 0 < x, y < 1, and the edges (0,0,1,y), (x,0,1,1) (lower bound 2 in the quadrant): self-play limit
     (0,1/2,1/2,0) and the 16 gains as l1 C1 + l2 C2 - mu with C1 = (1+v-u)y - (1+u+v), C2 = (1+v-u)x - 2v,
     l1, l2, mu >= 0 certified.

Every limit is the prover's limit_W: lower eps-orders of the numerators vanish identically and the leading
coefficient of the normaliser is certified positive on the open stratum, so each limit holds at every point.
"""
import itertools
import os
import sys
import time

import sympy as sp

sys.dont_write_bytecode = True
# the prover's directory: ../P1 when this script sits in check/ref3/eq3_boundary/<dir>/, else the project path
_P1_CAND = [os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'P1'),
            '/Users/martin/Documents/CL1/PRFull/check/ref3/eq3_boundary/P1']
P1 = next(d for d in _P1_CAND if os.path.isfile(os.path.join(d, 'prove_strata.py')))
sys.path.insert(0, P1)
import prove_strata as ps                                                    # noqa: E402

R, A, B, C = ps.R, ps.A, ps.B, ps.C
a, b, c = ps.a_s, ps.b_s, ps.c_s
u, v = ps.u, ps.v
FREE = [a, b, c]
SIGMAS = ps.SIGMAS
one, zero = R(1), R(0)


def pay(W):
    return W[0] + u * W[1] + (1 + v) * W[2]          # deviator's payoff, (R,S,T,P) = (1,u,1+v,0)


def limits(p, free):
    out = {}
    for s in SIGMAS:
        W, inf = ps.limit_W(ps.sig_tuple(s), p, free)
        assert inf['lead_ok'], (s, inf)
        out[s] = W
    Wpp, inf = ps.limit_W(p, p, free)
    assert inf['lead_ok']
    return Wpp, out


def facet_F2():
    p = (A, B, C, zero)
    Wpp, Ws = limits(p, FREE)
    assert [sp.simplify(w) for w in Wpp] == [0, 0, 0, 1]
    Is = (1 - a) * u + c
    IIs = (1 - b) * u + c * (1 + v)
    allok = True
    rows = {}
    for s in SIGMAS:
        g = sp.expand(sp.together(pay(Ws[s])))        # pi(p,p) = P = 0
        num, den = sp.fraction(sp.together(g))
        # g = A0 + B0 u + C0 v as rational functions of (a,b,c)
        A0 = sp.cancel(sp.together(g).subs({u: 0, v: 0}))
        B0 = sp.cancel(sp.diff(g, u))
        C0 = sp.cancel(sp.diff(g, v))
        l2 = sp.cancel(C0 / c)
        l1 = sp.cancel((B0 - l2 * (1 - b)) / (1 - a))
        mu = sp.cancel(l1 * c + l2 * c - A0)
        ident = sp.simplify(g - (l1 * Is + l2 * IIs - mu)) == 0
        oks = [ps.nonneg(q, FREE)[0] for q in (l1, l2, mu)]
        ok = ident and all(oks)
        allok = allok and ok
        rows[s] = (sp.factor(l1), sp.factor(l2), sp.factor(mu), ok)
    print('F2 p=(a,b,c,0): self-play (0,0,0,1); every gain = l1 I* + l2 II* - mu, l1,l2,mu >= 0 certified: %s' % allok)
    for s in SIGMAS:
        l1, l2, mu, ok = rows[s]
        print('   %s %-5s l1=%s  l2=%s  mu=%s  %s' % (s, ps.NAMES.get(s, ''), l1, l2, mu, ok))
    return allok


def facet_binding(name, p, free):
    Wpp, Ws = limits(p, free)
    pos = [ps.positive(w, free)[0] for w in Wpp]
    # payoff differences against ALLD
    P0 = pay(Ws['0000'])
    eqs = []
    for s in SIGMAS:
        if s == '0000':
            continue
        d = sp.together(pay(Ws[s]) - P0)
        num, _ = sp.fraction(d)
        num = sp.expand(num)
        if num == 0:
            continue
        poly = sp.Poly(num, *free)
        eqs.extend(poly.coeffs())                     # each affine in (u, v)
    # every coefficient vanishes <=> (u,v) solves eqs
    sol = sp.solve([sp.expand(e) for e in eqs], [u, v], dict=True)
    # also: is pi(sigma,p) - pi(p,p) identically 0 for all sigma at some game?  (the same system plus self-play)
    print('%s p=%s: self-play limit positive on all four states: %s; W(p,p) = %s' % (name, tuple(str(x) for x in p), all(pos), [sp.factor(w) for w in Wpp]))
    print('     games at which all 16 deviators tie identically on the open facet: %s  (%d coefficient equations)' % (sol, len(eqs)))
    return all(pos) and sol == []


def alternators():
    """The open alternating face (x,0,1,y) and its two edges that are still alternators, (0,0,1,y) and (x,0,1,1)."""
    cases = (('(x,0,1,y)', (A, zero, one, B), [a, b], a, b),
             ('(0,0,1,y)', (zero, zero, one, B), [b], sp.Integer(0), b),
             ('(x,0,1,1)', (A, zero, one, one), [a], a, sp.Integer(1)))
    allok_all = True
    for name, p, free, x, y in cases:
        Wpp, Ws = limits(p, free)
        assert [sp.simplify(w) for w in Wpp] == [0, sp.Rational(1, 2), sp.Rational(1, 2), 0], Wpp
        Ppp = (u + 1 + v) / 2
        C1 = (1 + v - u) * y - (1 + u + v)
        C2 = (1 + v - u) * x - 2 * v
        allok = True
        rows = {}
        for s in SIGMAS:
            g = sp.together(pay(Ws[s]) - Ppp)
            A0 = sp.cancel(g.subs({u: 0, v: 0}))
            B0 = sp.cancel(sp.diff(g, u))
            C0 = sp.cancel(sp.diff(g, v))
            # B0 = -(1+y) l1 - x l2 ;  C0 = (y-1) l1 + (x-2) l2   (determinant 2(1-x+y) > 0)
            l1s, l2s = sp.symbols('l1s l2s')
            sol = sp.solve([sp.Eq(B0, -(1 + y) * l1s - x * l2s), sp.Eq(C0, (y - 1) * l1s + (x - 2) * l2s)],
                           [l1s, l2s], dict=True)[0]
            l1, l2 = sp.cancel(sol[l1s]), sp.cancel(sol[l2s])
            mu = sp.cancel(l1 * (y - 1) + l2 * x - A0)
            ident = sp.simplify(g - (l1 * C1 + l2 * C2 - mu)) == 0
            oks = [ps.nonneg(q, free)[0] for q in (l1, l2, mu)]
            ok = ident and all(oks)
            allok = allok and ok
            rows[s] = (sp.factor(l1), sp.factor(l2), sp.factor(mu), ok)
        print('A  alternators %s: self-play (0,1/2,1/2,0); every gain = l1 C1 + l2 C2 - mu, l1,l2,mu >= 0 certified: %s'
              % (name, allok))
        if name == '(x,0,1,y)':
            for s in SIGMAS:
                l1, l2, mu, ok = rows[s]
                print('   %s %-5s l1=%s  l2=%s  mu=%s  %s' % (s, ps.NAMES.get(s, ''), l1, l2, mu, ok))
        allok_all = allok_all and allok
    return allok_all


def main():
    t0 = time.time()
    ok2 = facet_F2()
    facets = [('F3 p_CC=0', (zero, A, B, C)), ('F4 p_DD=1', (A, B, C, one)),
              ('F5 p_CD=0', (A, zero, B, C)), ('F6 p_CD=1', (A, one, B, C)),
              ('F7 p_DC=0', (A, B, zero, C)), ('F8 p_DC=1', (A, B, one, C))]
    okb = [facet_binding(n, p, FREE) for n, p in facets]
    okA = alternators()
    print('SUMMARY  F2 reduction: %s;  F3-F8 all states recur and no game makes all deviators tie identically: %s;'
          '  alternator reduction: %s   (%.1f s)' % (ok2, okb, okA, time.time() - t0))


if __name__ == '__main__':
    main()
