"""non-excluded strata that the census file st_eff.txt does not list"""
import sys, numpy as np, status as st
idx, qp, qm = st.load(sys.argv[1]); M, d = st.masks(idx)
def conv(g):
    r = 0
    for j in range(16):
        dg = g % 3; g //= 3; r += {0: 0, 1: 2, 2: 1}[dg] * 3**j
    return r
C = {'W': set(), 'N': set(), 'E': set(), 'ray+': set()}
for l in open('/Users/martin/Documents/CL1/PartnersRivalsM2f/check/st_eff.txt'):
    if not l.strip(): continue
    g, c, nP, nM, K = l.split(); i = conv(int(g))
    if c == 'C' and int(nP) > 0: C['W'].add(i)
    if c == 'A' and int(nP) > 0: C['N'].add(i)
    if c == 'A' and int(nM) > 0: C['E'].add(i)
    if int(nP) > 0: C['ray+'].add(i)
def mir(p): return ''.join({'0': '1', '1': '0', '*': '*'}[p[15 - s]] for s in range(16))
Epat = set(st.pattern(i) for i in C['E']); Npat = set(st.pattern(i) for i in C['N'])
for k in ('W', 'N', 'E', 'ray+'):
    m, sg, b = M[k]; q = qp if sg > 0 else qm
    ne = set(int(idx[t]) for t in np.where(m & (q < 0))[0]); x = sorted(ne - C[k])
    print(k, 'non-excluded', len(ne), 'census', len(C[k]), 'census subset:', C[k] <= ne, 'extra', len(x))
    for i in x:
        p = st.pattern(i); note = ''
        if k == 'N': note = 'mirror in census E: %s' % (mir(p) in Epat)
        if k == 'E': note = 'mirror in census N: %s' % (mir(p) in Npat)
        print('   ', p, 'dim', p.count('*'), note)
