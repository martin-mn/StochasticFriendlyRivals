"""consistency with the float census (PartnersRivalsM2f/check/st_eff.txt, read-only): every stratum in
which the census confirmed a friendly rival must (i) have the same exact efficiency class here and
(ii) NOT be excluded by the cover.  Also: the non-excluded efficient strata versus the census set."""
import sys, numpy as np, status as st
prefix = sys.argv[1]
idx, qp, qm = st.load(prefix); cls, dim = st.base()
pos = {int(i): t for t, i in enumerate(idx)}
CEN = {'-': 0, 'C': 1, 'A': 2, 'Z': 3}
def conv(g):        # census digits 0 -> 0, 1 -> free, 2 -> one  ->  ours 0 -> 0, 1 -> one, 2 -> free
    r = 0
    for j in range(16):
        dg = g % 3; g //= 3; r += {0: 0, 1: 2, 2: 1}[dg] * 3**j
    return r
lines = open('/Users/martin/Documents/CL1/PartnersRivalsM2f/check/st_eff.txt').read().split('\n')
mism = 0; excl = {'W': [], 'N': [], 'E': [], 'ray+': [], 'S': [], 'ray-': []}; cen = {k: set() for k in excl}
for l in lines:
    if not l.strip(): continue
    g, c, nP, nM, K = l.split(); i = conv(int(g)); nP = int(nP); nM = int(nM)
    if cls[i] != CEN[c]: mism += 1; continue
    t = pos[i]
    for k, cs, sg, n in (('W', {1}, +1, nP), ('S', {1}, -1, nM), ('N', {2}, +1, nP), ('E', {2}, -1, nM),
                          ('ray+', {1, 2, 3}, +1, nP), ('ray-', {1, 2, 3}, -1, nM)):
        # the census counts nP/nM for the stratum's own class; ray loci take every efficient class
        if cls[i] in cs and n > 0:
            cen[k].add(i)
            if (qp if sg > 0 else qm)[t] >= 0: excl[k].append(i)
print('census lines', len(lines) - 1, '; class mismatches with the exact classes here:', mism)
M, d = st.masks(idx)
for k in ('W', 'N', 'E', 'ray+', 'S', 'ray-'):
    m, sg, b = M[k]; q = qp if sg > 0 else qm
    notex = set(int(x) for x in idx[m & (q < 0)])
    c = cen[k]
    print('%-5s census FR strata %8d | of them excluded by the cover: %d | non-excluded %8d | census subset of non-excluded: %s | non-excluded not in census: %d'
          % (k, len(c), len(excl[k]), len(notex), c <= notex, len(notex - c)))
    if excl[k][:5]: print('   EXCLUDED CENSUS STRATA (contradiction):', [st.pattern(i) for i in excl[k][:10]])
