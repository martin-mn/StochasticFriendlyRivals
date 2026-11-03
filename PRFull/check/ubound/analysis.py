"""what is left after the exclusion: the non-excluded efficient strata at the top dimension of each
locus, their relation to the census and to the certified witnesses, and (wedge W) the split by
(p6, p8) = (p_{CD,DC}, p_{DC,CC})."""
import sys, numpy as np, status as st
prefix = sys.argv[1]
idx, qp, qm = st.load(prefix); M, d = st.masks(idx)
def conv(g):
    r = 0
    for j in range(16):
        dg = g % 3; g //= 3; r += {0: 0, 1: 2, 2: 1}[dg] * 3**j
    return r
CEN = {'W': set(), 'N': set(), 'E': set(), 'ray+': set()}
for l in open('/Users/martin/Documents/CL1/PartnersRivalsM2f/check/st_eff.txt'):
    if not l.strip(): continue
    g, c, nP, nM, K = l.split(); i = conv(int(g))
    if c == 'C' and int(nP) > 0: CEN['W'].add(i)
    if c == 'A' and int(nP) > 0: CEN['N'].add(i)
    if c == 'A' and int(nM) > 0: CEN['E'].add(i)
    if int(nP) > 0: CEN['ray+'].add(i)
CEN['W68'] = set(i for i in CEN['W'] if st.digit(i, 6) == 0 and st.digit(i, 8) == 0)
WIT = {'1*1*001*1*1**0*0': 'W witness, certified 7-dim box', '*10***00*1***010': 'N witness, certified 8-dim box',
       '1*1*0*0101*1*010': 'certified whole open 5-cube'}
for k in ('W', 'N', 'E', 'ray+', 'W68'):
    m, sg, b = M[k]; q = qp if sg > 0 else qm
    left = np.where(m & (q < 0))[0]; top = d[left].max()
    tl = [int(idx[t]) for t in left if d[t] == top]
    notex = set(int(idx[t]) for t in left)
    print('%-5s non-excluded efficient strata: %d (census FR strata %d, all among them: %s); top dimension %d: %d strata'
          % (k, len(notex), len(CEN[k]), CEN[k] <= notex, top, len(tl)))
    for i in tl:
        p = st.pattern(i); print('      %s  census FR: %-5s %s' % (p, i in CEN[k], WIT.get(p, '')))
# wedge W by (p6, p8)
m, sg, b = M['W']; left = np.where(m & (qp < 0))[0]
print('\nwedge W, non-excluded class-C strata by (p6, p8) and dimension   [census FR strata in brackets]')
lab = '01*'
for a in range(3):
    for c in range(3):
        s = [t for t in left if st.digit(int(idx[t]), 6) == a and st.digit(int(idx[t]), 8) == c]
        if not s: continue
        h = np.bincount(d[s], minlength=8)[:8]
        hc = np.bincount([d[t] for t in s if int(idx[t]) in CEN['W']], minlength=8)[:8]
        print('  p6=%s p8=%s: %4d strata, by dim %s   [%d: %s]' % (lab[a], lab[c], len(s), h.tolist(), int(hc.sum()), hc.tolist()))
