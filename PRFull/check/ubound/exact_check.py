"""exact (rational, gmpy2) cross-checks of the upper-bound certificate on random samples:
 (1) binary pairs (corner c, certificate co-player q): the row entry used by ub (E_c and the sign
     of w_DC - w_CD from the interval engine) against the exact limit;
 (2) target strata: the exact margin against the certificate q at the centre and at a random
     rational interior point is < 0 (what the exclusion asserts for EVERY interior point);
 (3) efficiency classes: exact self-play support at a random rational interior point vs cls.bin."""
import sys, os, subprocess, time, numpy as np
from fractions import Fraction
sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/check')
from cert import exact as ex
import status as st
prefix = sys.argv[1]; NP = int(sys.argv[2]); NS_ = int(sys.argv[3]); NC = int(sys.argv[4])
rng = np.random.default_rng(2026)
T = np.fromfile(prefix + '_tri.bin', dtype=np.int64).reshape(-1, 3)
cls, dim = st.base()
def digits(i):
    d = []
    for _ in range(16): d.append(i % 3); i //= 3
    return d
def randpoint(i, den=None):
    return [Fraction(int(rng.integers(1, den)), den) if x == 2 else Fraction(x) for x in digits(i)] if den else \
           [Fraction(1, 2) if x == 2 else Fraction(x) for x in digits(i)]
t0 = time.time()
# (1)
sel = T[rng.integers(0, len(T), NP)]
qs = sorted(set(int(q) for q in sel[:, 2])); rows = {}
os.makedirs('runs/rows', exist_ok=True)
for q in qs:
    fn = 'runs/rows/%d.bin' % q
    if not os.path.exists(fn): subprocess.run(['nice', '-n', '10', './ub', 'row', str(q), fn], check=True, capture_output=True)
    rows[q] = np.fromfile(fn, dtype=np.uint8)
bad1 = 0; hist = {}
for idx, sg, q in sel:
    d = digits(int(idx)); c = [int(rng.integers(0, 2)) if x == 2 else x for x in d]
    g = sum(b << j for j, b in enumerate(c))
    r = int(rows[int(q)][g]); E = r >> 2; f = r & 3
    L = ex.lead_exact([Fraction(b) for b in c], ex.g2p(int(q)))
    m = L['W'][2] - L['W'][1]
    fx = 2 if m > 0 else (1 if m < 0 else 0)
    hist[f] = hist.get(f, 0) + 1
    if E != L['E'] or f != fx: bad1 += 1; print('ROW MISMATCH', idx, q, g, E, f, L['E'], m)
print('(1) %d binary pairs (%d certificate co-players): row entry vs exact: %d mismatches; flags seen %s  [%.0f s]'
      % (NP, len(qs), bad1, hist, time.time() - t0), flush=True)
# (2)
t1 = time.time(); sel = T[rng.permutation(len(T))[:NS_]]; bad2 = 0; mx = None
for idx, sg, q in sel:
    for den in (None, 97):
        p = randpoint(int(idx), den); v = ex.margin(p, ex.g2p(int(q)), int(sg))
        if v >= 0: bad2 += 1; print('INTERIOR MARGIN NOT NEGATIVE', st.pattern(int(idx)), sg, q, den, v)
        mx = v if mx is None or v > mx else mx
print('(2) %d target strata x 2 interior points (centre, random k/97): exact margin < 0 at all but %d; largest %s = %.3g  [%.0f s]'
      % (NS_, bad2, mx, float(mx), time.time() - t1), flush=True)
# (3)
t2 = time.time(); bad3 = 0; cnt = {}
tg = T[rng.permutation(len(T))[:NC // 2], 0]; rnd = rng.integers(0, 3**16, NC - NC // 2)
for idx in list(tg) + list(rnd):
    p = randpoint(int(idx), 89); L = ex.lead_exact(p, p); s = L['supp']
    r = [0, 0, 0, 0]
    for k in range(16):
        if s >> k & 1: r[k >> 2] = 1
    c = 1 if not (r[1] or r[2] or r[3]) else (2 if not (r[0] or r[3]) else (3 if not r[3] else 0))
    cnt[c] = cnt.get(c, 0) + 1
    if c != cls[int(idx)]: bad3 += 1; print('CLASS MISMATCH', st.pattern(int(idx)), c, cls[int(idx)])
print('(3) %d strata (%d targets, %d uniform): exact self-play class at a random k/89 point vs cls.bin: %d mismatches; classes seen %s  [%.0f s]'
      % (NC, NC // 2, NC - NC // 2, bad3, cnt, time.time() - t2))
