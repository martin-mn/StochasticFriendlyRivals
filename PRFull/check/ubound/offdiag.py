"""wedge W: the non-excluded class-C strata whose (p6, p8) is not (0,0) or (1,1) -- a by-product for
the separation question (referee issue 2).  For each, an exhaustive single-co-player exclusion search
(check/cert/certify.py exclude_stratum, search_all: every q whose centre margin is certified negative
or unresolved)."""
import sys, time, json, numpy as np
sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/check')
from cert import engine as en, certify as ce
import status as st
idx, qp, qm = st.load(sys.argv[1]); M, d = st.masks(idx); m, sg, b = M['W']
left = [int(idx[t]) for t in np.where(m & (qp < 0))[0]]
off = [i for i in left if (st.digit(i, 6), st.digit(i, 8)) not in ((0, 0), (1, 1))]
off.sort(key=lambda i: -st.pattern(i).count('*'))
lim = int(sys.argv[2]) if len(sys.argv) > 2 else len(off)
out = []; t0 = time.time()
for i in off[:lim]:
    p = st.pattern(i); r = ce.exclude_stratum(p, wedge='W', search_all=True, nthreads=3)
    out.append(dict(pattern=p, d=p.count('*'), excluded=bool(r['excluded']), q=r['q'], n_centre_negative=r.get('n_centre_negative'),
                    tried=r.get('n_candidates_tried'), seconds=round(r['seconds'], 2), timed_out=r.get('timed_out', False)))
    print(out[-1], flush=True)
json.dump(out, open('runs/offdiag.json', 'w'), indent=1)
print('%d off-diagonal strata: excluded %d, not excludable by any single co-player %d  (%.0f s)'
      % (len(out), sum(o['excluded'] for o in out), sum(not o['excluded'] for o in out), time.time() - t0))
