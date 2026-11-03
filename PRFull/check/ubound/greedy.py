"""adaptive cover: pick not-yet-excluded target strata, propose co-players whose float margin at the
stratum centre is most negative, and run `ub cover` on the new ones (the proof is in ub.c; the float
scan only proposes candidates)."""
import sys, os, time, subprocess, numpy as np
sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/check')
from cert import engine as en
import status as st
PREFIX = sys.argv[1]; ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 10
PER = int(sys.argv[3]) if len(sys.argv) > 3 else 40          # strata sampled per locus per round
TOP = int(sys.argv[4]) if len(sys.argv) > 4 else 3            # candidates per stratum
rng = np.random.default_rng(11)
def mirror_g(g):
    return sum((1 - ((g >> (15 - s)) & 1)) << s for s in range(16))
tried = set()
logf = PREFIX + '.log'
if os.path.exists(logf):
    for l in open(logf): tried.add(int(l.split()[0]))
def run(qs):
    qs = [q for q in dict.fromkeys(qs) if q not in tried]
    if not qs: return 0
    fn = PREFIX + '.qnew'; open(fn, 'w').write('\n'.join(map(str, qs)) + '\n')
    r = subprocess.run(['nice', '-n', '10', './ub', 'cover', 'cls.bin', fn, PREFIX, '3'], capture_output=True, text=True)
    print(r.stdout.strip(), r.stderr.strip()); tried.update(qs); return len(qs)
t0 = time.time()
for rd in range(ROUNDS):
    idx, qp, qm = st.load(PREFIX); M, d = st.masks(idx)
    cand = []; nleft = 0
    for k, (m, sg, b) in M.items():
        if b > 15: continue
        q = qp if sg > 0 else qm
        left = np.where(m & (d >= b) & (q < 0))[0]; nleft += len(left)
        if not len(left): continue
        # highest dimensions first, then random
        order = left[np.lexsort((rng.random(len(left)), -d[left]))][:PER]
        for t in order:
            pat = en.parse_pattern(st.pattern(int(idx[t])))
            S = en.f_scan(en.centre(pat), nthreads=3)
            mg = sg * (S['W'][:, 2] - S['W'][:, 1])
            o = np.argsort(mg, kind='stable'); o = [int(x) for x in o if mg[x] < -1e-12 and int(x) not in tried][:TOP]
            cand += o
    print('round %d: %d targets left, %d candidates, %.0f s' % (rd, nleft, len(set(cand)), time.time() - t0), flush=True)
    if nleft == 0 or not cand: break
    run(cand)
st.report(PREFIX)
