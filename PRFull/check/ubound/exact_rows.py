"""EXACT rows for the certificate co-players: for every binary focal c (65536) against q, the exact
minimum tree exponent E and the exact sign of w_DC - w_CD (gmpy2 rationals, check/cert/exact.py),
written in ub's row format (E << 2 | 2 if > 0, 1 if < 0, 0 if = 0).  Then `ub rowcmp` checks
that the interval rows used by the DP are byte-identical."""
import sys, os, time, subprocess, numpy as np
from multiprocessing import Pool
sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/check')
from cert import exact as ex
def chunk(args):
    q, a, b = args; Q = ex.g2p(q); out = np.zeros(b - a, np.uint8)
    for c in range(a, b):
        p = ex.g2p(c); L = ex.lead_exact(p, Q, as_fraction=False); m = L['W'][2] - L['W'][1]
        f = 2 if m > 0 else (1 if m < 0 else 0); assert 0 <= L['E'] < 64
        out[c - a] = (L['E'] << 2) | f
    return q, a, out
if __name__ == '__main__':
    qs = [int(x) for x in open(sys.argv[1]).read().split()]; os.makedirs('runs/xrows', exist_ok=True)
    jobs = [(q, a, a + 2048) for q in qs for a in range(0, 65536, 2048)]; rows = {q: np.zeros(65536, np.uint8) for q in qs}
    t0 = time.time()
    with Pool(3) as P:
        for q, a, out in P.imap_unordered(chunk, jobs): rows[q][a:a + len(out)] = out
    print('exact rows for %d co-players in %.0f s' % (len(qs), time.time() - t0), flush=True)
    bad = 0
    for q in qs:
        fn = 'runs/xrows/%d.bin' % q; rows[q].tofile(fn)
        r = subprocess.run(['./ub', 'rowcmp', str(q), fn], capture_output=True, text=True); print(r.stdout.strip()); bad += r.returncode != 0
    print('rows differing from the interval engine:', bad)
