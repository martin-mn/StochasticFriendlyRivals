"""Direct exact checks of the structural lemma on the strata of the W-separation certificate, and the
two positive controls that make the separated set disconnected (both parts non-empty).

(1) L3 + multilinearity, exactly.  For sampled (stratum, q) -- the certifying q and a random q --
    the corner values N_q(c) are formed from the corners' OWN binary limits by the rule of L3
    (D*(c) x margin(c) if E_c = E, else 0, and E_c < E must not occur), and the multilinear
    interpolant of these values is compared, in exact rationals, with N_q(x) = D*(x) x margin(x)
    computed directly at random rational interior points x.
(2) finite eps: exact stationary margins of the clamped chain at eps = 1e-3, 1e-4, 1e-5
    (matrix-tree determinants, no leading-order arithmetic) converge to the limit value, and are
    negative, at interior points of certified strata.
(3) positive controls: the binary friendly rival 19079 (p6 = p8 = 0) and the exact point eq:certW
    (p6 = p8 = 1) are certified friendly rivals of W by cert.certify.certify_point.
usage: nice -n 10 python3 lemma_check.py
"""
import sys, os, random, itertools, time, json
from fractions import Fraction
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from cert import engine as en, exact as ex, certify as ce

def N_at(x, q):
    L = ex.lead_exact(x, ex.g2p(q)); return L['E'], L['D']*(L['W'][2] - L['W'][1])

def check_ml(args):
    pat, q, seed = args
    rng = random.Random(seed)
    fr = [j for j, ch in enumerate(pat) if ch == '*']; d = len(fr)
    base = [Fraction(1, 2) if ch == '*' else Fraction(int(ch)) for ch in pat]
    E, _ = N_at(base, q)
    Nc = {}; nhigh = 0
    for bits in itertools.product((0, 1), repeat=d):
        c = list(base)
        for j, b in zip(fr, bits): c[j] = Fraction(b)
        Ec, Nv = N_at(c, q)
        if Ec < E: return (pat, q, 'E_c<E')
        if Ec > E: Nv = Fraction(0); nhigh += 1
        Nc[bits] = Nv
    for t in range(3):
        x = list(base)
        xs = [Fraction(rng.randint(1, 99), 100) for _ in fr]
        for j, v in zip(fr, xs): x[j] = v
        Ex, Nx = N_at(x, q)
        if Ex != E: return (pat, q, 'E(x) != E')
        interp = Fraction(0)
        for bits, Nv in Nc.items():
            if Nv == 0: continue
            w = Fraction(1)
            for v, b in zip(xs, bits): w *= v if b else 1 - v
            interp += Nv*w
        if interp != Nx: return (pat, q, 'interpolation mismatch')
    return (pat, q, 'ok', nhigh, 1 << d)

def fe_check(args):
    pat, q, seed = args
    rng = random.Random(seed)
    x = [Fraction(rng.randint(1, 9), 10) if ch == '*' else Fraction(int(ch)) for ch in pat]
    lim = ex.margin(x, ex.g2p(q), +1)
    out = []
    for k in (3, 4, 5):
        pi, A = ex.finite_eps(x, ex.g2p(q), Fraction(1, 10**k))
        m = sum(pi[4*2:4*3]) - sum(pi[4*1:4*2])
        out.append(float(m))
    return pat, q, float(lim), out

def main():
    t0 = time.time(); rng = random.Random(2026)
    rows = []
    for name in ['c6', 'c8', 'c68', 'c9a', 'c9b']:
        R = [l.split('\t') for l in open(os.path.join(HERE, 'out', name + '.tsv'))]
        R = [r for r in R if int(r[2]) <= 6]
        rows += [(name, r[1], int(r[4])) for r in rng.sample(R, 40)]
    jobs = []
    for k, (name, pat, q) in enumerate(rows):
        jobs.append((pat, q, k)); jobs.append((pat, rng.randrange(65536), 10**6 + k))
    with Pool(3) as P:
        res = P.map(check_ml, jobs, chunksize=4)
    ok = sum(1 for r in res if r[2] == 'ok'); bad = [r for r in res if r[2] != 'ok']
    nh = sum(r[3] for r in res if r[2] == 'ok'); nc = sum(r[4] for r in res if r[2] == 'ok')
    rep = dict(multilinear_pairs=len(jobs), multilinear_ok=ok, multilinear_bad=bad[:5], corners=nc, corners_E_c_gt_E=nh,
               interior_points=3*len(jobs))
    print('multilinear/L3:', json.dumps(rep), round(time.time() - t0, 1), 's', flush=True)
    fej = [(pat, q, k) for k, (name, pat, q) in enumerate(rows[::20])]
    with Pool(3) as P:
        fe = P.map(fe_check, fej)
    fe_ok = all(o[-1] < 0 and abs(o[-1] - lim) < abs(o[0] - lim) + 1e-15 and abs(o[-1] - lim) < 1e-3 for _, _, lim, o in fe)
    rep['finite_eps'] = [dict(pattern=p, q=q, limit=lim, eps_1e3_1e4_1e5=o) for p, q, lim, o in fe]
    rep['finite_eps_ok'] = fe_ok
    print('finite eps ok:', fe_ok, flush=True)
    for p, q, lim, o in fe: print('  ', p, q, 'limit %.6g' % lim, ' eps=1e-3,1e-4,1e-5:', ['%.6g' % v for v in o])
    # positive controls
    P19079 = [Fraction((19079 >> j) & 1) for j in range(16)]
    PcW = [Fraction(v) for v in ['1', '1/2', '1', '1/2', '0', '0', '1', '1/2', '1', '1/4', '1', '1/4', '1/4', '0', '1/4', '0']]
    for nm, p in [('binary 19079', P19079), ('eq:certW', PcW)]:
        c = ce.certify_point(p, wedge='W')
        print('control', nm, 'p6,p8,p9 =', p[6], p[8], p[9], 'certified FR of W:', c['certified'],
              'class', c['efficiency_class'], 'neg', c['n_neg'], 'min margin', c['min_margin'], flush=True)
        rep['control_' + nm] = dict(p6=str(p[6]), p8=str(p[8]), p9=str(p[9]), certified=bool(c['certified']),
                                    efficiency_class=c['efficiency_class'], n_neg=int(c['n_neg']),
                                    n_pos=int(c['n_pos']), n_zero=int(c['n_zero']), min_margin=str(c['min_margin']))
    rep['seconds'] = round(time.time() - t0, 1)
    json.dump(rep, open(os.path.join(HERE, 'out', 'lemma_check.json'), 'w'), indent=1)

if __name__ == '__main__':
    main()
