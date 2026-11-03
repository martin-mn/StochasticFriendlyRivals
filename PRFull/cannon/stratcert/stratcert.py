"""stratcert: certify ONE exact friendly rival in every retained friendly-rival support stratum of the
memory-two cube, for wedges W, N, E and the ray T>S (referee 2, major 1 / P1).

A UNIT is (stratum, sign of T-S) from units.tsv (built by prep/mkstrata.py).  For a unit, the task tries
candidate points p in this order and stops at the first one that is certified:
  1. census  : the points the M2f census stored for this stratum and sign (points.tsv; exact doubles);
  2. mirror  : (class A only) the census points of the C/D-mirror unit, p_s -> 1 - p_{15-s} (exact rationals);
  3. round   : each census point rounded to the dyadic grid 2^-10, then 2^-16 (free coordinates kept in (0,1));
  4. random  : up to NRANDOM dyadic points n/2^20 in [0.02, 0.98] (seeded by the uid), screened in floats
               against the 1681 co-players of pool_q.txt and then against all 65536 before being certified.
A point is CERTIFIED when
  * its support pattern equals the unit's stratum (pinned coordinates exactly 0/1, free ones in (0,1));
  * the exact efficiency class agrees three ways: cls.bin (ub.c `classes`), en.efficiency_class(pattern)
    (lemma L1) and the exact self-play support at p (exact.py) -- and it meets every locus of the unit;
  * ce.certify_point(p, sign) finds no negative margin sign*(w_DC - w_CD) against any of the 65536 binary
    co-players (interval enclosure strictly signed, structural zero, or exact rational), so p is competitive
    against every co-player of any memory (SI Proposition).
The unit is then certified for all its loci: sign +, class C -> W and ray+; class A -> N and ray+;
class Z -> ray+; sign -, class A -> E.

usage: python3 stratcert.py TASK NTASKS OUTDIR [--nproc N] [--only UID,UID,...] [--nrandom K]
  units with uid % NTASKS == TASK; output OUTDIR/t_TASK.tsv, one line per unit, appended as units finish;
  on restart the uids already in the file are skipped (idempotent); OUTDIR/t_TASK.done at the end.
"""
import os, sys, time, argparse, random
os.environ.setdefault('CERT_MAX_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from fractions import Fraction
from cert import engine as en, exact as ex, certify as ce

COLS = ['uid', 'sign', 'pattern', 'dim', 'class', 'loci', 'certified', 'route', 'source', 'n_tried',
        'n_pos_interval', 'n_pos_exact', 'n_struct_zero', 'n_exact_zero', 'n_neg', 'n_unres_interval',
        'min_margin', 'selfplay', 'witness', 'seconds', 'reason']
CLSCODE = {0: '-', 1: 'C', 2: 'A', 3: 'Z'}
REQ = {'W': 'C', 'N': 'A', 'E': 'A', 'ray+': 'Z'}

def load_units():
    U = {}
    for l in open(os.path.join(HERE, 'units.tsv')):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        U[int(f[0])] = dict(uid=int(f[0]), sign=int(f[1]), pattern=f[2], ub_index=int(f[3]), census_code=int(f[4]),
                            dim=int(f[5]), cls=f[6], loci=f[7].split(','), npts=int(f[8]))
    return U

def load_points():
    P = {}
    for l in open(os.path.join(HERE, 'points.tsv')):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        P.setdefault(int(f[0]), []).append((f[2], [Fraction(float(x)) for x in f[3:19]]))
    return P

def mirror_pattern(pat): return ''.join({'0': '1', '1': '0', '*': '*'}[pat[15 - s]] for s in range(16))
def mirror_point(p): return [1 - p[15 - s] for s in range(16)]
def frac_str(x): return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)

_cls_bin = None
def cls_bin(i):
    global _cls_bin
    if _cls_bin is None:
        fn = os.path.join(HERE, 'cls.bin')
        _cls_bin = np.memmap(fn, dtype=np.uint8, mode='r') if os.path.exists(fn) else False
    return None if _cls_bin is False else CLSCODE[int(_cls_bin[i])]

_pool = None
def pool_q():
    global _pool
    if _pool is None:
        _pool = np.array([int(x) for x in open(os.path.join(HERE, 'pool_q.txt')).read().split()], dtype=np.int64)
    return _pool

def round_point(p, pat, bits):
    g = 1 << bits; out = []
    for x, c in zip(p, pat):
        if c != '*': out.append(x); continue
        n = round(x*g); n = min(max(n, 1), g - 1); out.append(Fraction(n, g))
    return out

def random_points(u, k):
    rng = random.Random(1000003*u['uid'] + 17*(u['sign'] + 1) + 5)
    g = 1 << 20; lo = -(-int(0.02*g*1000) // 1000); hi = int(0.98*g)
    for _ in range(k):
        yield [Fraction(int(c)) if c != '*' else Fraction(rng.randint(lo, hi), g) for c in u['pattern']]

def float_screen(p, sign):
    """cheap float pre-screen (proposes only; never used as proof): min margin over the pool, then all q"""
    x = np.array([float(v) for v in p])
    S = en.f_scan(x, pool_q(), nthreads=1)
    if (sign*(S['W'][:, 2] - S['W'][:, 1])).min() < -1e-9: return False
    S = en.f_scan(x, None, nthreads=1)
    return (sign*(S['W'][:, 2] - S['W'][:, 1])).min() >= -1e-9

def try_point(u, p, sign):
    """returns (ok, reason, result dict of certify_point or None)"""
    pat = u['pattern']
    if en.pattern_str(en.pattern_of(p)) != pat: return False, 'point not in stratum', None
    r = ce.certify_point(p, sign=sign, eff=None, nthreads=1, nproc=1)
    ecls = r['efficiency_class']
    if r['n_neg'] > 0: return False, 'negative margin (%d co-players; min %s at q=%s)' % (r['n_neg'], frac_str(r['min_margin']) if isinstance(r['min_margin'], Fraction) else r['min_margin'], r['argmin']), r
    if ecls != u['cls']: return False, 'class at point %s != stratum %s' % (ecls, u['cls']), r
    if not all(en.meets(ecls, REQ[l]) for l in u['loci']): return False, 'class %s does not meet loci' % ecls, r
    return True, '', r

def certify_unit(args):
    u, points, mirror_pts, nrandom = args
    t0 = time.time(); sign = u['sign']; pat = u['pattern']
    # exact class of the stratum, two routes (plus the point's own, inside try_point)
    c1 = en.efficiency_class(pat); c2 = cls_bin(u['ub_index'])
    base = dict(uid=u['uid'], sign='%+d' % sign, pattern=pat, dim=u['dim'], loci=','.join(u['loci']))
    if c1 != u['cls'] or (c2 is not None and c2 != u['cls']):
        return dict(base, **{'class': u['cls'], 'certified': 0, 'route': 'none', 'source': '', 'n_tried': 0,
                    'seconds': '%.2f' % (time.time() - t0), 'reason': 'class mismatch: units %s, L1 %s, cls.bin %s' % (u['cls'], c1, c2)})
    cands = []
    for k, (src, p) in enumerate(points): cands.append(('census:%d' % k, src, p))
    for k, (src, p) in enumerate(mirror_pts): cands.append(('mirror:%d' % k, 'mirror of ' + src, mirror_point(p)))
    for bits in (10, 16):
        for k, (src, p) in enumerate(points): cands.append(('round%d:%d' % (bits, k), src + ' rounded to 2^-%d' % bits, round_point(p, pat, bits)))
    ntried = 0; last = 'no candidate'; lastr = None; seen = set()
    def attempt(route, src, p):
        nonlocal ntried, last, lastr
        key = tuple(p)
        if key in seen: return None
        seen.add(key); ntried += 1
        ok, why, r = try_point(u, p, sign)
        if ok: return (route, src, p, r)
        last = '%s: %s' % (route, why); lastr = r
        return None
    hit = None
    for route, src, p in cands:
        hit = attempt(route, src, p)
        if hit: break
    if not hit and nrandom:
        nscr = 0
        for k, p in enumerate(random_points(u, nrandom)):
            if not float_screen(p, sign): continue
            nscr += 1
            hit = attempt('random:%d' % k, 'dyadic draw %d (seed uid)' % k, p)
            if hit or nscr >= 8: break
        if not hit: last += '; random: %d of %d draws passed the float screen' % (nscr, nrandom)
    out = dict(base, **{'class': u['cls'], 'n_tried': ntried})
    if hit:
        route, src, p, r = hit
        n_unres = r['n_unresolved_by_interval']; n_pos_exact = n_unres - r['n_exact_zero']
        mm = r['min_margin']
        out.update(certified=1, route=route, source=src, n_pos_interval=r['n_pos'] - n_pos_exact, n_pos_exact=n_pos_exact,
                   n_struct_zero=r['n_struct_zero'], n_exact_zero=r['n_exact_zero'], n_neg=r['n_neg'], n_unres_interval=n_unres,
                   min_margin=frac_str(Fraction(mm)), selfplay=','.join(frac_str(Fraction(w)) for w in r['selfplay']),
                   witness=','.join(frac_str(x) for x in p), reason='')
        assert out['n_pos_interval'] + n_pos_exact + r['n_struct_zero'] + r['n_exact_zero'] == en.NQ
    else:
        out.update(certified=0, route='none', source='', reason=last)
        if lastr is not None: out.update(n_neg=lastr['n_neg'], n_unres_interval=lastr['n_unresolved_by_interval'])
    out['seconds'] = '%.2f' % (time.time() - t0)
    return out

def fmt(d): return '\t'.join(str(d.get(c, '')) for c in COLS)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('task', type=int); ap.add_argument('ntasks', type=int); ap.add_argument('outdir')
    ap.add_argument('--nproc', type=int, default=1); ap.add_argument('--only', default='')
    ap.add_argument('--nrandom', type=int, default=4096)
    a = ap.parse_args()
    U = load_units(); P = load_points()
    byp = {(u['pattern'], u['sign']): u['uid'] for u in U.values()}
    uids = [i for i in sorted(U) if i % a.ntasks == a.task]
    if a.only: uids = [int(x) for x in a.only.split(',')]
    os.makedirs(a.outdir, exist_ok=True)
    fn = os.path.join(a.outdir, 't_%d.tsv' % a.task)
    done = set()
    if os.path.exists(fn):
        keep = []; drop = 0
        for l in open(fn):
            if l.startswith('#'): keep.append(l); continue
            f = l.rstrip('\n').split('\t')
            if l.endswith('\n') and len(f) == len(COLS) and int(f[0]) not in done: done.add(int(f[0])); keep.append(l)
            else: drop += 1                      # a line cut by a killed task (or a duplicate): redo that unit
        if drop:
            with open(fn, 'w') as fo: fo.writelines(keep)
            print('dropped %d incomplete line(s) of %s' % (drop, fn))
    todo = [i for i in uids if i not in done]
    jobs = []
    for i in todo:
        u = U[i]; mp = []
        if u['cls'] == 'A':
            j = byp.get((mirror_pattern(u['pattern']), -u['sign']))
            if j is not None: mp = P.get(j, [])
        jobs.append((u, P.get(i, []), mp, a.nrandom))
    t0 = time.time()
    newfile = not os.path.exists(fn)
    with open(fn, 'a') as fo:
        if newfile: fo.write('#' + '\t'.join(COLS) + '\n'); fo.flush()
        if a.nproc > 1:
            from multiprocessing import Pool
            with Pool(a.nproc) as pool:
                for r in pool.imap_unordered(certify_unit, jobs, chunksize=1):
                    fo.write(fmt(r) + '\n'); fo.flush()
        else:
            for j in jobs:
                fo.write(fmt(certify_unit(j)) + '\n'); fo.flush()
    # completeness check: every assigned uid has exactly one line
    have = {}
    for l in open(fn):
        if l.startswith('#') or not l.strip(): continue
        f = l.split('\t'); have[int(f[0])] = have.get(int(f[0]), 0) + 1
    missing = [i for i in uids if i not in have]; dup = [i for i, c in have.items() if c > 1]
    ncert = 0
    for l in open(fn):
        if l.startswith('#') or not l.strip(): continue
        ncert += int(l.split('\t')[6])
    msg = 'task %d/%d: %d units assigned, %d skipped (done before), %d run in %.1f s; certified %d; missing %d; duplicated %d' % (
        a.task, a.ntasks, len(uids), len(done), len(todo), time.time() - t0, ncert, len(missing), len(dup))
    print(msg)
    if not missing and not dup and not a.only:
        with open(os.path.join(a.outdir, 't_%d.done' % a.task), 'w') as f: f.write(msg + '\n')
    return 0 if not missing and not dup else 3

if __name__ == '__main__':
    sys.exit(main())
