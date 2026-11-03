"""Independent re-verification of the W-separation certificates written by wsep (out/<cls>.tsv).

For EVERY certified (stratum, q) row:
  (a) exact centre margin: ex.margin at the centre (all free coordinates 1/2), in rationals
      (gmpy2.mpq), must be < 0, lie inside the interval the C driver recorded, and the exact E
      must equal the recorded E;
  (b) corner table re-computed by the library's en.corner_scan (a separate code path: no early
      exit, x86_64 slice of libcert under Rosetta when run from the Anaconda python): must give
      neg_ok, no positive / unresolved / E_c<E corner, and the same counts nneg, nzs, nzi, nhigh.
For a random sample of strata of each class (any status), the efficiency and ALLD decisions are
recomputed exactly (ex.lead_exact at the centre, self-play and against ALLD) and compared with the
exponent-only decisions of the driver.  For a random sample of certificates, exact margins at
random rational interior points (including points 1/1000 from the faces) must all be < 0 -- a
direct check of the multilinear lemma (L3 + L4).

usage: nice -n 10 python3 verify.py [nproc<=3] [classes...]
"""
import sys, os, time, json, random
from fractions import Fraction
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from cert import engine as en, exact as ex

def centre_point(pat):
    return [Fraction(1, 2) if ch == '*' else Fraction(int(ch)) for ch in pat]

def check_row(row):
    cls, pat, d, E, q, how, nneg, nzs, nzi, nhigh, clo, chi, ntr = row
    d = int(d); E = int(E); q = int(q); nneg, nzs, nzi, nhigh = map(int, (nneg, nzs, nzi, nhigh))
    clo = float(clo); chi = float(chi)
    bad = []
    if q < 0: return (pat, None, ['NOT EXCLUDED'])
    c = centre_point(pat)
    L = ex.lead_exact(c, ex.g2p(q))
    m = L['W'][2] - L['W'][1]
    if not m < 0: bad.append('exact centre margin not negative')
    if not (Fraction(clo) <= m <= Fraction(chi)): bad.append('exact centre margin outside the C enclosure')
    if L['E'] != E: bad.append('E mismatch')
    C = en.corner_scan(pat, +1, [q], nthreads=1)
    cnt = C['cnt'][0]
    if not C['neg_ok'][0]: bad.append('library corner scan: not neg_ok')
    if int(C['E'][0]) != E: bad.append('library E mismatch')
    got = (int(cnt[en.C_NEG]), int(cnt[en.C_ZS]), int(cnt[en.C_ZI]), int(cnt[en.C_HIGH]))
    if got != (nneg, nzs, nzi, nhigh): bad.append('corner counts differ %s vs %s' % (got, (nneg, nzs, nzi, nhigh)))
    if cnt[en.C_POS] or cnt[en.C_LOW] or cnt[en.C_UNR]: bad.append('pos/low/unresolved corner')
    if nneg < 1: bad.append('no negative corner')
    return (pat, str(m), bad)

# ---- enumeration identical to wsep.c (for the sampled exact re-decisions)
CLASSES = {
    'c6':  ([0, 6], [(1, 2)]),
    'c8':  ([0, 6, 8], [(1, 0, 2), (1, 1, 2)]),
    'c68': ([0, 6, 8], [(1, 0, 1), (1, 1, 0)]),
    'c9a': ([0, 6, 8, 9], [(1, 0, 0, 0), (1, 0, 0, 2)]),
    'c9b': ([0, 6, 8, 9], [(1, 1, 1, 0), (1, 1, 1, 1)]),
}
TEMPLATES = {}
def parse_template(t):
    """'1...[1*]...........' -> per-coordinate lists of allowed values (0, 1, 2 = free)"""
    out = []; k = 0
    while k < len(t):
        ch = t[k]
        if ch == '[':
            e = t.index(']', k); out.append([2 if c == '*' else int(c) for c in t[k+1:e]]); k = e + 1
        elif ch == '.': out.append([0, 1, 2]); k += 1
        else: out.append([2 if ch == '*' else int(ch)]); k += 1
    assert len(out) == 16, t
    return out
def decode(name, i):
    if name in TEMPLATES:
        pat = []
        for vals in TEMPLATES[name]:
            pat.append(vals[i % len(vals)]); i //= len(vals)
        return ''.join('01*'[v] for v in pat)
    spec, assign = CLASSES[name]; rest = [j for j in range(16) if j not in spec]
    a = i % len(assign); r = i // len(assign); pat = [0]*16
    for t, j in enumerate(spec): pat[j] = assign[a][t]
    for j in rest: pat[j] = r % 3; r //= 3
    return ''.join('01*'[v] for v in pat)
def total(name):
    if name in TEMPLATES:
        n = 1
        for v in TEMPLATES[name]: n *= len(v)
        return n
    spec, assign = CLASSES[name]; return len(assign)*3**(16 - len(spec))

def exact_status(pat):
    """exact: ('ineff' | 'alld' | 'surv') from ex.lead_exact at the centre"""
    c = centre_point(pat)
    Wself = ex.lead_exact(c, c)['W']
    if Wself[0] != 1: return 'ineff'
    Wa = ex.lead_exact(c, ex.g2p(0))['W']
    assert Wa[2] == 0                        # ALLD: w_DC = 0 exactly in the limit
    return 'alld' if Wa[1] > 0 else 'surv'

def interior_check(args):
    pat, q, seed = args
    rng = random.Random(seed); bad = 0; worst = None
    fr = [j for j, ch in enumerate(pat) if ch == '*']
    for t in range(4):
        x = centre_point(pat)
        for j in fr:
            if t == 0: x[j] = Fraction(rng.choice([1, 999]), 1000)          # near a corner
            elif t == 1: x[j] = Fraction(rng.randint(1, 999), 1000)
            elif t == 2: x[j] = Fraction(rng.randint(1, 96), 97)
            else: x[j] = Fraction(1, 1000) if rng.random() < 0.5 else Fraction(rng.randint(1, 9), 10)
        m = ex.margin(x, ex.g2p(q), +1)
        if not m < 0: bad += 1
        worst = m if worst is None or m > worst else worst
    return pat, q, bad, float(worst)

def main():
    nproc = min(int(sys.argv[1]) if len(sys.argv) > 1 else 3, 3)
    names = sys.argv[2:] or ['c6', 'c8', 'c68', 'c9a', 'c9b']
    rep = {}
    for name in names:
        t0 = time.time()
        base = os.path.join(HERE, 'out', name)
        head = open(base + '.summary').readline().split()
        if head[1].startswith('t:'): TEMPLATES[name] = parse_template(head[1][2:])
        rows = [l.rstrip('\n').split('\t') for l in open(base + '.tsv')]
        with Pool(nproc) as P:
            res = P.map(check_row, rows, chunksize=500)
        nbad = sum(1 for r in res if r[2]); worst = max((Fraction(r[1]) for r in res if r[1] is not None), default=None)
        with open(base + '.exact.tsv', 'w') as f:
            for row, r in zip(rows, res):
                f.write('\t'.join(row[:5] + [r[1] or 'NA', ';'.join(r[2]) or 'ok']) + '\n')
        t1 = time.time()
        # sampled exact re-decision of efficiency / ALLD for strata of every status
        rng = random.Random(sum(map(ord, name)))
        surv = set(r[1] for r in rows)
        idx = [rng.randrange(total(name)) for _ in range(4000)]
        pats = [decode(name, i) for i in idx]
        with Pool(nproc) as P:
            st = P.map(exact_status, pats, chunksize=100)
        mism = 0
        for p, s in zip(pats, st):
            dr = 'surv' if p in surv else None
            if s == 'surv' and p not in surv: mism += 1
            if s != 'surv' and p in surv: mism += 1
            if s == 'ineff' and en.efficiency_class(p) == 'C': mism += 1
            if s != 'ineff' and en.efficiency_class(p) != 'C': mism += 1
        from collections import Counter
        stc = Counter(st)
        # interior points of a sample of certificates
        samp = rng.sample(rows, min(1500, len(rows)))
        with Pool(nproc) as P:
            ic = P.map(interior_check, [(r[1], int(r[4]), k) for k, r in enumerate(samp)], chunksize=50)
        nib = sum(r[2] for r in ic); iworst = max((r[3] for r in ic), default=None)
        rep[name] = dict(rows=len(rows), rows_failing=nbad, max_exact_centre_margin=str(worst),
                         exact_resample=len(pats), exact_status_counts=dict(stc), exact_resample_mismatches=mism,
                         interior_points=4*len(samp), interior_nonnegative=nib, interior_max_margin=iworst,
                         seconds_rows=round(t1 - t0, 1), seconds_total=round(time.time() - t0, 1))
        print(name, json.dumps(rep[name]), flush=True)
        if nbad:
            for r in res:
                if r[2]: print('  FAIL', r)
    json.dump(rep, open(os.path.join(HERE, 'out', 'verify_%s.json' % '_'.join(n.replace('/', '-') for n in names)), 'w'), indent=1)

if __name__ == '__main__':
    main()
