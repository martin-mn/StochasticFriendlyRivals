"""Harvest of the stratcert kit: python3 collect.py OUTDIR [NTASKS=64] [--no-recheck]
Checks that every unit of units.tsv has exactly one well-formed line (over all t_*.tsv), that every task of the
partition wrote its .done marker, and reports per locus (W, its two families p6=p8=0 / p6=p8=1, N, E, ray+):
strata, certified, failed (with reasons), the routes used and the co-player counts by decision route; the per-unit
cost by locus and dimension.  Writes OUTDIR/stratcert_all.tsv (all lines, sorted by uid).

Audit of every certified line (independent of the task's own checks; review 2026-09-27):
  * the witness has 16 coordinates in [0,1] and lies in the OPEN stratum of its unit (0/1 exactly where the
    pattern pins them, strictly inside (0,1) where it says free);
  * the co-player counts add up to 65536, n_neg = 0, n_unres_interval = n_pos_exact + n_exact_zero and the exact
    minimum margin is >= 0;
  * the recorded exact self-play (w_CC, w_CD, w_DC, w_DD) sums to 1, its class equals the unit's class and meets
    every locus, and the sign matches the loci;
  * unless --no-recheck: the exact self-play is recomputed from the witness with cert/exact.py (gmpy2 if present,
    else fractions.Fraction) and must equal the recorded one.
The verdict is True only if nothing is missing, duplicated or malformed, every .done marker is present, the audit
finds nothing, and every retained stratum of every locus is certified."""
import sys, os, glob, collections, statistics
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
argv = [a for a in sys.argv[1:] if not a.startswith('--')]
recheck = '--no-recheck' not in sys.argv
out = argv[0]; ntasks = int(argv[1]) if len(argv) > 1 else 64
COLS = ['uid', 'sign', 'pattern', 'dim', 'class', 'loci', 'certified', 'route', 'source', 'n_tried',
        'n_pos_interval', 'n_pos_exact', 'n_struct_zero', 'n_exact_zero', 'n_neg', 'n_unres_interval',
        'min_margin', 'selfplay', 'witness', 'seconds', 'reason']          # = stratcert.COLS
U = {}
for l in open(os.path.join(HERE, 'units.tsv')):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t'); U[int(f[0])] = f
rows = {}; dup = 0; malformed = 0; H = COLS
for fn in sorted(glob.glob(os.path.join(out, 't_*.tsv'))):
    for l in open(fn):
        if l.startswith('#'):
            H = l[1:].rstrip('\n').split('\t')
            if H != COLS: print('header of %s differs from the kit columns' % fn); malformed += 1
            continue
        if not l.strip(): continue
        f = l.rstrip('\n').split('\t')
        if not l.endswith('\n') or len(f) != len(COLS):
            malformed += 1; print('malformed line in %s: %r' % (fn, l[:80])); continue
        u = int(f[0])
        if u in rows: dup += 1; print('duplicate uid %d (again in %s)' % (u, fn))
        rows[u] = dict(zip(COLS, f))
missing = sorted(set(U) - set(rows))
print('units %d, lines %d, duplicated %d, malformed %d, missing %d %s' % (len(U), len(rows), dup, malformed, len(missing), missing[:20]))
nd = [t for t in range(ntasks) if not os.path.exists(os.path.join(out, 't_%d.done' % t))]
print('tasks without .done (of %d): %d %s' % (ntasks, len(nd), nd[:20]))
for u, r in rows.items():            # the line must describe its unit
    f = U[u]; assert r['pattern'] == f[2] and r['sign'] == f[1] and r['class'] == f[6] and r['loci'] == f[7], u

# ---------------------------------------------------------------- audit of the certified lines
REQ = {'W': 'C', 'N': 'A', 'E': 'A', 'ray+': 'Z'}
SIGN = {'W': '+1', 'N': '+1', 'ray+': '+1', 'E': '-1'}
def meets(c, e): return c == e or (e == 'Z' and c in ('C', 'A'))
def class_of(W):
    """class of exact limit frequencies (w_CC, w_CD, w_DC, w_DD); a limit weight is > 0 exactly on the support"""
    if W[1] == 0 and W[2] == 0 and W[3] == 0: return 'C'
    if W[0] == 0 and W[3] == 0: return 'A'
    if W[3] == 0: return 'Z'
    return '-'
ex = None
if recheck:
    try:
        sys.path.insert(0, HERE); from cert import exact as ex
    except Exception as e:
        print('exact self-play recheck NOT possible (%s): verdict cannot be True' % e)
bad = collections.Counter(); badu = collections.defaultdict(list); nrecheck = 0
def flag(u, what): bad[what] += 1; badu[what].append(u)
for u, r in sorted(rows.items()):
    if r['certified'] != '1': continue
    try:
        w = [Fraction(x) for x in r['witness'].split(',')]
        W = [Fraction(x) for x in r['selfplay'].split(',')]
        c = [int(r[k]) for k in ('n_pos_interval', 'n_pos_exact', 'n_struct_zero', 'n_exact_zero', 'n_neg')]
        mm = Fraction(r['min_margin']); nun = int(r['n_unres_interval'])
    except Exception:
        flag(u, 'unparsable field'); continue
    if len(w) != 16 or len(W) != 4: flag(u, 'wrong number of coordinates'); continue
    if not all(0 <= x <= 1 for x in w): flag(u, 'coordinate outside [0,1]')
    if ''.join('0' if x == 0 else ('1' if x == 1 else '*') for x in w) != r['pattern']: flag(u, 'witness not in the open stratum')
    if sum(c[:4]) != 65536: flag(u, 'co-player counts do not add up to 65536')
    if c[4] != 0: flag(u, 'negative margin')
    if nun != c[1] + c[3]: flag(u, 'n_unres_interval != n_pos_exact + n_exact_zero')
    if mm < 0: flag(u, 'minimum margin < 0')
    if sum(W) != 1: flag(u, 'self-play does not sum to 1')
    if class_of(W) != r['class']: flag(u, 'self-play class != unit class')
    for loc in r['loci'].split(','):
        if not meets(r['class'], REQ[loc]): flag(u, 'class does not meet locus')
        if r['sign'] != SIGN[loc]: flag(u, 'sign does not match locus')
    if ex is not None:
        Wx = ex.lead_exact(w, w)['W']; nrecheck += 1
        if [Fraction(x) for x in Wx] != W: flag(u, 'recomputed exact self-play differs')
ncert_lines = sum(1 for r in rows.values() if r['certified'] == '1')
print('audit of %d certified lines (exact self-play recomputed for %d): %s' % (
    ncert_lines, nrecheck, 'no problem' if not bad else dict(bad)))
for k, v in badu.items(): print('   %s: uids %s' % (k, v[:20]))

def locus_rows(loc):
    if loc == 'W00': return [r for r in rows.values() if 'W' in r['loci'].split(',') and r['pattern'][6] == '0' and r['pattern'][8] == '0']
    if loc == 'W11': return [r for r in rows.values() if 'W' in r['loci'].split(',') and r['pattern'][6] == '1' and r['pattern'][8] == '1']
    return [r for r in rows.values() if loc in r['loci'].split(',')]
print('\n%-5s %6s %9s %6s  %-40s %s' % ('locus', 'strata', 'certified', 'failed', 'routes', 'co-players by route (sum over strata): interval>0 / exact>0 / structural 0 / exact 0 / neg'))
for loc in ('W', 'W00', 'W11', 'N', 'E', 'ray+'):
    R = locus_rows(loc); c = [r for r in R if r['certified'] == '1']
    routes = collections.Counter(r['route'].split(':')[0] for r in c)
    s = [sum(int(r[k]) for r in c) for k in ('n_pos_interval', 'n_pos_exact', 'n_struct_zero', 'n_exact_zero', 'n_neg')]
    print('%-5s %6d %9d %6d  %-40s %s' % (loc, len(R), len(c), len(R) - len(c), dict(routes), ' / '.join(map(str, s))))
    for r in R:
        if r['certified'] != '1': print('      FAILED uid %s %s: %s' % (r['uid'], r['pattern'], r['reason']))
EXPECT = {'W': 364, 'W00': 126, 'W11': 238, 'N': 5022, 'E': 5022, 'ray+': 7726}   # prep/mkstrata.out
def unit_locus(f, loc):
    lo = f[7].split(',')
    if loc in ('W00', 'W11'): return 'W' in lo and f[2][6] == loc[1] and f[2][8] == loc[2]
    return loc in lo
print('\ncertified strata by dimension (certified / retained, 0..8) and verdict')
badset = {u for v in badu.values() for u in v}
allproved = (not missing and dup == 0 and malformed == 0 and not nd and not bad
             and (ex is not None and nrecheck == ncert_lines))
for loc in ('W', 'W00', 'W11', 'N', 'E', 'ray+'):
    ret = collections.Counter(int(f[5]) for f in U.values() if unit_locus(f, loc))
    cer = collections.Counter(int(r['dim']) for r in locus_rows(loc) if r['certified'] == '1' and int(r['uid']) not in badset)
    nret = sum(ret.values()); ncer = sum(cer.values())
    good = (nret == EXPECT[loc] and ncer == nret); allproved &= good
    print('  %-5s %s  %d/%d %s' % (loc, ' '.join('%d/%d' % (cer[d], ret[d]) for d in range(9)), ncer, nret,
                                   'ALL CERTIFIED' if good else ('incomplete' if missing else 'NOT all certified')))
print('verdict: every retained stratum of W, N, E and ray+ has a certified exact friendly rival:', allproved)
ok = [r for r in rows.values() if r['certified'] == '1']
mm = collections.Counter(r['min_margin'] for r in ok)
print('\nexact minimum margin over the 65536 co-players at the certified witnesses:', dict(mm.most_common(5)))
print('witness = the first stored census point in %d of %d certified units' % (sum(r['route'] == 'census:0' for r in ok), len(ok)))
print('\nseconds per unit (single core) by locus and dimension: mean [max]')
for loc in ('W', 'N', 'E', 'ray+'):
    R = locus_rows(loc); byd = collections.defaultdict(list)
    for r in R: byd[int(r['dim'])].append(float(r['seconds']))
    print('  %-5s ' % loc + '  '.join('d%d %.1f [%.1f]' % (d, statistics.mean(v), max(v)) for d, v in sorted(byd.items())))
tot = sum(float(r['seconds']) for r in rows.values())
print('total %.0f core-s = %.2f core-h over %d units (%.2f s per unit)' % (tot, tot/3600, len(rows), tot/max(1, len(rows))))
with open(os.path.join(out, 'stratcert_all.tsv'), 'w') as fo:
    fo.write('#' + '\t'.join(COLS) + '\n')
    for u in sorted(rows): fo.write('\t'.join(rows[u].get(k, '') for k in COLS) + '\n')
print('wrote', os.path.join(out, 'stratcert_all.tsv'))
