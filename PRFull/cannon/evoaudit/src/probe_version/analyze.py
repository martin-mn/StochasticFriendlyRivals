"""analyze.py -- confusion tables and population-share estimates of the audit (run locally, after collect.sh).

  python3 analyze.py <outdir> [<resultsdir>]      (default outdir ../out, resultsdir ../results)

Reads sample/sample.npy, the design (sample.py, recomputed from the dumps in data/), and every <outdir>/chunk_*.jsonl.
Writes <resultsdir>/summary.json, confusion.txt and records.tsv.gz (one line per audited strategy).

Estimates (Poisson sampling, inclusion probabilities pi_i exactly known, sample.py):
  share of the population (samples; every game weighted equally) with property y:
      hat s = sum_{i in sample} y_i (c_i / C) / pi_i,   Var-hat = sum (1 - pi_i)/pi_i^2 (y_i c_i / C)^2,
  C = 40 000 x games of the group.  When no audited record of a group has y, the 95% upper bound
      s <= -ln(0.05) / kappa,   kappa = min over the group's records with pi_i < 1 of -ln(1 - pi_i) C / c_i,
  holds whichever records carry y, because P(none of them drawn) = prod (1 - pi_i) <= exp(-kappa s).
"""
import os, sys, json, gzip, math, collections, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); KIT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import sample as sm

ATOMS = ['000', '001', '010', '011', '100', '110', '111']
L95 = -math.log(0.05)

def load_out(outdir):
    R = {}; dup = 0; bad = 0
    for fn in sorted(os.listdir(outdir)):
        if not (fn.startswith('chunk_') and fn.endswith('.jsonl')): continue
        for line in open(os.path.join(outdir, fn)):
            try: r = json.loads(line)
            except Exception: bad += 1; continue
            if r['sid'] in R: dup += 1
            R[r['sid']] = r
    return R, dup, bad

def main(outdir=os.path.join(KIT, 'out'), resdir=os.path.join(KIT, 'results')):
    S = np.load(os.path.join(KIT, 'sample', 'sample.npy'))
    R, dup, bad = load_out(outdir)
    os.makedirs(resdir, exist_ok=True)
    G, A, h, c, Cr, hs, Nh, Ch, piG, piA, piP, pi = sm.population()
    # the sample must be the one the design draws
    rng = np.random.default_rng(sm.SEED); sel = np.where(rng.random(len(h)) < pi)[0]
    assert len(sel) == len(S) and (A['task'][sel] == S['task']).all() and (A['row'][sel] == S['row']).all()
    assert np.allclose(pi[sel], S['pi'], rtol=0, atol=0)
    games_w = collections.Counter(g['wedge'] for g in G)
    def C_of(reg, wedge=None, atom=None):
        if atom is not None: return Ch[sm.stratum_id(reg, sm.WEDGE_ID[wedge], int(atom, 2))]
        return sm.SPR*(512 if wedge is None else games_w[wedge])
    # kappa per group
    wid = A['wedge'].astype(int); reg = A['reg'].astype(int); atm = A['atom'].astype(int)
    lg = np.where(pi < 1, -np.log1p(-np.minimum(pi, np.nextafter(1.0, 0.0))), np.inf)
    def kappa(mask, Cg):
        m = mask & (pi < 1)
        return float((lg[m]*Cg/c[m]).min()) if m.any() else float('inf')
    missing = [int(s) for s in S['sid'] if int(s) not in R]
    rows = []
    for s in S:
        r = R.get(int(s['sid']))
        rows.append(dict(sid=int(s['sid']), reg=int(s['reg']), wedge='WSNE'[int(s['wedge'])], task=int(s['task']),
                         count=int(s['count']), pi=float(s['pi']), prod=format(int(s['atom']), '03b'),
                         val=None if r is None else r['val'], r=r))
    def groups():
        for rg in (1, 2):
            yield ('F%d' % rg, None), (lambda x, rg=rg: x['reg'] == rg), (reg == rg), C_of(rg)
            for w in 'WSNE':
                yield ('F%d' % rg, w), (lambda x, rg=rg, w=w: x['reg'] == rg and x['wedge'] == w), \
                      ((reg == rg) & (wid == sm.WEDGE_ID[w])), C_of(rg, w)
    def est(sel_fn, y_fn, Cg):
        s = 0.0; v = 0.0; n = 0; ny = 0
        for x in rows:
            if not sel_fn(x): continue
            n += 1; y = y_fn(x)
            if y: ny += 1; t = x['count']/Cg/x['pi']; s += t; v += (1 - x['pi'])*t*t
        return s, math.sqrt(v), n, ny
    props = {
        'changed atom': lambda x: x['r'] is not None and x['val'] != x['prod'],
        'unresolved (not audited)': lambda x: x['r'] is None,
        'changed E': lambda x: x['r'] is not None and 'E' in x['r']['changed'],
        'changed S': lambda x: x['r'] is not None and 'S' in x['r']['changed'],
        'changed C': lambda x: x['r'] is not None and 'C' in x['r']['changed'],
        'production not reproduced': lambda x: x['r'] is not None and x['r']['prod_repro'] != x['prod'],
        'production C decided by the 1e-9 tolerance': lambda x: x['r'] is not None and x['r']['pr']['n_neg'] == 0 and x['r']['pr']['n_neg_tol'] > 0,
        'production S decided by the 1e-9 tolerance': lambda x: x['r'] is not None and not (x['r']['pr']['eff'] and x['r']['pr']['riv'])
                                                         and x['r']['pr']['n_gain'] == 0 and x['r']['pr']['n_gain_tol'] > 0,
        'C needed exact rationals': lambda x: x['r'] is not None and x['r']['C']['how'] == 'exact',
        'S needed exact rationals': lambda x: x['r'] is not None and x['r']['S']['how'] == 'exact',
    }
    summary = dict(outdir=outdir, sample=len(S), audited=len(R), missing=len(missing), duplicates=dup, bad_lines=bad,
                   missing_sids=missing[:50], groups={})
    for (rg, w), sf, mask, Cg in groups():
        key = rg + ('' if w is None else ' ' + w)
        k = kappa(mask, Cg); d = dict(C=Cg, kappa=k, bound95_if_none=L95/k, props={})
        for name, yf in props.items():
            s, se, n, ny = est(sf, yf, Cg)
            d['props'][name] = dict(n=n, n_y=ny, share=s, se=se, **({'upper95': L95/k} if ny == 0 else {}))
        # confusion table (unweighted counts and HT population shares)
        conf = collections.Counter(); confw = collections.Counter()
        for x in rows:
            if not sf(x): continue
            v = x['val'] or 'unres'; conf[(x['prod'], v)] += 1; confw[(x['prod'], v)] += x['count']/Cg/x['pi']
        d['confusion'] = {'%s->%s' % kv: n for kv, n in sorted(conf.items())}
        d['confusion_share'] = {'%s->%s' % kv: s for kv, s in sorted(confw.items())}
        summary['groups'][key] = d
    # the theorem (efficient and competitive => stable), and how the bits were decided
    au = [x['r'] for x in rows if x['r'] is not None]
    summary['validated_101'] = sum(1 for r in au if r['val'] == '101')
    summary['how'] = {b: dict(collections.Counter(r[b]['how'] for r in au)) for b in 'ECS'}
    summary['L1_mismatch'] = sum(1 for r in au if r['E'].get('L1_mismatch'))
    cost = collections.defaultdict(list)
    for x in rows:
        if x['r'] is not None: cost[x['prod']].append((x['r']['sec'], x['r']['n_exact_calls']))
    summary['cost'] = {a: dict(n=len(v), sec_mean=float(np.mean([t for t, _ in v])), sec_max=float(np.max([t for t, _ in v])),
                               exact_mean=float(np.mean([e for _, e in v]))) for a, v in sorted(cost.items())}
    summary['sec_total'] = float(sum(r['sec'] for r in au))
    summary['changed_records'] = [dict(sid=x['sid'], reg=x['reg'], wedge=x['wedge'], task=x['task'], prod=x['prod'],
                                       val=x['val'], changed=x['r']['changed'], responsible=x['r'].get('responsible'))
                                  for x in rows if x['r'] is not None and x['val'] != x['prod']]
    summary['not_reproduced'] = [dict(sid=x['sid'], prod=x['prod'], repro=x['r']['prod_repro'], val=x['val'])
                                 for x in rows if x['r'] is not None and x['r']['prod_repro'] != x['prod']]
    json.dump(summary, open(os.path.join(resdir, 'summary.json'), 'w'), indent=1)
    # human-readable confusion tables
    with open(os.path.join(resdir, 'confusion.txt'), 'w') as f:
        f.write('audit: %d of %d sampled strategies audited, %d missing, %d changed atoms, %d not reproduced\n\n' %
                (len(R), len(S), len(missing), len(summary['changed_records']), len(summary['not_reproduced'])))
        for key, d in summary['groups'].items():
            if ' ' not in key: continue
            f.write('%s   (rows: production atom, columns: validated atom; sampled strategies)\n' % key)
            cols = ATOMS + ['101', 'unres']
            f.write('        ' + ' '.join('%6s' % a for a in cols) + '\n')
            for a in ATOMS:
                cells = [d['confusion'].get('%s->%s' % (a, b), 0) for b in cols]
                if sum(cells): f.write('  %s  ' % a + ' '.join('%6d' % n for n in cells) + '\n')
            f.write('\n')
        for key, d in summary['groups'].items():
            p = d['props']['changed atom']; u = d['props']['unresolved (not audited)']
            f.write('%-6s changed: %d of %d, share %.3g +- %.2g%s; unresolved share %.3g\n' %
                    (key, p['n_y'], p['n'], p['share'], p['se'], ('' if p['n_y'] else ' (95%% upper bound %.2g)' % p['upper95']),
                     u['share']))
    with gzip.open(os.path.join(resdir, 'records.tsv.gz'), 'wt') as f:
        f.write('sid\treg\twedge\ttask\tgid\tipt\tirep\trow\tcount\tpi\tprod\tprod_repro\tval\tchanged\thowC\thowS\t'
                'C_min\tC_qmin\tS_max\tS_qmax\tprod_mmin\tprod_gmax\tn_neg_tol\tn_gain_tol\tn_exact\tsec\n')
        for x in rows:
            r = x['r']
            if r is None:
                f.write('%d\tF%d\t%s\t%d\t\t\t\t\t%d\t%.6g\t%s\t\tUNRESOLVED\n' % (x['sid'], x['reg'], x['wedge'], x['task'], x['count'], x['pi'], x['prod'])); continue
            f.write('\t'.join(str(v) for v in [x['sid'], r['reg'], r['wedge'], r['task'], r['gid'], r['ipt'], r['irep'], r['row'],
                    r['count'], '%.6g' % r['pi'], r['prod'], r['prod_repro'], r['val'], r['changed'] or '-', r['C']['how'],
                    r['S']['how'], r['C']['min'], r['C']['q_min'], r['S']['max'], r['S']['q_max'], repr(r['pr']['mmin']),
                    repr(r['pr']['gmax']), r['pr']['n_neg_tol'], r['pr']['n_gain_tol'], r['n_exact_calls'], r['sec']]) + '\n')
    print(open(os.path.join(resdir, 'confusion.txt')).read())

if __name__ == '__main__':
    main(*sys.argv[1:])
