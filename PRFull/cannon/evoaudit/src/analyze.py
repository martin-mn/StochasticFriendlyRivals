"""analyze.py -- confusion tables and population-share estimates of the audit (run locally, after collect.sh).

  python3 analyze.py <outdir>[,<outdir>...] [<resultsdir>]   (default out_probe,out and results)

Several output directories are merged (the probe's chunks 0-3 in out_probe, the full run's 4-127 in out).  When fewer
than all NCHUNKS chunks are complete (a chunk is complete when its chunk_<k>.done exists), the analysis is PARTIAL:
it is restricted to the sampled strategies of the complete chunks (chunk = sid mod NCHUNKS, a systematic 1-in-NCHUNKS
subsample of the sample), and every inclusion probability is multiplied by (complete chunks)/NCHUNKS.  Shares from a
partial analysis are indicative only (the subsample is systematic, not Poisson).

Reads sample/sample.npy, the design (sample.py, recomputed from the dumps in data/), and every <outdir>/chunk_*.jsonl.
Writes <resultsdir>/summary.json, confusion.txt and records.tsv.gz (one line per audited strategy).

Estimates (Poisson sampling, inclusion probabilities pi_i exactly known, sample.py):
  share of the population (samples; every game weighted equally) with property y:
      hat s = sum_{i in sample} y_i (c_i / C) / pi_i,   Var-hat = sum (1 - pi_i)/pi_i^2 (y_i c_i / C)^2,
  C = 40 000 x games of the group (Horvitz-Thompson; share_hajek = the same sum divided by the estimated total, the
  ratio form, which stays in [0, 1]).  When no audited record of a group has y, the 95% upper bound
      s <= -ln(0.05) / kappa,   kappa = min over the group's records with pi_i < 1 of -ln(1 - pi_i) C / c_i,
  holds whichever records carry y, because P(none of them drawn) = prod (1 - pi_i) <= exp(-kappa s).
"""
import os, sys, json, gzip, math, collections, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); KIT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import sample as sm

ATOMS = ['000', '001', '010', '011', '100', '110', '111']
L95 = -math.log(0.05)

NCHUNKS = 128

def load_out(outdirs):
    R = {}; dup = 0; bad = 0; done = set()
    for outdir in outdirs:
        if not os.path.isdir(outdir): continue
        for fn in sorted(os.listdir(outdir)):
            if fn.startswith('chunk_') and fn.endswith('.done'): done.add(int(fn[6:10]))
            if not (fn.startswith('chunk_') and fn.endswith('.jsonl')): continue
            for line in open(os.path.join(outdir, fn)):
                try: r = json.loads(line)
                except Exception: bad += 1; continue
                if r['sid'] in R: dup += 1
                R[r['sid']] = r
    return R, dup, bad, done

def main(outdir=os.path.join(KIT, 'out_probe') + ',' + os.path.join(KIT, 'out'), resdir=os.path.join(KIT, 'results')):
    S = np.load(os.path.join(KIT, 'sample', 'sample.npy'))
    outdirs = outdir.split(',')
    R, dup, bad, done = load_out(outdirs)
    frac = len(done)/NCHUNKS; partial = len(done) < NCHUNKS
    os.makedirs(resdir, exist_ok=True)
    G, A, h, c, Cr, hs, Nh, Ch, piG, piA, piP, pi = sm.population()
    # the sample must be the one the design draws
    rng = np.random.default_rng(sm.SEED); sel = np.where(rng.random(len(h)) < pi)[0]
    assert len(sel) == len(S) and (A['task'][sel] == S['task']).all() and (A['row'][sel] == S['row']).all()
    assert np.allclose(pi[sel], S['pi'], rtol=0, atol=0)
    if partial:                              # the complete chunks only; inclusion probabilities scaled
        keep = np.isin(np.arange(len(S)) % NCHUNKS, sorted(done)); S = S[keep]
        R = {k: v for k, v in R.items() if k % NCHUNKS in done}
        pi = pi*frac
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
        if r is not None and r.get('val') is None: r = None          # an error record: unresolved
        rows.append(dict(sid=int(s['sid']), reg=int(s['reg']), wedge='WSNE'[int(s['wedge'])], task=int(s['task']),
                         count=int(s['count']), pi=float(s['pi'])*frac, prod=format(int(s['atom']), '03b'),
                         val=None if r is None else r['val'], r=r))
    def groups():
        for rg in (1, 2):
            yield ('F%d' % rg, None), (lambda x, rg=rg: x['reg'] == rg), (reg == rg), C_of(rg)
            for w in 'WSNE':
                yield ('F%d' % rg, w), (lambda x, rg=rg, w=w: x['reg'] == rg and x['wedge'] == w), \
                      ((reg == rg) & (wid == sm.WEDGE_ID[w])), C_of(rg, w)
    def est(sel_fn, y_fn, Cg):
        s = 0.0; v = 0.0; n = 0; ny = 0; tot = 0.0
        for x in rows:
            if not sel_fn(x): continue
            n += 1; y = y_fn(x); t = x['count']/Cg/x['pi']; tot += t
            if y: ny += 1; s += t; v += (1 - x['pi'])*t*t
        return s, math.sqrt(v), n, ny, (s/tot if tot > 0 else float('nan'))
    props = {
        'changed atom': lambda x: x['r'] is not None and x['val'] != x['prod'],
        'unresolved (not audited or error)': lambda x: x['r'] is None,
        'changed E': lambda x: x['r'] is not None and 'E' in x['r']['changed'],
        'changed S': lambda x: x['r'] is not None and 'S' in x['r']['changed'],
        'changed C': lambda x: x['r'] is not None and 'C' in x['r']['changed'],
        'C 0->1 (becomes competitive)': lambda x: x['r'] is not None and x['prod'][2] == '0' and x['val'][2] == '1',
        'C 1->0 (loses competitiveness)': lambda x: x['r'] is not None and x['prod'][2] == '1' and x['val'][2] == '0',
        'S 0->1': lambda x: x['r'] is not None and x['prod'][1] == '0' and x['val'][1] == '1',
        'S 1->0': lambda x: x['r'] is not None and x['prod'][1] == '1' and x['val'][1] == '0',
        'E 0->1': lambda x: x['r'] is not None and x['prod'][0] == '0' and x['val'][0] == '1',
        'E 1->0': lambda x: x['r'] is not None and x['prod'][0] == '1' and x['val'][0] == '0',
        'production not reproduced': lambda x: x['r'] is not None and x['r']['prod_repro'] != x['prod'],
        'production C decided by the 1e-9 tolerance': lambda x: x['r'] is not None and x['r']['pr']['n_neg'] == 0 and x['r']['pr']['n_neg_tol'] > 0,
        'production S decided by the 1e-9 tolerance': lambda x: x['r'] is not None and not (x['r']['pr']['eff'] and x['r']['pr']['riv'])
                                                         and x['r']['pr']['n_gain'] == 0 and x['r']['pr']['n_gain_tol'] > 0,
        'C needed exact rationals': lambda x: x['r'] is not None and x['r']['C']['how'] == 'exact',
        'S needed exact rationals': lambda x: x['r'] is not None and x['r']['S']['how'] == 'exact',
    }
    errors = [int(k) for k, v in R.items() if v.get('val') is None]
    summary = dict(outdir=outdirs, partial=partial, chunks_complete=len(done), nchunks=NCHUNKS, pi_scale=frac,
                   sample=len(S), audited=len(R) - len(errors), errors=len(errors), error_sids=errors[:50],
                   error_msgs=sorted({R[k]['error'] for k in errors})[:20], missing=len(missing), duplicates=dup,
                   bad_lines=bad, missing_sids=missing[:50], groups={})
    for (rg, w), sf, mask, Cg in groups():
        key = rg + ('' if w is None else ' ' + w)
        k = kappa(mask, Cg); d = dict(C=Cg, kappa=k, bound95_if_none=L95/k, props={})
        # production shares of the group (exact, from the whole population of dumps)
        d['production_share'] = {a: float(c[mask & (atm == int(a, 2))].sum()/Cg) for a in ATOMS}
        d['production_share']['competitive (**1)'] = float(c[mask & ((atm & 1) == 1)].sum()/Cg)
        for name, yf in props.items():
            s, se, n, ny, hj = est(sf, yf, Cg)
            d['props'][name] = dict(n=n, n_y=ny, share=s, se=se, share_hajek=hj, **({'upper95': L95/k} if ny == 0 else {}))
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
    # wedge W: the paper's share of friendly rivals comes from production's exact count of individual-generations
    # (wfc.c ex111/settle111), which tests only strategies with p_CC,CC = 1 and p_DD,DD = 0 and uses the same double
    # tests as the dumped atoms.  Report, over the sampled W strategies, production's and the audit's friendly rivals
    # (E = C = 1) and any validated one outside that filter (which the exact count could not have seen).
    pofs = {int(s['sid']): s['p'] for s in S}
    wfr = dict(n_prod=0, n_val=0, prod_only=[], val_only=[], val_outside_ex111_filter=[])
    for x in rows:
        if x['wedge'] != 'W' or x['r'] is None: continue
        fp = x['prod'][0] == '1' and x['prod'][2] == '1'; fv = x['val'][0] == '1' and x['val'][2] == '1'
        wfr['n_prod'] += fp; wfr['n_val'] += fv
        if fp and not fv: wfr['prod_only'].append(x['sid'])
        if fv and not fp: wfr['val_only'].append(x['sid'])
        p = pofs[x['sid']]
        if fv and not (p[0] == 1.0 and p[15] == 0.0): wfr['val_outside_ex111_filter'].append(x['sid'])
    summary['W_friendly_rivals'] = wfr
    summary['how'] = {b: dict(collections.Counter(r[b]['how'] for r in au)) for b in 'ECS'}
    summary['L1_mismatch'] = sum(1 for r in au if r['E'].get('L1_mismatch'))
    # how close production's double values came to its tolerances, and how far the exact values are from them
    cv = [r for r in au if r['C']['val'] == 1]; cn = [r for r in au if r['C']['val'] == 0]
    sv = [r for r in au if r['S']['val'] == 1]; sn = [r for r in au if r['S']['val'] == 0]
    summary['tolerance_margins'] = dict(
        C_valid_most_negative_double_margin=min((r['pr']['mmin'] for r in cv), default=None),
        C_invalid_least_negative_exact_min_margin=max((r['C']['min_f'] for r in cn), default=None),
        S_valid_largest_double_gain=max((r['pr']['gmax'] for r in sv), default=None),
        S_invalid_smallest_exact_max_gain=min((r['S']['max_f'] for r in sn), default=None),
        n_C_tolerance_engaged=sum(1 for r in au if r['pr']['n_neg'] == 0 and r['pr']['n_neg_tol'] > 0),
        n_S_tolerance_engaged=sum(1 for r in au if not (r['pr']['eff'] and r['pr']['riv']) and r['pr']['n_gain'] == 0
                                  and r['pr']['n_gain_tol'] > 0))
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
        if partial:
            f.write('PARTIAL: %d of %d chunks complete; restricted to their strategies, inclusion probabilities x %.4f; '
                    'shares indicative only\n' % (len(done), NCHUNKS, frac))
        f.write('audit: %d of %d sampled strategies audited, %d errors, %d missing, %d changed atoms, %d not reproduced\n' %
                (summary['audited'], len(S), len(errors), len(missing), len(summary['changed_records']), len(summary['not_reproduced'])))
        f.write('validated 101 (efficient, competitive, not stable: would contradict the theorem): %d\n' % summary['validated_101'])
        w = summary['W_friendly_rivals']
        f.write('W friendly rivals (E = C = 1) among the sampled W strategies: production %d, validated %d, production only %d, '
                'validated only %d, validated outside the ex111 filter (p0 = 1, p15 = 0) %d\n' %
                (w['n_prod'], w['n_val'], len(w['prod_only']), len(w['val_only']), len(w['val_outside_ex111_filter'])))
        f.write('tolerances: %s\n\n' % json.dumps(summary['tolerance_margins']))
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
            p = d['props']['changed atom']; u = d['props']['unresolved (not audited or error)']
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
