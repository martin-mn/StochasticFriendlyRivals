"""Collect the counts of every wsep run (out/*.summary, out/extra/*.summary) and the verification
results (out/verify_*.json, out/lemma_check.json) into out/summary.json, and print a table."""
import os, json, glob, collections
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'out')

def parse(path):
    r = dict(coplayers={}); byd = {}
    for line in open(path):
        a = line.split()
        if not a: continue
        if a[0] == 'class': r.update(cls=a[1], total=int(a[3]), processed=int(a[7]))
        elif a[0] == 'q': r['coplayers'][int(a[1])] = int(a[2])
        elif a[0] == 'seconds': r['seconds'] = float(a[1])
        elif a[0] == 'by_d':
            for t in a[2:]:
                d, v = t.split(':'); n, e, s = map(int, v.split('/')); byd[int(d)] = (n, e, s)
        elif a[0] == 'coplayers:': pass
        elif len(a) == 2:
            try: r[a[0]] = int(a[1])
            except ValueError: pass
    r['by_d'] = byd
    r['max_d_efficient'] = max((d for d, (n, e, s) in byd.items() if e), default=None)
    r['max_d_survivor'] = max((d for d, (n, e, s) in byd.items() if s), default=None)
    r['n_distinct_coplayers'] = len(r['coplayers'])
    r['top_coplayers'] = sorted(r['coplayers'].items(), key=lambda t: -t[1])[:8]
    return r

def main():
    S = {}
    for p in sorted(glob.glob(os.path.join(OUT, '*.summary')) + glob.glob(os.path.join(OUT, 'extra', '*.summary'))):
        name = os.path.relpath(p, OUT)[:-len('.summary')]
        S[name] = parse(p)
    V = {}
    for p in glob.glob(os.path.join(OUT, 'verify_*.json')): V.update(json.load(open(p)))
    allq = collections.Counter()
    for name, r in S.items():
        if name in ('c6', 'c8', 'c68'): allq.update(r['coplayers'])
        if name in V: r['verify'] = V[name]
    L = json.load(open(os.path.join(OUT, 'lemma_check.json'))) if os.path.exists(os.path.join(OUT, 'lemma_check.json')) else None
    extra = {}
    for f in ['exact_table.json', 'exact_corners.json', os.path.join('cover', 'sep_cover.json'), 'lemma_check.json']:
        fp = os.path.join(OUT, f)
        if os.path.exists(fp): extra[f] = json.load(open(fp))
    summary = dict(runs=S, checks=extra, separation_coplayers_distinct=len(allq),
                   separation_top_coplayers=allq.most_common(10), lemma_check=L)
    json.dump(summary, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1, default=str)
    hdr = '%-12s %10s %10s %10s %8s %8s %6s %6s %5s %5s %6s %8s %s'
    print(hdr % ('class', 'strata', 'efficient', 'by ALLD', 'other q', 'by pool', 'full', 'NONE', 'dEff', 'dSrv', '#q', 'seconds', 'verified'))
    for name, r in S.items():
        v = r.get('verify')
        vs = '-' if v is None else ('rows ok %d/%d, resample mism %d, interior>=0 %d' %
                                    (v['rows'] - v['rows_failing'], v['rows'], v['exact_resample_mismatches'], v['interior_nonnegative']))
        print(hdr % (name, r['total'], r.get('efficient'), r.get('excluded_by_ALLD'), r.get('survivors'),
                     r.get('excluded_pool'), r.get('excluded_full'), r.get('NOT_EXCLUDED'), r['max_d_efficient'],
                     r['max_d_survivor'], r['n_distinct_coplayers'], r.get('seconds'), vs))
    print('separation (c6+c8+c68): %d distinct co-players besides ALLD; top: %s' % (len(allq), allq.most_common(6)))

if __name__ == '__main__':
    main()
