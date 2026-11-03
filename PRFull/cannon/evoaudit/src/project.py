"""project.py -- per-strategy cost from the probe and the projected cost of the full run.

  python3 project.py [probe_outdir] [nchunks]          (default ../out_probe, 128)

Cost model: the audit time of a strategy ('sec', measured inside audit(), one core on a Cannon sapphire node) is drawn
from the probe's empirical distribution of its production atom (the cost is driven by how many co-players go to exact
rationals, which depends on the atom far more than on the regime or wedge).  Per chunk: the sum over its strategies plus
the probe's measured per-chunk overhead (chunk wall time from chunk_<k>.done minus the sum of 'sec').  1000 bootstrap
replicates give the distribution of every chunk's time; reported: the mean total, the mean and the 99.9th percentile of
the slowest chunk.
"""
import os, sys, json, glob, collections, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); KIT = os.path.dirname(HERE)

def main(probe=os.path.join(KIT, 'out_probe'), nchunks=128, nboot=1000, seed=7):
    nchunks = int(nchunks)
    S = np.load(os.path.join(KIT, 'sample', 'sample.npy'))
    R = []; over = []
    for fn in sorted(glob.glob(os.path.join(probe, 'chunk_*.jsonl'))):
        k = int(os.path.basename(fn)[6:10]); rr = [json.loads(l) for l in open(fn)]; R += rr
        dn = fn[:-6] + '.done'
        if os.path.exists(dn):
            wall = float(open(dn).read().split(' in ')[1].split()[0]); over.append(wall - sum(r['sec'] for r in rr))
    by = collections.defaultdict(list)
    for r in R: by[r['prod']].append(r['sec'])
    by = {a: np.array(v) for a, v in by.items()}
    out = dict(probe_records=len(R), per_atom={}, overhead_per_chunk=dict(mean=float(np.mean(over)), max=float(np.max(over))))
    for a in sorted(by):
        v = by[a]; out['per_atom'][a] = dict(n=int(len(v)), mean=float(v.mean()), median=float(np.median(v)),
                                             p99=float(np.quantile(v, 0.99)), max=float(v.max()))
    atoms = np.array([format(int(x), '03b') for x in S['atom']]); chunk = np.arange(len(S)) % nchunks
    n_by = collections.Counter(atoms.tolist())
    out['sample_by_atom'] = dict(n_by)
    out['expected_total_sec'] = float(sum(n_by[a]*by[a].mean() for a in n_by)) + nchunks*float(np.mean(over))
    rng = np.random.default_rng(seed)
    T = np.zeros((nboot, nchunks))
    for a in n_by:
        m = atoms == a; ck = chunk[m]
        draws = rng.choice(by[a], size=(nboot, int(m.sum())), replace=True)
        for b in range(nboot): T[b] += np.bincount(ck, weights=draws[b], minlength=nchunks)
    T += np.mean(over)
    mx = T.max(1)
    out['chunk_sec'] = dict(mean=float(T.mean()), sd=float(T.std()), max_mean=float(mx.mean()),
                            max_p999=float(np.quantile(mx, 0.999)))
    out['core_hours'] = out['expected_total_sec']/3600
    print(json.dumps(out, indent=1))
    return out

if __name__ == '__main__':
    main(*sys.argv[1:])
