"""sample.py -- the stratified Poisson sample of stored strategies for the audit (writes sample/sample.npy, sample/design.json).

Population: every record of the production strategy dumps (one record = one distinct sampled strategy of one replicate
at one game; 'count' = its samples in that replicate), both regimes: F1 (N = 1000) data/F1/S_1o3_<task>_1.bin and
F2 (N = 100) data/F2/S_1o3_<task>_1.bin, task = games.dat line, 512 games each, 40 000 samples per game.
Strata h = (regime, wedge, production atom).

Each record i is included independently (Poisson sampling) with the exactly known probability
    pi_i = 1 - (1 - piG_i)(1 - piA_h)(1 - piP_i),
the union of three components:
    G  population share: piG_i = min(1, bG_w c_i / C_w),   C_w = all samples of the regime's games of wedge w
       (40 000 x games of w) -- random individuals of the population of each wedge (probability proportional to the
       samples of the strategy), bG_w = 3000 in W, N and E and 1000 in S;
    A  per strategy:     piA_h = min(1, aH / N_h),          N_h = records of the stratum (equal probability);
    P  per sample:       piP_i = min(1, bH c_i / C_h),      C_h = samples of the stratum,
with aH = bH = OVER for the production atoms 110, 111, 011, 001 (oversampled) and BASE otherwise.
A strategy's share of its game's population is c_i / 40 000, so the Horvitz-Thompson estimate of the population share
of records with property y (pooled over the games of a regime or a wedge, every game weighted equally) is
    sum_{i in sample} y_i (c_i / C) / pi_i,   C = 40 000 x (number of games).
"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); KIT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import load

SEED = 20260927
BG = {'W': 3000, 'S': 1000, 'N': 3000, 'E': 3000}   # expected G draws per regime and wedge
OVER, BASE = 150, 50           # aH = bH per stratum
OVERSAMPLED = {6, 7, 3, 1}     # 110, 111, 011, 001
WEDGE_ID = {'W': 0, 'S': 1, 'N': 2, 'E': 3}
SPR = 40000                    # samples per game

SDT = np.dtype([('sid', '<i8'), ('reg', '<i1'), ('task', '<i4'), ('row', '<i4'), ('irep', '<i4'), ('count', '<i4'),
                ('atom', '<i4'), ('wedge', '<i1'), ('stratum', '<i4'), ('p', '<f8', (16,)),
                ('pi', '<f8'), ('piG', '<f8'), ('piA', '<f8'), ('piP', '<f8')])

def stratum_id(reg, wedge, atom): return (reg - 1)*32 + wedge*8 + atom

def population():
    """every record of both regimes with its stratum and inclusion probabilities (deterministic)"""
    G = load.games(os.path.join(KIT, 'games.dat'))
    cols = {k: [] for k in ('reg', 'task', 'row', 'irep', 'count', 'atom', 'wedge')}
    for reg in (1, 2):
        for task in range(1, 513):
            d = load.dump('F%d' % reg, task); n = len(d)
            assert d['count'].sum() == SPR, (reg, task)
            assert ((d['atom'] >= 0) & (d['atom'] <= 7) & (d['atom'] != 5)).all()
            cols['reg'].append(np.full(n, reg, np.int8)); cols['task'].append(np.full(n, task, np.int32))
            cols['row'].append(np.arange(n, dtype=np.int32)); cols['irep'].append(d['irep'].astype(np.int32))
            cols['count'].append(d['count'].astype(np.int32)); cols['atom'].append(d['atom'].astype(np.int32))
            cols['wedge'].append(np.full(n, WEDGE_ID[G[task - 1]['wedge']], np.int8))
    A = {k: np.concatenate(v) for k, v in cols.items()}
    n = len(A['reg']); h = stratum_id(A['reg'].astype(np.int64), A['wedge'].astype(np.int64), A['atom'].astype(np.int64))
    c = A['count'].astype(np.float64)
    Cr = {1: 512*SPR, 2: 512*SPR}
    for r in (1, 2): assert c[A['reg'] == r].sum() == Cr[r]
    hs = np.unique(h); Nh = {int(k): int((h == k).sum()) for k in hs}; Ch = {int(k): int(c[h == k].sum()) for k in hs}
    aH = np.array([OVER if (k % 8) in OVERSAMPLED else BASE for k in h], dtype=np.float64)
    ngw = {w: sum(1 for g in G if g['wedge'] == w) for w in 'WSNE'}
    bgA = np.array([BG['WSNE'[int(w)]] for w in A['wedge']], dtype=np.float64)
    CwA = np.array([SPR*ngw['WSNE'[int(w)]] for w in A['wedge']], dtype=np.float64)
    piG = np.minimum(1.0, bgA*c/CwA)
    NhA = np.array([Nh[int(k)] for k in h], dtype=np.float64); ChA = np.array([Ch[int(k)] for k in h], dtype=np.float64)
    piA = np.minimum(1.0, aH/NhA); piP = np.minimum(1.0, aH*c/ChA)
    pi = 1.0 - (1.0 - piG)*(1.0 - piA)*(1.0 - piP)
    return G, A, h, c, Cr, hs, Nh, Ch, piG, piA, piP, pi

def main(out=os.path.join(KIT, 'sample')):
    G, A, h, c, Cr, hs, Nh, Ch, piG, piA, piP, pi = population(); n = len(h)
    rng = np.random.default_rng(SEED); U = rng.random(n)
    sel = np.where(U < pi)[0]
    S = np.zeros(len(sel), dtype=SDT)
    S['sid'] = np.arange(len(sel))
    for k in ('reg', 'task', 'row', 'irep', 'count', 'atom', 'wedge'): S[k] = A[k][sel]
    S['stratum'] = h[sel]; S['pi'] = pi[sel]; S['piG'] = piG[sel]; S['piA'] = piA[sel]; S['piP'] = piP[sel]
    # the strategies themselves, from the dumps
    for reg in (1, 2):
        for task in np.unique(S['task'][S['reg'] == reg]):
            m = np.where((S['reg'] == reg) & (S['task'] == task))[0]
            d = load.dump('F%d' % reg, int(task))
            assert (d['atom'][S['row'][m]] == S['atom'][m]).all() and (d['count'][S['row'][m]] == S['count'][m]).all()
            S['p'][m] = d['p'][S['row'][m]]
    os.makedirs(out, exist_ok=True)
    np.save(os.path.join(out, 'sample.npy'), S)
    strata = {}
    for k in hs:
        k = int(k); m = h == k; ms = S['stratum'] == k
        strata[str(k)] = dict(reg='F%d' % (k // 32 + 1), wedge='WSNE'[(k % 32) // 8], atom=format(k % 8, '03b'),
                              N=Nh[k], C=Ch[k], expected=float(pi[m].sum()), drawn=int(ms.sum()),
                              drawn_samples=int(S['count'][ms].sum()))
    design = dict(seed=SEED, BG=BG, OVER=OVER, BASE=BASE, oversampled=['110', '111', '011', '001'],
                  samples_per_game=SPR, games_per_wedge={w: sum(1 for g in G if g['wedge'] == w) for w in 'WSNE'},
                  records=int(n), records_by_regime={('F%d' % r): int((A['reg'] == r).sum()) for r in (1, 2)},
                  expected_size=float(pi.sum()), drawn=int(len(sel)),
                  drawn_by_regime={('F%d' % r): int((S['reg'] == r).sum()) for r in (1, 2)},
                  strata=strata, doc=__doc__)
    json.dump(design, open(os.path.join(out, 'design.json'), 'w'), indent=1)
    print('records %d, expected %.1f, drawn %d (F1 %d, F2 %d)' % (n, pi.sum(), len(sel), (S['reg'] == 1).sum(), (S['reg'] == 2).sum()))
    for k in sorted(strata, key=int):
        s = strata[k]; print('  %s %s %s  N=%8d C=%9d  expected %7.1f drawn %5d' % (s['reg'], s['wedge'], s['atom'], s['N'], s['C'], s['expected'], s['drawn']))

if __name__ == '__main__':
    main()
