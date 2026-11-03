"""test_exact_full.py -- brute-force reference: for a few sampled strategies, the exact margin and the exact payoff
difference against ALL 65536 binary co-players (exact rationals only, no interval, no structural rule), compared with
the audit's bits and counts.  About 25 s per strategy on one core.
Usage: python3 test_exact_full.py <n> [seed]"""
import os, sys, time, numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); KIT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(KIT, 'lib'))
import audit as au, load
from cert import engine as en, exact as ex

def reference(p, cs, ct):
    P = [Fraction(float(x)) for x in p]; S, T = Fraction(cs), Fraction(ct); sign = 1 if T > S else -1
    below = S + T < 2
    Ls = ex.lead_exact(P, P); W = Ls['W']; pp = W[0] + S*W[1] + T*W[2]
    E = int(en.meets(en.support_class(Ls['supp']), 'C' if below else 'A'))
    nneg = nzero = ngain = ntie = 0
    for q in range(65536):
        Wq = ex.lead_exact(P, ex.g2p(q))['W']
        m = sign*(Wq[2] - Wq[1]); g = Wq[0] + T*Wq[1] + S*Wq[2] - pp
        nneg += m < 0; nzero += m == 0; ngain += g > 0; ntie += g == 0
    return E, int(ngain == 0), int(nneg == 0), dict(nneg=nneg, nzero=nzero, ngain=ngain, ntie=ntie)

def main(n=4, seed=3):
    S = np.load(os.path.join(KIT, 'sample', 'sample.npy')); G = load.games(os.path.join(KIT, 'games.dat'))
    rng = np.random.default_rng(seed)
    want = [6, 7, 3, 1, 2, 6, 7, 0]
    for k in range(n):
        idx = np.where(S['atom'] == want[k % len(want)])[0]; i = int(rng.choice(idx)); r = S[i]; g = G[int(r['task']) - 1]
        t0 = time.time(); E, Sb, C, cnt = reference(r['p'], g['cs'], g['ct']); t1 = time.time() - t0
        o = au.audit(r['p'], g['cs'], g['ct'], prod_atom=int(r['atom']))
        ref = '%d%d%d' % (E, Sb, C)
        # audit counts that are complete: C counts when the audit resolved every margin, S counts likewise
        c_neg = o['C']['n_neg_i'] + o['C']['n_neg_x']; s_gain = o['S']['n_gain_i'] + o['S']['n_struct_gain'] + o['S']['n_gain_x']
        okC = (o['C']['n_exact'] < o['C']['n_unres']) or (c_neg == cnt['nneg'] and
               o['C']['n_struct0'] + o['C']['n_zero_x'] == cnt['nzero'])
        okS = (o['S']['n_exact'] < o['S']['n_unres']) or (s_gain == cnt['ngain'] and
               o['S']['n_struct_tie'] + o['S']['n_tie_x'] == cnt['ntie'])
        print('sid %d F%d %s task %d prod %s audit %s reference %s  counts ok C %s S %s  ref %s  (%.0f s)' %
              (r['sid'], r['reg'], g['wedge'], r['task'], o['prod'], o['val'], ref, okC, okS, cnt, t1), flush=True)

if __name__ == '__main__':
    main(*(int(x) for x in sys.argv[1:]))
