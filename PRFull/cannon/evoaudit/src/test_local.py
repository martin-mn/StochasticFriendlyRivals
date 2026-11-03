"""test_local.py -- local validation of the audit on a handful of sampled strategies (< 15 min on 3 threads).

  (1) the margin classes formed from i_scanW equal the library's own en.i_scan, co-player by co-player;
  (2) the mirror identity used for pi(q, p): for random binary q, R w_CC + T w_CD + S w_DC of the (p, q) chain equals,
      exactly, q's payoff R w'_CC + S w'_CD + T w'_DC in the (q, p) chain (ex.lead_exact both ways);
  (3) every margin and payoff-difference class decided by an interval or a structural rule agrees with the exact sign
      at random co-players (and at every co-player the audit sent to exact arithmetic, by construction);
  (4) the audit of one record per stratum (up to NPER strata): production atom, reproduced production atom, validated
      atom, time per record.
Usage: python3 test_local.py [nstrata] [seed]
"""
import os, sys, time, json, numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); KIT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(KIT, 'lib'))
import audit as au, load
from cert import engine as en, exact as ex

def check_point(p, cs, ct, nq=60, rng=None):
    P = [Fraction(float(x)) for x in p]; S, T = Fraction(cs), Fraction(ct); sign = 1 if T > S else -1
    plo, phi = en.point_box(P)
    Z = en.i_scanW(plo, phi, nthreads=au.NTH); I = en.i_scan(plo, phi, sign, nthreads=au.NTH)
    supp = Z['supp'].astype(np.int64); on = (supp & 0x0FF0) != 0
    Wlo, Whi = Z['Wlo'], Z['Whi']
    if sign > 0: mlo = au.dn(Wlo[:, 2] - Whi[:, 1]); mhi = au.up(Whi[:, 2] - Wlo[:, 1])
    else: mlo = au.dn(Wlo[:, 1] - Whi[:, 2]); mhi = au.up(Whi[:, 1] - Wlo[:, 2])
    mcls = np.where(~on, 0, np.where(mlo > 0, 1, np.where(mhi < 0, -1, 2)))
    ok1 = bool((mcls == I['cls']).all())
    Ls = ex.lead_exact(P, P); pp = Ls['W'][0] + S*Ls['W'][1] + T*Ls['W'][2]
    tlo, thi = au.scal(float(ct), Wlo[:, 1], Whi[:, 1]); slo, shi = au.scal(float(cs), Wlo[:, 2], Whi[:, 2])
    qlo = au.dn(au.dn(Wlo[:, 0] + tlo) + slo); qhi = au.up(au.up(Whi[:, 0] + thi) + shi)
    pplo, pphi = en.enclose(pp); glo = au.dn(qlo - pphi); ghi = au.up(qhi - pplo)
    qs = rng.choice(65536, nq, replace=False)
    ok2 = ok3 = True
    for q in qs:
        q = int(q); W = ex.lead_exact(P, ex.g2p(q))['W']; Wq = ex.lead_exact(ex.g2p(q), P)['W']
        a = W[0] + T*W[1] + S*W[2]; b = Wq[0] + S*Wq[1] + T*Wq[2]
        if a != b: ok2 = False
        m = sign*(W[2] - W[1]); g = a - pp
        if mcls[q] == 1 and not m > 0: ok3 = False
        if mcls[q] == -1 and not m < 0: ok3 = False
        if mcls[q] == 0 and m != 0: ok3 = False
        if glo[q] > 0 and not g > 0: ok3 = False
        if ghi[q] <= 0 and not g <= 0: ok3 = False
        if not (glo[q] <= float(g) <= ghi[q]): ok3 = False
    return ok1, ok2, ok3

def main(nstrata=24, seed=1):
    S = np.load(os.path.join(KIT, 'sample', 'sample.npy')); G = load.games(os.path.join(KIT, 'games.dat'))
    rng = np.random.default_rng(seed)
    hs = np.unique(S['stratum']); rng.shuffle(hs); hs = hs[:nstrata]
    tot = 0; bad = 0
    print('reg wedge task  prod repro  val  changed  sec   nexact  checks(i_scan, mirror, signs)  tol(neg,gain)')
    for h in sorted(hs):
        idx = np.where(S['stratum'] == h)[0]; i = int(rng.choice(idx)); r = S[i]; g = G[int(r['task']) - 1]
        t0 = time.time()
        o = au.audit(r['p'], g['cs'], g['ct'], prod_atom=int(r['atom']))
        chk = check_point(r['p'], g['cs'], g['ct'], rng=rng)
        tot += o['sec']
        flag = '' if (o['prod'] == o['prod_repro'] and all(chk)) else '  <--'
        if flag: bad += 1
        print('F%d  %s  %4d  %s  %s  %s  %-4s %6.2f %6d   %s  (%d,%d)%s' % (r['reg'], g['wedge'], r['task'], o['prod'], o['prod_repro'],
              o['val'], o['changed'] or '-', o['sec'], o['n_exact_calls'], chk, o['pr']['n_neg_tol'], o['pr']['n_gain_tol'], flag))
        if o.get('responsible'): print('     responsible:', json.dumps(o['responsible']))
    print('audit time %.1f s for %d records (%.2f s each, %d threads); flagged %d' % (tot, len(hs), tot/len(hs), au.NTH, bad))

if __name__ == '__main__':
    main(*(int(x) for x in sys.argv[1:]))
