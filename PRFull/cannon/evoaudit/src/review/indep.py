"""Independent recomputation of the atom of probe strategies with the certificate library (not audit.py).
E: exact self-play frequencies (W_CC == 1 below the switch line, W_CD + W_DC == 1 above).
C: ce.certify_point (i_scan + exact margins on unresolved), sign from exact T vs S.
S: q-FOCAL chain (no mirror identity): i_eval(q, q, p) enclosures of pi(q,p) = R W'_CC + S W'_CD + T W'_DC,
   outward-rounded with Fractions at the ends (enclosure endpoints are floats; combine in exact rationals),
   exact lead_exact(q, P) wherever the enclosure of pi(q,p) - pi(p,p) is not strictly signed or touches 0.
mode 'full': in addition exact rationals for ALL 65536 q, margin from lead_exact(P, q), pi(q,p) from lead_exact(q, P)."""
import sys, os, json, time, numpy as np
from fractions import Fraction
from multiprocessing import Pool
sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/check')
from cert import engine as en, exact as ex, certify as ce
KIT = '/Users/martin/Documents/CL1/PRFull/cannon/evoaudit'

def games():
    G = []
    for line in open(KIT + '/games.dat'):
        if line.startswith('#'): continue
        t = line.split()
        if len(t) < 8: continue
        G.append(dict(R=Fraction(t[2]), S=Fraction(float(t[3])), T=Fraction(float(t[4])), tag=t[7]))
    return G

def full_job(args):
    Pl, q = args
    P = [Fraction(x) for x in Pl]; Q = ex.g2p(q)
    Wpq = ex.lead_exact(P, Q)['W']; Wqp = ex.lead_exact(Q, P)['W']
    return q, Wpq, Wqp

def check(sid, full=False, nproc=3):
    S = np.load(KIT + '/sample/sample.npy'); r = S[sid]; assert int(r['sid']) == sid
    g = games()[int(r['task']) - 1]; R, Sg, T = g['R'], g['S'], g['T']
    p = np.array(r['p'], dtype=float); P = [Fraction(float(x)) for x in p]
    below = Sg + T < 2*R; sign = 1 if T > Sg else -1
    wedge = {(True, 1): 'W', (True, -1): 'S', (False, 1): 'N', (False, -1): 'E'}[(below, sign)]
    t0 = time.time()
    Ls = ex.lead_exact(P, P); Ws = Ls['W']
    E = int(Ws[0] == 1) if below else int(Ws[1] + Ws[2] == 1)
    pp = R*Ws[0] + Sg*Ws[1] + T*Ws[2]
    # C by the library's point certificate
    cp = ce.certify_point(P, sign=sign, nthreads=3, nproc=1)
    C = int(cp['n_neg'] == 0)
    # S: q-focal interval + exact
    n_gain = n_tie = n_loss = n_exact = 0; gmax = None
    for q in range(65536):
        Q = en.g2p(q)
        z = en.i_eval(Q, Q, p)                       # focal q (binary), co-player p (point)
        # NB: q's payoff at q-focal outcome CD (q cooperates, p defects) is S, at DC it is T
        lo2 = R*Fraction(z['Wlo'][0]) + (Sg*Fraction(z['Wlo'][1]) if Sg >= 0 else Sg*Fraction(z['Whi'][1])) + \
              (T*Fraction(z['Wlo'][2]) if T >= 0 else T*Fraction(z['Whi'][2]))
        hi2 = R*Fraction(z['Whi'][0]) + (Sg*Fraction(z['Whi'][1]) if Sg >= 0 else Sg*Fraction(z['Wlo'][1])) + \
              (T*Fraction(z['Whi'][2]) if T >= 0 else T*Fraction(z['Wlo'][2]))
        glo, ghi = lo2 - pp, hi2 - pp
        if glo > 0: n_gain += 1; val = None
        elif ghi < 0: n_loss += 1; val = None
        else:
            W2 = ex.lead_exact(Q, P)['W']; val = R*W2[0] + Sg*W2[1] + T*W2[2] - pp; n_exact += 1
            if val > 0: n_gain += 1
            elif val == 0: n_tie += 1
            else: n_loss += 1
        est = val if val is not None else (glo + ghi)/2
        if gmax is None or est > gmax[0]: gmax = (est, q)
    Sb = int(n_gain == 0)
    out = dict(sid=sid, reg='F%d' % r['reg'], wedge=wedge, tagwedge=g['tag'][-1], task=int(r['task']), prod=format(int(r['atom']), '03b'),
               indep='%d%d%d' % (E, Sb, C), E_selfW=[str(x) for x in Ws], C=dict(n_neg=cp['n_neg'], n_zero=cp['n_zero'],
               n_struct=cp['n_struct_zero'], n_exact_zero=cp['n_exact_zero'], n_unres_int=cp['n_unresolved_by_interval'],
               min=str(cp['min_margin']), min_f=float(cp['min_margin'])),
               S=dict(n_gain=n_gain, n_tie=n_tie, n_loss=n_loss, n_exact=n_exact, gmax_f=float(gmax[0]), q_gmax=gmax[1]),
               sec=round(time.time() - t0, 1))
    if full:
        t1 = time.time(); Pl = [str(x) for x in P]
        with Pool(nproc) as pool:
            res = pool.map(full_job, [(Pl, q) for q in range(65536)], chunksize=512)
        nneg = nzero = ngain = ntie = 0; mmin = None; gmx = None; mirror_bad = 0
        for q, Wpq, Wqp in res:
            m = sign*(Wpq[2] - Wpq[1]); gq = R*Wqp[0] + Sg*Wqp[1] + T*Wqp[2] - pp
            if not (Wqp[0] == Wpq[0] and Wqp[1] == Wpq[2] and Wqp[2] == Wpq[1] and Wqp[3] == Wpq[3]): mirror_bad += 1
            nneg += m < 0; nzero += m == 0; ngain += gq > 0; ntie += gq == 0
            if mmin is None or m < mmin: mmin = m
            if gmx is None or gq > gmx: gmx = gq
        out['full'] = dict(atom='%d%d%d' % (E, int(ngain == 0), int(nneg == 0)), nneg=nneg, nzero=nzero, ngain=ngain, ntie=ntie,
                           mmin=float(mmin), gmax=float(gmx), mirror_mismatch=mirror_bad, sec=round(time.time() - t1, 1))
    return out

if __name__ == '__main__':
    full = sys.argv[1] == 'full'
    for s in sys.argv[2:]:
        print(json.dumps(check(int(s), full=full)), flush=True)
