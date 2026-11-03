"""audit.py -- direct interval / exact audit of the production atoms of stored memory-two evolutionary strategies.

For a stored strategy p (the doubles of the production dump, taken as exact rationals) at its game
(R, S, T, P) = (1, cs, ct, 0) of games.dat (the exact doubles), in the limit eps -> 0:

  E  efficient: the exact support class of the self-play limit (exact rationals, ex.lead_exact(p, p); cross-checked
     against en.efficiency_class of p's support stratum, lemma L1), required class 'C' (w_CC = 1) below the switch line
     S + T < 2 and 'A' (w_CD + w_DC = 1) above it (no game of games.dat lies on the line or on T = S);
  C  competitive: sign*(w_DC - w_CD) >= 0 against every binary co-player q, sign = +1 iff T > S;
  S  stable: pi(q, p) <= pi(p, p) against every binary co-player q, with pi(p, p) = R w_CC + S w_CD + T w_DC of the
     exact self-play limit and pi(q, p) = R w_CC + T w_CD + S w_DC + P w_DD in the focal-p frequencies of the (p, q)
     chain (the (q, p) chain is its mirror image, state j <-> mirror(j)).
Weak inequalities: a tie is not a gain and satisfies competitiveness.

Per co-player the frequencies come from one outward-rounded interval evaluation (en.i_scanW); the margin enclosure is
formed exactly as certeng.c's imargin does, the payoff-difference enclosure with one outward step (nextafter) after
every correctly rounded operation.  A margin is a structural zero when no support state has recent outcome CD or DC;
a payoff is structural when the limit support holds a single outcome (then pi(q, p) is that outcome's payoff exactly).
Every enclosure that contains 0 is decided by exact rationals (ex.lead_exact).  Nothing is decided by a tolerance.

The production classifier (src/prodcls.c + production's m2f.c) is re-run on the same doubles for the reproduction
check and to name the co-players that decided each production bit.

Usage:
  python3 audit.py run <sample.npy> <chunk> <nchunks> <outdir>     one chunk (rows i with i % nchunks == chunk);
                                                                    idempotent: records already in the output are skipped
  python3 audit.py one <regime F1|F2> <task> <row>                  one dump record, printed
"""
import os, sys, json, time, math, ctypes
import numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); KIT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(KIT, 'lib')); sys.path.insert(0, HERE)
from cert import engine as en, exact as ex
import load

NTH = en.MAX_THREADS
_P = ctypes.CDLL(os.path.join(HERE, 'libprod.dylib' if sys.platform == 'darwin' else 'libprod.so'))
_dp = ctypes.POINTER(ctypes.c_double); _lp = ctypes.POINTER(ctypes.c_long)
_P.prod_classify.argtypes = [_dp, ctypes.c_double, ctypes.c_double, _dp, _lp]; _P.prod_classify.restype = ctypes.c_int
_P.prod_pair.argtypes = [_dp, ctypes.c_double, ctypes.c_double, ctypes.c_long, _dp]

def prod_classify(p, cs, ct):
    p = np.ascontiguousarray(p, dtype=np.float64); d = np.zeros(8); i = np.zeros(13, dtype=np.int64)
    assert ctypes.sizeof(ctypes.c_long) == 8
    _P.prod_classify(p.ctypes.data_as(_dp), cs, ct, d.ctypes.data_as(_dp), i.ctypes.data_as(_lp))
    return dict(atom=int(i[0]), eff=int(i[1]), riv=int(i[2]), nash=int(i[3]), nash_scan=int(i[4]),
                W=[float(x) for x in d[:4]], pp=float(d[4]), mmin=float(d[5]), gmax=float(d[6]),
                q_mmin=int(i[5]), q_gmax=int(i[6]), n_neg=int(i[7]), n_neg_tol=int(i[8]), n_gain=int(i[9]),
                n_gain_tol=int(i[10]), q_first_neg=int(i[11]), q_first_gain=int(i[12]))

def prod_pair(p, cs, ct, g):
    p = np.ascontiguousarray(p, dtype=np.float64); d = np.zeros(7)
    _P.prod_pair(p.ctypes.data_as(_dp), cs, ct, int(g), d.ctypes.data_as(_dp))
    return dict(margin=float(d[0]), payq=float(d[1]), pp=float(d[2]))

def dn(x): return np.nextafter(x, -np.inf)
def up(x): return np.nextafter(x, np.inf)
def scal(c, lo, hi):
    """enclosure of c*[lo, hi] for a double constant c"""
    if c == 1.0: return lo, hi
    if c >= 0: return dn(c*lo), up(c*hi)
    return dn(c*hi), up(c*lo)

def fr(x):
    x = Fraction(x); return '%d/%d' % (x.numerator, x.denominator)

def atom_str(a): return format(int(a), '03b')

def audit(p, cs, ct, prod_atom=None, want_prod=True):
    """validated atom of the strategy p (16 doubles) at the game (1, cs, ct, 0); see the module docstring"""
    t0 = time.time()
    p = np.asarray(p, dtype=np.float64)
    assert all((x == 0.0 or x == 1.0 or 0.0 < x < 1.0) for x in p)
    R, S, T = Fraction(1), Fraction(cs), Fraction(ct)
    below = (S + T < 2); sign = 1 if T > S else -1; req = 'C' if below else 'A'
    assert S + T != 2 and T != S
    P = [Fraction(float(x)) for x in p]
    # ---------------------------------------------------------------- E: exact self-play limit
    Ls = ex.lead_exact(P, P)
    ecls = en.support_class(Ls['supp']); ecls_L1 = en.efficiency_class(en.pattern_of(p))
    E = int(en.meets(ecls, req))
    Wself = Ls['W']; pp = R*Wself[0] + S*Wself[1] + T*Wself[2]
    # ---------------------------------------------------------------- one interval scan over the 65536 co-players
    plo, phi = en.point_box(P)
    Z = en.i_scanW(plo, phi, nthreads=NTH)
    Wlo, Whi, supp = Z['Wlo'], Z['Whi'], Z['supp'].astype(np.int64)
    onCDDC = (supp & 0x0FF0) != 0                       # states 4..11: recent outcome CD or DC
    if sign > 0: mlo = dn(Wlo[:, 2] - Whi[:, 1]); mhi = up(Whi[:, 2] - Wlo[:, 1])
    else:        mlo = dn(Wlo[:, 1] - Whi[:, 2]); mhi = up(Whi[:, 1] - Wlo[:, 2])
    mcls = np.where(~onCDDC, 0, np.where(mlo > 0, 1, np.where(mhi < 0, -1, 2)))
    mlo = np.where(onCDDC, mlo, 0.0); mhi = np.where(onCDDC, mhi, 0.0)
    # pi(q, p) = R w_CC + T w_CD + S w_DC  (+ P w_DD, P = 0)
    tlo, thi = scal(float(ct), Wlo[:, 1], Whi[:, 1]); slo, shi = scal(float(cs), Wlo[:, 2], Whi[:, 2])
    qlo = dn(dn(Wlo[:, 0] + tlo) + slo); qhi = up(up(Whi[:, 0] + thi) + shi)
    pplo, pphi = en.enclose(pp)
    glo = dn(qlo - pphi); ghi = up(qhi - pplo)
    gcls = np.where(glo > 0, 1, np.where(ghi <= 0, -1, 2))      # +1 gain, -1 no gain (gain <= 0), 2 unresolved
    # structural payoffs: a single outcome in the support
    pres = np.stack([((supp >> (4*o)) & 0xF) != 0 for o in range(4)], 1)
    single = pres.sum(1) == 1
    outpay = [R, T, S, Fraction(0)]                 # q's payoff at outcome CC, CD, DC, DD (focal p first)
    gstruct = np.zeros(en.NQ, dtype=np.int8)       # 0 not structural, else 1 gain / 2 tie / 3 loss
    for o in range(4):
        m = single & pres[:, o]
        if m.any():
            d = outpay[o] - pp; gstruct[m] = 1 if d > 0 else (2 if d == 0 else 3)
            # sanity: an interval that already decides must agree with the structural value
            if d > 0: assert not (gcls[m] == -1).any()
            if d < 0: assert not (gcls[m] == 1).any()
            if d == 0: assert not (gcls[m] == 1).any() and not ((ghi[m] < 0).any())
    gcls = np.where(gstruct == 1, 1, np.where(gstruct == 2, 0, np.where(gstruct == 3, -1, gcls)))
    # ---------------------------------------------------------------- exact rationals where an enclosure contains 0
    # an unresolved enclosure is computed exactly only when the bit depends on it: C when no margin is certified
    # negative, S when no gain is certified positive (otherwise the bit is already decided and the rest is skipped)
    mun = np.where(mcls == 2)[0]; gun = np.where(gcls == 2)[0]
    needC = not (mcls == -1).any(); needS = not (gcls == 1).any()
    need = sorted(set(mun.tolist() if needC else []) | set(gun.tolist() if needS else []))
    mx = {}; gx = {}
    for q in need:
        W = ex.lead_exact(P, ex.g2p(int(q)))['W']
        if mcls[q] == 2 and needC: mx[int(q)] = sign*(W[2] - W[1])
        if gcls[q] == 2 and needS: gx[int(q)] = R*W[0] + T*W[1] + S*W[2] - pp
    mneg_x = [q for q, v in mx.items() if v < 0]; mzero_x = sum(1 for v in mx.values() if v == 0)
    gpos_x = [q for q, v in gx.items() if v > 0]; gzero_x = sum(1 for v in gx.values() if v == 0)
    n_neg_i = int((mcls == -1).sum()); n_gain_i = int((gcls == 1).sum()) - int((gstruct == 1).sum())
    C = int(n_neg_i == 0 and not mneg_x)
    Sb = int(n_gain_i == 0 and not gpos_x and not (gstruct == 1).any())
    atom = 4*E + 2*Sb + C
    # ---------------------------------------------------------------- extremes, exact at the interval arg-extreme
    def exact_W(q): return ex.lead_exact(P, ex.g2p(int(q)))['W']
    qm = int(np.argmin(mlo)); Wm = exact_W(qm); mmin_x = sign*(Wm[2] - Wm[1])
    if mneg_x:
        qq = min(mneg_x, key=lambda q: mx[q])
        if mx[qq] < mmin_x: qm, mmin_x = qq, mx[qq]
    gup = ghi.copy()
    for o in range(4):
        m = single & pres[:, o]
        if m.any(): gup[m] = float(outpay[o] - pp)
    qg = int(np.argmax(gup)); Wg = exact_W(qg); gmax_x = R*Wg[0] + T*Wg[1] + S*Wg[2] - pp
    if gpos_x:
        qq = max(gpos_x, key=lambda q: gx[q])
        if gx[qq] > gmax_x: qg, gmax_x = qq, gx[qq]
    how_C = ('interval' if n_neg_i else 'exact') if not C else ('interval+structural' if not mx else 'exact')
    how_S = ('interval' if n_gain_i else ('structural' if (gstruct == 1).any() else 'exact')) if not Sb else \
            ('interval+structural' if not gx else 'exact')
    out = dict(
        game=dict(cs=cs, ct=ct, below=bool(below), sign=sign, req=req),
        val=atom_str(atom),
        E=dict(val=E, cls=ecls, cls_L1=ecls_L1, supp=int(Ls['supp']), w_self=[fr(x) for x in Wself], how='exact support'),
        C=dict(val=C, how=how_C, n_pos_i=int((mcls == 1).sum()), n_neg_i=n_neg_i, n_struct0=int((mcls == 0).sum()),
               n_unres=int(len(mun)), n_exact=len(mx), n_pos_x=len(mx) - mzero_x - len(mneg_x), n_zero_x=mzero_x, n_neg_x=len(mneg_x),
               q_min=qm, min=fr(mmin_x), min_f=float(mmin_x)),
        S=dict(val=Sb, how=how_S, n_loss_i=int((gcls == -1).sum()) - int((gstruct == 3).sum()), n_gain_i=n_gain_i,
               n_struct_gain=int((gstruct == 1).sum()), n_struct_tie=int((gstruct == 2).sum()),
               n_struct_loss=int((gstruct == 3).sum()), n_unres=int(len(gun)), n_exact=len(gx), n_gain_x=len(gpos_x), n_tie_x=gzero_x,
               n_loss_x=len(gx) - gzero_x - len(gpos_x), q_max=qg, max=fr(gmax_x), max_f=float(gmax_x),
               pp=fr(pp)),
        n_exact_calls=int(len(need)) + 3)
    if ecls != ecls_L1: out['E']['L1_mismatch'] = True
    if want_prod:
        pr = prod_classify(p, cs, ct); out['prod_repro'] = atom_str(pr['atom']); out['pr'] = pr
        if prod_atom is not None:
            out['prod'] = atom_str(prod_atom); pa = int(prod_atom)
            bits = {'E': (pa >> 2) & 1, 'S': (pa >> 1) & 1, 'C': pa & 1}
            ch = [b for b in 'ESC' if bits[b] != out[b]['val']]
            out['changed'] = ''.join(ch)
            resp = {}
            if 'E' in ch: resp['E'] = dict(exact_cls=ecls, prod_W=pr['W'])
            if 'C' in ch:
                if bits['C'] == 1:            # production accepted, a margin is negative exactly
                    resp['C'] = dict(q=qm, exact=fr(mmin_x), exact_f=float(mmin_x), prod=prod_pair(p, cs, ct, qm)['margin'])
                else:                          # production rejected at a q whose margin is >= 0 exactly
                    qs = sorted({pr['q_first_neg'], pr['q_mmin']} - {-1})
                    resp['C'] = [dict(q=q, prod=prod_pair(p, cs, ct, q)['margin'],
                                      exact=fr(sign*(exact_W(q)[2] - exact_W(q)[1]))) for q in qs]
            if 'S' in ch:
                if bits['S'] == 1:            # production had no gain (or skipped the scan by the theorem), exact gain > 0
                    pq = prod_pair(p, cs, ct, qg)
                    resp['S'] = dict(q=qg, exact=fr(gmax_x), exact_f=float(gmax_x), prod=pq['payq'] - pq['pp'],
                                     prod_by_theorem=bool(pr['eff'] and pr['riv']))
                else:
                    qs = sorted({pr['q_first_gain'], pr['q_gmax']} - {-1}); L = []
                    for q in qs:
                        W = exact_W(q); pq = prod_pair(p, cs, ct, q)
                        L.append(dict(q=q, prod=pq['payq'] - pq['pp'], exact=fr(R*W[0] + T*W[1] + S*W[2] - pp)))
                    resp['S'] = L
            if resp: out['responsible'] = resp
    out['sec'] = round(time.time() - t0, 4)
    return out

# ---------------------------------------------------------------------------------------------- chunk runner
def run_chunk(samplefile, chunk, nchunks, outdir):
    S = np.load(samplefile)
    G = load.games(os.path.join(KIT, 'games.dat'))
    rows = np.arange(len(S))[np.arange(len(S)) % nchunks == chunk]
    os.makedirs(outdir, exist_ok=True)
    fn = os.path.join(outdir, 'chunk_%04d.jsonl' % chunk)
    done = set()
    if os.path.exists(fn):
        good = []
        for line in open(fn):
            try: r = json.loads(line); done.add(int(r['sid'])); good.append(line)
            except Exception: pass                      # a line cut by a kill: dropped and redone
        with open(fn + '.tmp', 'w') as f: f.writelines(good)
        os.replace(fn + '.tmp', fn)
    t0 = time.time(); n = 0
    with open(fn, 'a') as f:
        for i in rows:
            r = S[i]; sid = int(r['sid'])
            if sid in done: continue
            g = G[int(r['task']) - 1]
            o = audit(r['p'], g['cs'], g['ct'], prod_atom=int(r['atom']))
            o.update(sid=sid, reg='F%d' % int(r['reg']), task=int(r['task']), gid=g['gid'], ipt=g['ipt'], wedge=g['wedge'],
                     irep=int(r['irep']), row=int(r['row']), count=int(r['count']), stratum=int(r['stratum']),
                     pi=float(r['pi']))
            f.write(json.dumps(o, separators=(',', ':')) + '\n'); f.flush(); n += 1
    todo = len(rows); have = len(done) + n
    with open(os.path.join(outdir, 'chunk_%04d.done' % chunk), 'w') as f:
        f.write('%d of %d records, %d new in %.1f s, threads %d\n' % (have, todo, n, time.time() - t0, NTH))
    print('chunk %d: %d of %d records, %d new in %.1f s' % (chunk, have, todo, n, time.time() - t0))

if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run_chunk(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    elif sys.argv[1] == 'one':
        reg, task, row = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
        d = load.dump(reg, task); G = load.games(os.path.join(KIT, 'games.dat')); g = G[task - 1]
        print(json.dumps(audit(d['p'][row], g['cs'], g['ct'], prod_atom=int(d['atom'][row])), indent=1))
