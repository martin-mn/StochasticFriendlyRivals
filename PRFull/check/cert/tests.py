#!/usr/bin/env python3
"""Tests of the certification library.  Run:   cd PRFull/check && nice -n 10 python3 cert/tests.py
(at most 3 cores: C threads and exact-arithmetic pools are capped by engine.MAX_THREADS).

  (a) SI Table 2 / M2f Table 3: exact zero-margin counts and self-play of the W, N, S, E witnesses,
      by the interval scan + exact rationals on the unresolved co-players, and for W and N also by a
      full exact scan of all 65536 co-players (no floating point at all);
  (b) float vs exact (and interval containment, mpq vs Fraction, bit-identity with the frozen M2f
      engine) on 200 random (stratum point, co-player) pairs;
  (c) the L3 corner rule: exact multilinear extrapolation of N_q = D* margin and of every A*_k from
      a 2-level grid to the corners and to random face points, against 'the corner's own margin
      times D* if E_c = E, else 0'; the same at points 1e-9 from a corner; E and D* against the
      exact finite-eps matrix-tree determinants;
  (d) 10^5 random binary pairs exactly: the '1/37 gap' claim (REFUTED: 1/963 occurs) and its
      rigorous replacement, integrality of D* (w_DC - w_CD) (binary_pair_exact == exact);
  (e) the certificate primitives: positive and negative controls and sampled soundness checks;
  (t) timings per call.
"""
import sys, os, time, random, itertools, math
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from fractions import Fraction as F
from multiprocessing import Pool
from cert import engine as en, exact as ex, certify as ce

M2F = '/Users/martin/Documents/CL1/PartnersRivalsM2f/check'
RES = {}
def report(key, ok, **kw):
    RES[key] = dict(ok=bool(ok), **kw)
    print(('PASS ' if ok else 'FAIL ') + key + '  ' + '  '.join('%s=%s' % (k, v) for k, v in kw.items()), flush=True)

def skipped(key, why):
    """a check that cannot run here: said, and counted in the summary (no silent skips)"""
    RES[key] = dict(skipped=True, why=why)
    print('SKIPPED ' + key + '  (' + why + ')', flush=True)

WIT = {
 'W': [1, F(1,2), 1, F(1,2), 0, 0, 1, F(1,2), 1, F(1,4), 1, F(1,4), F(1,4), 0, F(1,4), 0],
 'S': [1, F(1,2), F(1,2), F(1,2), 1, F(1,2), F(1,2), F(1,2), 1, F(1,2), F(1,2), 1, 0, F(1,2), F(1,2), F(1,2)],
 'N': [F(2,3), 1, 0, F(1,3), F(1,3), F(1,3), 0, 0, F(2,3), 1, F(2,3), F(1,3), F(1,3), 0, 1, 0],
 'E': [1, 0, 1, F(2,3), F(2,3), F(1,3), 0, F(1,3), 1, 1, F(2,3), F(2,3), F(2,3), 1, 0, F(1,3)],
}
# SI Table 2 (ms.tex:515-518) and M2f Text/tables.tex: stratum, self-play, #co-players at margin 0
TABLE2 = {'W': ('1*1*001*1*1**0*0', (1, 0, 0, 0), 37996), 'S': ('1***1***1**10***', (1, 0, 0, 0), 30850),
          'N': ('*10***00*1***010', (0, F(1,2), F(1,2), 0), 31024), 'E': ('101***0*11***10*', (0, F(1,2), F(1,2), 0), 31024)}
SPLIT = {'W': (33260, 4736), 'N': (13436, 17588)}      # structural + identities (triage report 1)

# ------------------------------------------------------------------ (a)
def _full_exact_chunk(args):
    name, a, b = args
    sign = en.WEDGES[name][0]; p = WIT[name]; nz = 0; worst = None; nstruct = 0
    for g in range(a, b):
        r = ex.lead_exact(p, ex.g2p(g), as_fraction=False); W = r['W']
        d = sign*(W[2] - W[1])
        if d == 0:
            nz += 1
            if W[1] == 0 and W[2] == 0: nstruct += 1
        if worst is None or d < worst: worst = d
    return nz, nstruct, F(int(worst.numerator), int(worst.denominator))

def test_a():
    t0 = time.time()
    for name in ('W', 'S', 'N', 'E'):
        pat, sp, nzero = TABLE2[name]
        r = ce.certify_point(WIT[name], wedge=name)
        ok = (en.pattern_str(en.pattern_of(WIT[name])) == pat and tuple(r['selfplay']) == tuple(F(x) for x in sp)
              and r['n_zero'] == nzero and r['n_neg'] == 0 and r['min_margin'] == 0 and r['certified'])
        if name in SPLIT: ok = ok and (r['n_struct_zero'], r['n_exact_zero']) == SPLIT[name]
        report('a.%s interval+exact' % name, ok, stratum=pat, selfplay='(' + ','.join(map(str, r['selfplay'])) + ')',
               zero=r['n_zero'], struct=r['n_struct_zero'], identities=r['n_exact_zero'], pos=r['n_pos'],
               neg=r['n_neg'], sec='%.1f' % r['seconds'])
    for name in ('W', 'N'):
        t = time.time()
        chunks = [(name, a, min(a + 2048, 65536)) for a in range(0, 65536, 2048)]
        with Pool(en.MAX_THREADS) as pool: res = pool.map(_full_exact_chunk, chunks)
        nz = sum(r[0] for r in res); ns = sum(r[1] for r in res); worst = min(r[2] for r in res)
        ok = nz == TABLE2[name][2] and worst == 0 and (ns, nz - ns) == SPLIT[name]
        report('a.%s full exact scan (65536 co-players, no floats)' % name, ok, min=worst, zero=nz, struct=ns,
               identities=nz - ns, sec='%.1f' % (time.time() - t))
    # float engine: same structural split, no negative margin
    for name in ('W', 'N'):
        S = en.f_scan(np.array([float(x) for x in WIT[name]]))
        m = en.WEDGES[name][0]*(S['W'][:, 2] - S['W'][:, 1]); st = (S['W'][:, 1] == 0) & (S['W'][:, 2] == 0)
        ok = st.sum() == SPLIT[name][0] and (np.abs(m[~st]) <= 1e-12).sum() == SPLIT[name][1] and m.min() > -1e-15
        report('a.%s float engine split' % name, ok, struct=int(st.sum()), ties=int((np.abs(m[~st]) <= 1e-12).sum()),
               min='%.2e' % m.min())
    RES['a.seconds'] = time.time() - t0

# ------------------------------------------------------------------ (b)
def rand_dyadic(rng): return F(rng.randint(1, 2**20 - 1), 2**20)
def rand_point(rng, pfree=None):
    pfree = rng.choice([0.2, 0.4, 0.6, 0.8]) if pfree is None else pfree
    return [rand_dyadic(rng) if rng.random() < pfree else F(rng.randint(0, 1)) for _ in range(16)]

def test_b():
    import ctypes
    L0 = ctypes.CDLL(os.path.join(M2F, '..', 'cannon', 'libm2f.dylib'))
    dp = ctypes.POINTER(ctypes.c_double); L0.c_freqs.argtypes = [dp, dp, dp]
    rng = random.Random(20260927); t0 = time.time()
    bad_int = 0; bad_exp = 0; werr = 0.0; derr = 0.0; bitdiff = 0; mpq_frac = 0; nencl = 0
    for it in range(200):
        p = rand_point(rng); g = rng.randrange(65536); q = ex.g2p(g)
        X = ex.lead_exact(p, q); pf = np.array([float(x) for x in p])
        Fz = en.f_eval(pf, g); Iv = en.i_eval(pf, pf, g)
        if not (Fz['E'] == X['E'] == Iv['E'] and list(Fz['Ek']) == X['Ek'] and Fz['supp'] == X['supp'] == Iv['supp']):
            bad_exp += 1
        werr = max(werr, max(abs(float(a) - b) for a, b in zip(X['w16'], Fz['w16'])))
        derr = max(derr, abs(float(X['D']) - Fz['D'])/float(X['D']))
        for k in range(16):
            nencl += 1
            if not (F(Iv['wlo'][k]) <= X['w16'][k] <= F(Iv['whi'][k])): bad_int += 1
        for o in range(4):
            if not (F(Iv['Wlo'][o]) <= X['W'][o] <= F(Iv['Whi'][o])): bad_int += 1
        if not (F(Iv['Dlo']) <= X['D'] <= F(Iv['Dhi'])): bad_int += 1
        W0 = np.zeros(4); qf = np.array([float(x) for x in q])
        L0.c_freqs(pf.ctypes.data_as(dp), qf.ctypes.data_as(dp), W0.ctypes.data_as(dp))
        if not np.array_equal(W0, Fz['W'].copy()) and not np.array_equal(W0, Fz['w16'].reshape(4, 4).sum(1)): bitdiff += 1
        if it < 40 and ex.lead_exact(p, q, num='fraction')['w16'] != X['w16']: mpq_frac += 1
    report('b.float vs exact, 200 random pairs', bad_exp == 0 and werr < 1e-9 and derr < 1e-9,
           exponent_or_support_mismatch=bad_exp, max_abs_w_err='%.1e' % werr, max_rel_Dstar_err='%.1e' % derr)
    report('b.interval encloses exact', bad_int == 0, violations=bad_int, checked=nencl + 200*5)
    report('b.float W identical to frozen M2f libm2f', bitdiff == 0, differing=bitdiff)
    report('b.mpq == Fraction', mpq_frac == 0, differing=mpq_frac, checked=40)
    # non-dyadic rational points: interval point enclosures contain the exact values
    bad = 0
    for it in range(60):
        p = [F(rng.randint(1, 8), 9) if rng.random() < 0.5 else F(rng.randint(0, 1)) for _ in range(16)]
        g = rng.randrange(65536); X = ex.lead_exact(p, ex.g2p(g)); lo, hi = en.point_box(p)
        Iv = en.i_eval(lo, hi, g)
        bad += sum(not (F(Iv['Wlo'][o]) <= X['W'][o] <= F(Iv['Whi'][o])) for o in range(4))
    report('b.interval encloses exact at non-dyadic points (k/9)', bad == 0, violations=bad, pairs=60)
    RES['b.seconds'] = time.time() - t0

# ------------------------------------------------------------------ (c)
def Nval(p, q, sign):
    """exact N_q = D* sign (w_DC - w_CD), D*, the A*_k and E at a point"""
    r = ex.lead_exact(p, q)
    return r['D']*sign*(r['W'][2] - r['W'][1]), r['D'], [r['D']*w for w in r['w16']], r['E']

def extrap(vals, fr, y, a=F(1, 3), b=F(2, 3)):
    """multilinear interpolation from the grid {a,b}^d (bit t of m selects b for coordinate fr[t])"""
    tot = 0
    for m, v in enumerate(vals):
        if v == 0: continue
        w = F(1)
        for t, j in enumerate(fr):
            w *= (y[j] - a)/(b - a) if (m >> t) & 1 else (b - y[j])/(b - a)
        tot += v*w
    return tot

def test_c():
    rng = random.Random(37); t0 = time.time()
    ncase = 0; nsame = 0; nhigh = 0; nfaces = 0; bad = []; badml = 0; nlow = 0; nface_high = 0; nD0 = 0
    for it in range(260):
        d = rng.randint(2, 5); fr = sorted(rng.sample(range(16), d))
        pat = [F(rng.randint(0, 1)) for _ in range(16)]
        g = rng.randrange(65536); q = ex.g2p(g); sign = rng.choice([1, -1])
        grid = []
        for m in range(1 << d):
            x = pat[:]
            for t, j in enumerate(fr): x[j] = F(2, 3) if (m >> t) & 1 else F(1, 3)
            grid.append(Nval(x, q, sign))
        E = grid[0][3]
        if any(gv[3] != E for gv in grid): bad.append(('E not constant', it)); continue
        GN = [gv[0] for gv in grid]; GD = [gv[1] for gv in grid]; GA = [[gv[2][k] for gv in grid] for k in range(16)]
        # multilinearity: a random interior rational point
        y = pat[:]
        for j in fr: y[j] = F(rng.randint(1, 96), 97)
        Ny, Dy, Ay, Ey = Nval(y, q, sign)
        if extrap(GN, fr, y) != Ny or extrap(GD, fr, y) != Dy or any(extrap(GA[k], fr, y) != Ay[k] for k in range(16)):
            badml += 1
        # corners and random face points
        targets = []
        for m in range(1 << d):
            c = pat[:]
            for t, j in enumerate(fr): c[j] = F((m >> t) & 1)
            targets.append(c)
        for _ in range(3):
            c = pat[:]
            for j in fr: c[j] = rng.choice([F(0), F(1), F(rng.randint(1, 96), 97)])
            targets.append(c); nfaces += 1
        for c in targets:
            Nc, Dc, Ac, Ec = Nval(c, q, sign); ncase += 1
            Npred = extrap(GN, fr, c); Dpred = extrap(GD, fr, c); Apred = [extrap(GA[k], fr, c) for k in range(16)]
            if Ec < E: nlow += 1; bad.append(('E_c<E', it)); continue
            if Ec == E:
                nsame += 1
                if not (Npred == Nc and Dpred == Dc and Apred == Ac): bad.append(('E_c=E mismatch', it))
            else:
                nhigh += 1; nface_high += any(0 < c[j] < 1 for j in fr)
                if not (Npred == 0 and Dpred == 0 and all(a == 0 for a in Apred)): bad.append(('E_c>E nonzero', it))
            if Dpred == 0: nD0 += 1
    report('c.L3 exact extrapolation (corners + faces)', not bad and badml == 0 and nlow == 0,
           strata_q=260, targets=ncase, same_E=nsame, higher_E=nhigh, higher_E_on_open_faces=nface_high,
           multilinearity_failures=badml, failures=len(bad), first=bad[:3])
    # the prompt's version: exact evaluation 1e-9 from a corner, d up to 8
    delta = F(1, 10**9); worst = 0.0; nbad = 0; nn = 0; nh = 0
    for it in range(150):
        d = rng.randint(3, 8); fr = sorted(rng.sample(range(16), d))
        pat = [F(rng.randint(0, 1)) for _ in range(16)]; g = rng.randrange(65536); q = ex.g2p(g); sign = rng.choice([1, -1])
        x = pat[:]
        for j in fr: x[j] = F(1, 2)
        E = ex.lead_exact(x, q)['E']
        for _ in range(4):
            c = pat[:]
            for j in fr: c[j] = F(rng.randint(0, 1))
            y = c[:]
            for j in fr: y[j] = delta if c[j] == 0 else 1 - delta
            Ny, Dy, _, Ey = Nval(y, q, sign); Nc, Dc, _, Ec = Nval(c, q, sign)
            pred = Nc if Ec == E else 0; predD = Dc if Ec == E else 0
            nh += Ec > E; nn += 1
            err = max(abs(float(Ny - pred)), abs(float(Dy - predD)))/float(delta)
            worst = max(worst, err)
            if err > 1e4 or Ey != E: nbad += 1
    report('c.near-corner (1e-9) values match the L3 rule', nbad == 0, pairs=nn, higher_E=nh,
           max_abs_err_over_delta='%.1f' % worst)
    # E and D* against exact finite-eps determinants (matrix-tree theorem, no leading-order arithmetic)
    nb = 0; worst = 0.0
    for it in range(40):
        p = [rng.choice([F(0), F(1), F(rng.randint(1, 7), 8)]) for _ in range(16)]
        g = rng.randrange(65536); q = ex.g2p(g); L = ex.lead_exact(p, q)
        e1, e2 = F(1, 10**6), F(1, 10**7)
        pi1, A1 = ex.finite_eps(p, q, e1); pi2, A2 = ex.finite_eps(p, q, e2)
        r1 = sum(A1)/e1**L['E']; r2 = sum(A2)/e2**L['E']
        # D(eps)/eps^E -> D*, with error O(eps): the two errors must shrink ~10x
        err1 = abs(float(r1 - L['D'])); err2 = abs(float(r2 - L['D']))
        wd = max(abs(float(a - b)) for a, b in zip(pi2, L['w16']))
        rel = err2/float(L['D']); worst = max(worst, rel)
        if not (rel < 1e-4 and (err2 <= err1/5 or err1 == 0) and wd < 1e-4): nb += 1
    report('c.E and D* vs exact finite-eps tree sums', nb == 0, cases=40, max_rel_err_eps1e7='%.1e' % worst)
    RES['c.seconds'] = time.time() - t0

# ------------------------------------------------------------------ (d)
def _gap_chunk(seed):
    rng = random.Random(seed); out = []
    for _ in range(2500):
        a = rng.randrange(65536); b = rng.randrange(65536)
        W = ex.lead_exact(ex.g2p(a), ex.g2p(b), as_fraction=False)['W']
        m = W[2] - W[1]; out.append((a, b, F(int(m.numerator), int(m.denominator))))
    return out

def test_d():
    t0 = time.time()
    with Pool(en.MAX_THREADS) as pool: res = pool.map(_gap_chunk, range(1000, 1040))
    pairs = [x for r in res for x in r]
    nz = [(abs(m), a, b) for a, b, m in pairs if m != 0]
    mn, ma, mb = min(nz); nzero = len(pairs) - len(nz)
    nbelow = sum(1 for m, a, b in nz if m < F(1, 37))
    # the claim of facts_m2.md:114 / ms.tex:326 / M2f methods.tex:76 -- informational, refuted
    RES['d.claim_1_37'] = dict(info=True, holds=(mn >= F(1, 37)), min_nonzero=str(mn), at=(ma, mb), n_below=nbelow)
    print('INFO d.claim "smallest nonzero |w_DC-w_CD| over binary pairs is 1/37": %s -- min nonzero in the sample'
          ' %s at (p,q)=(%d,%d); %d of %d nonzero sampled margins are below 1/37'
          % ('holds' if mn >= F(1, 37) else 'REFUTED', mn, ma, mb, nbelow, len(nz)))
    # the rigorous replacement: integrality, D* (w_DC - w_CD) in Z, via the interval engine alone
    ferr = 0.0; ibad = 0; intbad = 0; maxD = 0; nchk = 0
    for a, b, m in pairs[:20000]:
        pa = en.g2p(a); nchk += 1
        W = en.f_eval(pa, b)['W']; ferr = max(ferr, abs(float(m) - (W[2] - W[1])))
        Iv = en.i_eval(pa, pa, b); lo = Iv['Wlo'][2] - Iv['Whi'][1]; hi = Iv['Whi'][2] - Iv['Wlo'][1]
        if not (lo - 1e-15 <= float(m) <= hi + 1e-15): ibad += 1
        A, Dn = en.binary_pair_exact(a, b); maxD = max(maxD, Dn)
        if F(sum(A[8:12]) - sum(A[4:8]), Dn) != m: intbad += 1
    report('d.integrality: binary_pair_exact (interval + integer rounding) == exact', intbad == 0 and ibad == 0,
           pairs=nchk, mismatches=intbad, interval_violations=ibad, max_Dstar=maxD, max_float_err='%.1e' % ferr)
    report('d.exact margins on 10^5 random binary pairs', True, pairs=len(pairs), exact_zeros=nzero,
           min_nonzero=str(mn), min_nonzero_float='%.2e' % float(mn), sec='%.1f' % (time.time() - t0))
    RES['d.seconds'] = time.time() - t0

# ------------------------------------------------------------------ (e)
def test_e():
    t0 = time.time(); rng = np.random.default_rng(5)
    # ratio route, S (report 1 / certS_ratio.py): 15-dimensional box
    c = [F(n, 32) for n in [32, 27, 15, 25, 29, 28, 9, 6, 21, 14, 27, 15, 6, 12, 30, 22]]
    r = ce.certify_ratio_box(c, F(1, 1024), wedge='S')
    report('e.ratio box S (d=15, h=1/1024)', r['certified'] and r['rho_min'] == F(419698, 307347) and r['rho_argmin'] == 2902
           and r['n_struct_zero'] == 25159 and r['n_unfav_struct_zero'] == 3710,
           U_over_L='%.6f' % float(r['U_over_L']), rho_min=str(r['rho_min']), argmin=r['rho_argmin'],
           struct=r['n_struct_zero'], unfav_struct=r['n_unfav_struct_zero'], exact_needed=r['n_ratio_exact'],
           sec='%.2f' % r['seconds'])
    r2 = ce.certify_ratio_box(c, F(1, 64), wedge='S')
    report('e.ratio box S negative control (h=1/64 too wide)', not r2['certified'], U_over_L='%.3f' % float(r2['U_over_L']),
           n_fail=r2['n_fail'])
    # whole open 5-cube of the binary-anchored family
    r = ce.certify_open_stratum_all_q('1*1*0*0101*1*010', wedge='W')
    report('e.open stratum 1*1*0*0101*1*010 (d=5) certified for all q', r['certified'],
           corner_counts=r['corner_counts'], sec='%.2f' % r['seconds'])
    # soundness sample: random interior points of the 5-cube, including near faces
    pat = en.parse_pattern('1*1*0*0101*1*010'); fr = en.free_coords(pat); worst = 1.0
    for k in range(24):
        x = en.centre(pat); u = rng.random(len(fr)); u[rng.random(len(fr)) < 0.3] = rng.choice([1e-6, 1 - 1e-6])
        x[fr] = u; S = en.f_scan(x, nthreads=3); worst = min(worst, (S['W'][:, 2] - S['W'][:, 1]).min())
    report('e.  sampled check: 24 random points x 65536 q', worst > -1e-12, min_float_margin='%.2e' % worst)
    # negative controls: the W witness stratum is NOT wholly competitive, and is NOT excludable
    r = ce.certify_open_stratum_all_q('1*1*001*1*1**0*0', wedge='W')
    report('e.W witness stratum not certified as a whole (control)', not r['certified'], failing_q=r['failing_q'][:3],
           sec='%.2f' % r['seconds'])
    for name in ('W', 'N', 'S', 'E'):
        pt = en.pattern_of(WIT[name])
        r = ce.exclude_stratum(pt, wedge=name, search_all=True)
        report('e.%s witness stratum is not excluded (it holds a friendly rival)' % name, not r['excluded'],
               candidates_tried=r['n_candidates_tried'], sec='%.1f' % r['seconds'])
    r = ce.exclude_stratum('1*1*0*0101*1*010', wedge='W')
    report('e.5-cube not excluded (control)', not r['excluded'], sec='%.2f' % r['seconds'])
    # boxes by L2
    PW = WIT['W']; patW = en.pattern_of(PW); frW = en.free_coords(patW)
    r = ce.certify_box(patW, [PW[j] - F(1, 16) for j in frW], [PW[j] + F(1, 16) for j in frW], wedge='W')
    report('e.box W witness +-1/16 (d=7)', r['certified'], by_corners=r['n_by_corners'], vertex_q=r['n_vertex_q'],
           exact=r.get('n_vertex_exact'), min_vertex_lower='%.3g' % r.get('min_vertex_margin_lower_bound', float('nan')),
           volume='%.3g' % r['volume'], sec='%.1f' % r['seconds'])
    wbox = ([float(PW[j] - F(1, 16)) for j in frW], [float(PW[j] + F(1, 16)) for j in frW])
    r = ce.certify_box(patW, [PW[j] - F(1, 8) for j in frW], [PW[j] + F(1, 8) for j in frW], wedge='W')
    report('e.box W witness +-1/8 fails (control)', not r['certified'], counterexample_q=r.get('counterexample_q'),
           sec='%.1f' % r['seconds'])
    PN = WIT['N']; patN = en.pattern_of(PN); frN = en.free_coords(patN)
    r = ce.certify_box(patN, [PN[j] - F(1, 12) for j in frN], [PN[j] + F(1, 12) for j in frN], wedge='N')
    report('e.box N witness +-1/12 (d=8)', r['certified'], by_corners=r['n_by_corners'], vertex_q=r['n_vertex_q'],
           exact=r.get('n_vertex_exact'), min_vertex_lower='%.3g' % r.get('min_vertex_margin_lower_bound', float('nan')),
           volume='%.3g' % r['volume'], sec='%.1f' % r['seconds'])
    # soundness sample inside the certified W box: random points, all co-players
    worst = 1.0
    for k in range(24):
        x = np.array([float(v) for v in PW]); x[frW] = rng.uniform(wbox[0], wbox[1])
        S = en.f_scan(x, nthreads=3); worst = min(worst, (S['W'][:, 2] - S['W'][:, 1]).min())
    report('e.  sampled check inside the W box: 24 points x 65536 q', worst > -1e-12, min_float_margin='%.2e' % worst)
    # identities: the exact zeros of the W witness are identities on the whole stratum
    Iw = en.i_scan(*en.point_box(PW), +1)
    ids = np.where(Iw['cls'] == en.CLS_UNRES)[0]; pos = np.where(Iw['cls'] == 1)[0]
    pick = rng.choice(ids, 60, replace=False); ok_id = sum(ce.is_identity(patW, +1, int(q)) for q in pick)
    ok_pos = sum(not ce.is_identity(patW, +1, int(q)) for q in rng.choice(pos, 60, replace=False))
    report('e.is_identity: W-witness ties are identities, positives are not', ok_id == 60 and ok_pos == 60,
           ties_identity='%d/60' % ok_id, positives_not_identity='%d/60' % ok_pos)
    # exclusion on 300 efficient strata with p8 free (slab probe P8) + sampled soundness
    SP = '/private/tmp/claude-501/-Users-martin-Documents-CL1-PRFull/5fece454-5fef-4210-8fd4-4a077a942c23/scratchpad/slab/effstrata8.npy'
    if os.path.exists(SP):
        D = np.load(SP); idx = rng.choice(len(D), 300, replace=False); t = time.time(); nex = 0; sound = 0; nchk = 0; dims = []
        for i in idx:
            pat = ''.join('0*1'[v] for v in D[i]); r = ce.exclude_stratum(pat, wedge='W'); nex += r['excluded']
            dims.append(r['d'])
            if r['excluded'] and nchk < 60:
                pp = en.parse_pattern(pat); frr = en.free_coords(pp); nchk += 1
                X = np.repeat(en.centre(pp)[None, :], 40, 0); U = rng.random((40, len(frr)))
                U[rng.random(U.shape) < 0.2] = 1e-5; X[:, frr] = U
                m = np.array([en.f_eval(x, r['q'])['W'] for x in X]); mm = (m[:, 2] - m[:, 1]).max()
                sound += mm < 0
        report('e.exclude_stratum on 300 random efficient P8 strata', nex == 300 and sound == nchk,
               excluded=nex, d_range='%d-%d' % (min(dims), max(dims)), sampled_soundness='%d/%d' % (sound, nchk),
               sec_per_stratum='%.3f' % ((time.time() - t)/300))
    else:
        skipped('e.exclude_stratum on 300 random efficient P8 strata',
                'its input, effstrata8.npy, a sample of efficient strata with p8 free kept in a session scratchpad of the '
                'author, no longer exists and is not deposited; the exclusion certificates themselves are rerun and '
                'checked by the step "exclusions" of reproduce_certificates.sh')
    RES['e.seconds'] = time.time() - t0

# ------------------------------------------------------------------ timings
def test_t():
    PW = np.array([float(x) for x in WIT['W']]); T = {}
    def tm(name, f, n):
        t = time.time()
        for _ in range(n): f()
        T[name] = (time.time() - t)/n
    tm('f_eval (python call)', lambda: en.f_eval(PW, 1234), 2000)
    tm('i_eval (python call)', lambda: en.i_eval(PW, PW, 1234), 2000)
    tm('f_scan 65536 q, 1 thread', lambda: en.f_scan(PW, nthreads=1), 5)
    tm('f_scan 65536 q, 3 threads', lambda: en.f_scan(PW, nthreads=3), 5)
    tm('i_scan 65536 q, 1 thread', lambda: en.i_scan(PW, PW, 1, nthreads=1), 5)
    tm('i_scan 65536 q, 3 threads', lambda: en.i_scan(PW, PW, 1, nthreads=3), 5)
    tm('corner_scan d=5 x 65536 q, 3 threads', lambda: en.corner_scan('1*1*0*0101*1*010', 1, None, nthreads=3), 2)
    tm('corner_scan d=7 x 65536 q, 3 threads', lambda: en.corner_scan('1*1*001*1*1**0*0', 1, None, nthreads=3), 1)
    pq = [F(x) for x in WIT['N']]
    tm('exact lead_exact (mpq)', lambda: ex.lead_exact(pq, ex.g2p(4321)), 300)
    tm('exact lead_exact (Fraction)', lambda: ex.lead_exact(pq, ex.g2p(4321), num='fraction'), 300)
    for k, v in T.items():
        print('   %-42s %10.3f ms' % (k, v*1e3))
    T['corner eval, us per (q, corner), 3 threads'] = T['corner_scan d=7 x 65536 q, 3 threads']/(128*65536)*1e6
    print('   %-42s %10.3f us' % ('corner_scan per (q,corner), 3 threads', T['corner eval, us per (q, corner), 3 threads']))
    RES['timings'] = T

if __name__ == '__main__':
    t0 = time.time()
    which = sys.argv[1:] or ['a', 'b', 'c', 'd', 'e', 't']
    for w in which:
        print('==== (%s)' % w, flush=True); globals()['test_' + w]()
    nfail = sum(1 for k, v in RES.items() if isinstance(v, dict) and 'ok' in v and not v['ok'])
    for k, v in RES.items():
        if isinstance(v, dict) and v.get('info'): print('INFO', k, v)
    nskip = sum(1 for v in RES.values() if isinstance(v, dict) and v.get('skipped'))
    print('==== %d checks, %d failed, %d skipped, %.0f s' % (sum(1 for v in RES.values() if isinstance(v, dict) and 'ok' in v),
                                                          nfail, nskip, time.time() - t0), flush=True)
