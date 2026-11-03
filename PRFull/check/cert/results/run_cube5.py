"""Referee P2 (d): the binary-anchored family of wedge W.  Certifies the WHOLE open 5-dimensional stratum
1*1*0*0101*1*010 by the corner test (L3 + L4: 32 corners x 65536 co-players), and then decides every
one of the 3^5 = 243 strata of its closed cube:
  * p1 in {*,1}, p5 in {*,0}, p3,p10,p12 in {*,0,1} (108 strata): every point a friendly rival (corner test
    + exact efficiency class);
  * p1 = 0 (81 strata): not efficient (exact class, L1);
  * p5 = 1, p1 in {*,1} (54 strata): excluded -- one co-player beats every point of the open stratum
    (L3/L4 negative direction with a strictly negative centre margin).
(Free coordinates: p1 = p_{CC,CD}, p3 = p_{CC,DD}, p5 = p_{CD,CD}, p10 = p_{DC,DC}, p12 = p_{DD,CC}.)
Run: cd PRFull/check && nice -n 10 python3 cert/results/run_cube5.py"""
import sys, os, itertools; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
PAT = '1*1*0*0101*1*010'
if __name__ == '__main__':
    t0 = time.time(); T = {}; L = []; sign, eff = en.WEDGES['W']
    pat = en.parse_pattern(PAT); fcs = en.free_coords(pat)          # [1,3,5,10,12]
    # 1. the open 5-cube --------------------------------------------------------------
    t = time.time()
    R = ce.certify_open_stratum_all_q(PAT, wedge='W', nthreads=NTH)
    Cs = en.corner_scan(PAT, sign, None, nthreads=NTH); cnt = Cs['cnt']
    ident = (cnt[:, en.C_POS] == 0) & (cnt[:, en.C_NEG] == 0) & (cnt[:, en.C_UNR] == 0) & (cnt[:, en.C_LOW] == 0)
    cen = [F(1, 2) if j in fcs else F(int(pat[j])) for j in range(16)]
    Sc = en.i_scan(*en.point_box(cen), sign, nthreads=NTH)
    n_struct = int((Sc['cls'] == 0).sum())
    spc = ex.selfplay(cen)
    T['open_cube'] = time.time() - t
    # exact spot check: margins exactly 0 at random rational interior points
    t = time.time(); rng = random.Random(5); nz = 0; nt = 0
    for _ in range(6):
        x = list(cen)
        for j in fcs: x[j] = F(rng.randint(1, 96), 97)
        for q in rng.sample(range(NQ), 150):
            nt += 1; nz += ex.margin(x, ex.g2p(q), sign) == 0
    fmax = 0.0; nr = np.random.default_rng(3)
    for _ in range(100):
        x = np.array([float(v) for v in cen]); x[fcs] = nr.random(5)
        W = en.f_scan(x, nthreads=NTH)['W']; fmax = max(fmax, float(np.abs(W[:, 2] - W[:, 1]).max()))
    T['spot_checks'] = time.time() - t
    # 2. every stratum of the closed cube ---------------------------------------------------
    t = time.time()
    fr_ok = []; fr_bad = []; ineff = []; ineff_bad = []; excl = []; excl_bad = []; pool = []
    for vals in itertools.product('*01', repeat=5):
        s = list(PAT)
        for j, v in zip(fcs, vals): s[j] = v
        s = ''.join(s); p1, p3, p5, p10, p12 = vals
        if p1 == '0':
            c = en.efficiency_class(s); (ineff if not en.meets(c, 'C') else ineff_bad).append((s, c)); continue
        if p5 == '1':
            r = ce.exclude_stratum(s, wedge='W', candidates=pool, nthreads=NTH)
            if r['excluded']:
                excl.append((s, r['q'], str(r['centre_margin_exact']) if r['centre_margin_exact'] is not None
                             else r['centre_margin_enclosure']))
                if r['q'] not in pool: pool.append(r['q'])
            else: excl_bad.append(s)
            continue
        r = ce.certify_open_stratum_all_q(s, wedge='W', nthreads=NTH)
        (fr_ok if r['certified'] else fr_bad).append((s, r['d'], r['efficiency_class']))
    T['closed_cube_strata'] = time.time() - t
    out = dict(name='cube5', stratum=PAT, wedge='W', open_stratum=dict(
        certified=R['certified'], efficiency_class=R['efficiency_class'], selfplay_centre_exact=[S(x) for x in spc],
        n_q_checked=R['n_q_checked'], n_q_failing=R['n_q_failing'], corner_counts=R['corner_counts'],
        n_q_identically_zero=int(ident.sum()), n_struct_zero=n_struct, n_nonstructural_identities=int((ident & (Sc['cls'] != 0)).sum()),
        exact_spot_zero=f"{nz}/{nt}", float_max_abs_margin=fmax),
        faces=dict(n_friendly_rival_strata=len(fr_ok), friendly_rival_strata=fr_ok, not_certified=fr_bad,
                   n_inefficient_strata=len(ineff), inefficient=ineff, efficient_on_p1_0=ineff_bad,
                   n_excluded_strata=len(excl), excluded=excl, not_excluded=excl_bad, exclusion_pool=pool),
        timing_seconds=T, seconds_total=time.time() - t0)
    out['certified'] = bool(R['certified'] and len(fr_ok) == 108 and not fr_bad)
    out['closed_cube_decided'] = bool(out['certified'] and len(ineff) == 81 and len(excl) == 54)
    with open(os.path.join(HERE, 'cube5_stratum.json'), 'w') as f: json.dump(jsonable(out), f, indent=1)
    o = out['open_stratum']; fz = out['faces']
    L += [f"# cube5: open stratum {'CERTIFIED' if R['certified'] else 'NOT CERTIFIED'}; closed cube "
          f"{'fully decided' if out['closed_cube_decided'] else 'NOT fully decided'}",
          "STATEMENT.  (i) Every point of the open 5-dimensional stratum 1*1*0*0101*1*010 (p1,p3,p5,p10,p12 in (0,1),",
          "the other eleven coordinates pinned as shown) is a friendly rival of wedge W in the eps -> 0 limit: exact",
          "self-play class C, and against every binary memory-two co-player (hence every memory-two co-player) the",
          "margin w_DC - w_CD is IDENTICALLY ZERO on the stratum -- N_q = 0 at all 32 corners for all 65536 q.",
          "(ii) The friendly rivals in the closed cube are exactly the points with p1 > 0 and p5 < 1:",
          "all 108 strata with p1 in {*,1}, p5 in {*,0} are certified (every point a friendly rival), the 81 strata",
          "with p1 = 0 are not efficient, and each of the 54 strata with p5 = 1 (p1 > 0) is excluded by one co-player",
          "that strictly beats every point of it.",
          "",
          f"open stratum: efficiency class {o['efficiency_class']} (exact, L1); exact self-play at the centre ({', '.join(o['selfplay_centre_exact'])})",
          f"corner test: {o['n_q_checked']} co-players checked, {o['n_q_failing']} failing; corner counts {o['corner_counts']}",
          f"N_q identically 0 for {o['n_q_identically_zero']} of 65536 co-players = {o['n_struct_zero']} structural zeros "
          f"(w_CD = w_DC = 0) + {o['n_nonstructural_identities']} non-structural identities (w_CD = w_DC > 0)",
          f"exact spot check: margin exactly 0 at {o['exact_spot_zero']} (random rational interior point, random q) pairs;",
          f"float sanity: max |w_DC - w_CD| over 100 random interior points x 65536 = {o['float_max_abs_margin']:.3g}",
          "",
          f"closed cube: {fz['n_friendly_rival_strata']}/108 friendly-rival strata certified (not certified: {fz['not_certified']});",
          f"   {fz['n_inefficient_strata']}/81 strata with p1 = 0 inefficient (exceptions: {fz['efficient_on_p1_0']});",
          f"   {fz['n_excluded_strata']}/54 strata with p5 = 1 excluded (not excluded: {fz['not_excluded']}); co-players used {fz['exclusion_pool']}",
          "   certified strata (pattern, dimension, class):"]
    for s, d, c in fr_ok: L.append(f"      {s}  {d}  {c}")
    L.append("   excluded strata (pattern, co-player, centre margin):")
    for s, q, m in excl: L.append(f"      {s}  q={q}  {m}")
    L.append("   inefficient strata (pattern, class):  " + '  '.join(f"{s}:{c}" for s, c in ineff))
    L.append("")
    L.append("TIMING (s, <= 3 threads, nice 10): " + ', '.join(f"{k} {v:.1f}" for k, v in T.items()) + f"; total {out['seconds_total']:.1f}")
    txt = '\n'.join(L) + '\n'
    open(os.path.join(HERE, 'cube5_stratum.txt'), 'w').write(txt); print(txt)
