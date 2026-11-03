"""Driver for the positive-volume certificates of referee point P2 (friendly rivals, memory two,
eps -> 0 before any inequality).  Uses only PRFull/check/cert (engine, exact, certify).

box_certificate(...) proves that a closed box B = prod [lo_s, hi_s] over the free coordinates of a
stratum (strictly inside the open stratum) consists of friendly rivals of a wedge, by the VERTEX
ROUTE (lemma L2), independently of the corner lemma L3:

  * efficiency: exact (lemma L1): the self-play support is decided by eps-exponents, which depend
    only on the pattern, so the class of the open stratum is that of any one of its points; the
    exact rational self-play is also computed at the centre and at every vertex;
  * competitiveness against each of the 65536 binary co-players q: N_q = D* x margin is
    multilinear in the free coordinates with D* > 0 on the open stratum, so min_B N_q is attained
    at a vertex.  Every one of the 2^d x 65536 vertex margins is (i) a strictly positive
    outward-rounded interval enclosure, or (ii) a structural zero (w_CD = w_DC = 0 by exponents),
    or (iii) an exact rational (gmpy2.mpq) when the enclosure contains 0.  exact_all=False
    replaces (iii) by the L3/L4 corner identity test (N_q = 0 identically) for those q whose
    every corner value is 0 -- used only where stated.

A certified-negative value anywhere is a counterexample and the certificate fails.
Also: the L3/L4 corner table for every q (independent second route), an exact multilinearity
check of N_q on the box for a sample of q, and a float scan at random interior points (sanity
only, not part of the proof).  Local rule: at most 3 threads/processes, run under nice.
"""
import sys, os, time, json, math, random
HERE = os.path.dirname(os.path.abspath(__file__))
CHECK = os.path.dirname(os.path.dirname(HERE))
if CHECK not in sys.path: sys.path.insert(0, CHECK)
import numpy as np
from fractions import Fraction as F
from multiprocessing import Pool
from cert import engine as en, exact as ex, certify as ce

NQ = en.NQ
NPROC = min(3, en.MAX_THREADS)
NTH = min(3, en.MAX_THREADS)


# ------------------------------------------------------------------ small helpers
def S(x): return str(F(x))
def mirror_point(p):
    """C <-> D exchange: pbar_s = 1 - p_{15-s} (commutes with the clamp eps + (1-2eps)p)"""
    return [1 - F(p[15 - s]) for s in range(16)]
def mirror_genome(g):
    return sum((1 - ((int(g) >> (15 - s)) & 1)) << s for s in range(16))
def mirror_pattern(pat):
    pat = en.parse_pattern(pat)
    return en.pattern_str([{0: 1, 1: 0, 2: 2}[int(pat[15 - s])] for s in range(16)])
def eff_ok(W, eff):
    """exact efficiency test of a self-play limit W=(w_CC,w_CD,w_DC,w_DD)"""
    if eff == 'C': return W[0] == 1
    if eff == 'A': return W[1] == W[2] == F(1, 2)
    if eff == 'Z': return W[3] == 0
    raise ValueError(eff)
def jsonable(o):
    if isinstance(o, F): return str(o)
    if isinstance(o, dict): return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [jsonable(v) for v in o]
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, np.ndarray): return jsonable(o.tolist())
    return o
def ptstr(p): return '(' + ', '.join(S(x) for x in p) + ')'


# ------------------------------------------------------------------ exact workers
def _exact_chunk(args):
    m, pt, qs, sign = args
    P = [F(x) for x in pt]
    pos = []; neg = []; nz = 0; mn = None; am = None
    for q in qs:
        v = ex.margin(P, ex.g2p(int(q)), sign)
        if v == 0: nz += 1
        elif v < 0: neg.append((int(q), str(v)))
        else:
            pos.append(int(q))
            if mn is None or v < mn: mn, am = v, int(q)
    return m, len(qs), nz, pos, neg, (None if mn is None else str(mn)), am

def _lemma_job(args):
    """exact multilinearity check of N_q = D* sign (w_DC - w_CD) on a box, for one co-player"""
    q, sign, fcs, lo, hi, C, xs = args
    d = len(fcs); nv = 1 << d
    def pt(vals):
        x = [F(v) for v in C]
        for j, v in zip(fcs, vals): x[j] = F(v)
        return x
    def N(x):
        r = ex.lead_exact(x, ex.g2p(q)); return r['D']*sign*(r['W'][2] - r['W'][1]), r['E']
    Nv = []; Es = set()
    for m in range(nv):
        v, e = N(pt([hi[t] if (m >> t) & 1 else lo[t] for t in range(d)])); Nv.append(v); Es.add(e)
    res = []
    for xv in xs:
        lam = [(F(x) - F(a))/(F(b) - F(a)) for x, a, b in zip(xv, lo, hi)]
        interp = F(0)
        for m in range(nv):
            w = F(1)
            for t in range(d): w *= lam[t] if (m >> t) & 1 else 1 - lam[t]
            interp += w*Nv[m]
        val, e = N(pt(xv)); Es.add(e)
        res.append((interp == val, val >= 0, str(val)))
    return q, min(Nv), len(Es) == 1, res


# ------------------------------------------------------------------ the box certificate
def box_certificate(name, wedge, pat, centre, lo, hi, exact_all=True, n_lemma=24, n_float=200,
                    seed=20260927, log=print, max_exact_seconds=780):
    t00 = time.time(); T = {}
    sign, eff = en.WEDGES[wedge]
    pat = en.parse_pattern(pat); fcs = en.free_coords(pat); d = len(fcs); nv = 1 << d
    C = [F(x) for x in centre]
    if not (en.pattern_of(C) == pat).all(): raise ValueError("centre not in the stratum")
    lo = [F(x) for x in lo]; hi = [F(x) for x in hi]
    for j, a, b in zip(fcs, lo, hi):
        if not (0 < a <= C[j] <= b < 1): raise ValueError("box must contain the centre and lie in the open stratum")
    vol = F(1)
    for a, b in zip(lo, hi): vol *= (b - a)
    def vpt(m):
        x = list(C)
        for t, j in enumerate(fcs): x[j] = hi[t] if (m >> t) & 1 else lo[t]
        return x
    V = [vpt(m) for m in range(nv)]
    out = dict(name=name, wedge=wedge, sign=sign, efficiency_required=eff, stratum=en.pattern_str(pat),
               d=d, free_coordinates=fcs, centre=[S(x) for x in C],
               box={f"p{j}": [S(a), S(b)] for j, a, b in zip(fcs, lo, hi)}, volume=str(vol),
               volume_float=float(vol), n_vertices=nv)
    log(f"[{name}] stratum {out['stratum']} d={d}, {nv} vertices, volume {float(vol):.3g}")

    # 1. efficiency (exact, L1) ------------------------------------------------------------
    t = time.time()
    ecls = en.efficiency_class(pat)
    spc = ex.selfplay(C)
    spv_ok = all(eff_ok(ex.selfplay(v), eff) for v in V)
    out['efficiency'] = dict(class_of_open_stratum=ecls, required=eff, meets=en.meets(ecls, eff),
                             selfplay_centre_exact=[S(x) for x in spc], selfplay_centre_efficient=eff_ok(spc, eff),
                             selfplay_all_vertices_exact_efficient=spv_ok)
    T['efficiency'] = time.time() - t
    log(f"[{name}] efficiency class {ecls}, self-play at centre {[S(x) for x in spc]}, all vertices efficient: {spv_ok}")

    # 2. the centre itself (exact point certificate; the SI Table 2 entry) --------------------
    t = time.time()
    P = ce.certify_point(C, wedge=wedge, nthreads=NTH)
    out['centre_point'] = dict(certified=P['certified'], n_pos=P['n_pos'], n_struct_zero=P['n_struct_zero'],
                               n_exact_zero_nonstructural=P['n_exact_zero'], n_neg=P['n_neg'],
                               min_margin=S(P['min_margin']), n_attaining_min=P['n_zero'] if P['min_margin'] == 0 else None,
                               n_unresolved_by_interval=P['n_unresolved_by_interval'])
    T['centre_point'] = time.time() - t
    log(f"[{name}] centre: pos {P['n_pos']}, struct 0 {P['n_struct_zero']}, exact 0 {P['n_exact_zero']}, neg {P['n_neg']} ({T['centre_point']:.1f}s)")

    # 3. corner table (L3/L4) for every co-player: second, independent route ------------------
    t = time.time()
    Cs = en.corner_scan(pat, sign, None, nthreads=NTH)
    cnt = Cs['cnt']
    if (cnt[:, en.C_LOW] > 0).any(): raise RuntimeError("E_c < E: lemma violated")
    ident = (cnt[:, en.C_POS] == 0) & (cnt[:, en.C_NEG] == 0) & (cnt[:, en.C_UNR] == 0)
    T['corner_scan'] = time.time() - t

    # 4. interval scan at every vertex -----------------------------------------------------------
    t = time.time()
    # streamed: cls (nv x NQ int8) is kept; of the enclosures only the per-q minimum lower end over
    # the certified-positive vertices (the exact minimum is recovered in step 6)
    cls = np.empty((nv, NQ), dtype=np.int8); minpos = np.full(NQ, np.inf)
    E0 = S0 = None; Econst = True; Sconst = True
    for m, v in enumerate(V):
        plo, phi = en.point_box(v)
        R = en.i_scan(plo, phi, sign, nthreads=NTH)
        cls[m] = R['cls']
        if E0 is None: E0, S0 = R['E'].copy(), R['supp'].copy()
        else:
            Econst &= bool((R['E'] == E0).all()); Sconst &= bool((R['supp'] == S0).all())
        pm = R['cls'] == 1
        minpos[pm] = np.minimum(minpos[pm], R['mlo'][pm])
    Econst &= bool((Cs['E'] == E0).all())
    T['vertex_interval_scan'] = time.time() - t
    struct_any = (cls == 0).any(0); struct_all = (cls == 0).all(0)
    if not (struct_any == struct_all).all(): raise RuntimeError("structural zeros not constant over the box")
    struct = struct_all
    negpairs = np.argwhere(cls == -1)
    unres = np.argwhere(cls == 2)             # rows (m, q)
    allpos = (cls == 1).all(0)
    has_unres = (cls == 2).any(0)
    log(f"[{name}] vertex interval scan {T['vertex_interval_scan']:.1f}s: struct {int(struct.sum())}, "
        f"all-positive {int(allpos.sum())}, q with unresolved vertices {int(has_unres.sum())}, "
        f"unresolved pairs {len(unres)}, certified-negative pairs {len(negpairs)}")
    out['vertex_route'] = dict(E_constant_over_box_and_corners=Econst, support_constant_over_box=Sconst,
                               n_struct_zero=int(struct.sum()), n_positive_by_interval=int(allpos.sum()),
                               n_q_with_unresolved_vertex=int(has_unres.sum()),
                               n_unresolved_pairs=int(len(unres)), n_certified_negative_pairs=int(len(negpairs)))
    if len(negpairs):
        m, q = negpairs[0]
        out.update(certified=False, counterexample=dict(q=int(q), vertex=[S(x) for x in V[m]]))
        out['timing_seconds'] = T; out['seconds_total'] = time.time() - t00
        return out

    # 5. exact resolution of the unresolved (vertex, q) pairs ------------------------------------
    t = time.time()
    if exact_all: todo = unres
    else: todo = unres[~ident[unres[:, 1]]]
    skipped_ident_q = np.unique(unres[ident[unres[:, 1]], 1]) if not exact_all else np.array([], dtype=int)
    est = len(todo)*0.33e-3/NPROC
    log(f"[{name}] exact pairs {len(todo)} (est. {est:.0f}s on {NPROC} procs)")
    if est > max_exact_seconds:
        out.update(certified=False, refused=True, reason=f"exact stage estimated {est:.0f}s")
        return out
    jobs = []
    order = np.lexsort((todo[:, 1], todo[:, 0])) if len(todo) else []
    todo = todo[order] if len(todo) else todo
    CH = 1500
    for m in range(nv):
        qs = todo[todo[:, 0] == m, 1] if len(todo) else []
        for a in range(0, len(qs), CH):
            jobs.append((m, [S(x) for x in V[m]], [int(x) for x in qs[a:a + CH]], sign))
    n_exact = 0; n_exact_zero = 0; posq = set(); negs = []; ex_min = None; ex_arg = None
    if jobs:
        with Pool(NPROC) as pool:
            for (m, n, nz, pos, neg, mn, am) in pool.imap_unordered(_exact_chunk, jobs):
                n_exact += n; n_exact_zero += nz; posq.update(pos); negs += [(m, q, v) for q, v in neg]
                if mn is not None and (ex_min is None or F(mn) < ex_min): ex_min, ex_arg = F(mn), (am, m)
    T['exact_resolution'] = time.time() - t
    q_exact = np.unique(todo[:, 1]) if len(todo) else np.array([], dtype=int)
    exact_rate = T['exact_resolution']/max(1, n_exact)*NPROC
    log(f"[{name}] exact: {n_exact} pairs in {T['exact_resolution']:.1f}s ({exact_rate*1e3:.3f} ms/pair/proc): "
        f"zero {n_exact_zero}, positive {n_exact - n_exact_zero - len(negs)}, negative {len(negs)}")
    out['exact_resolution'] = dict(n_pairs=int(n_exact), n_pairs_zero=int(n_exact_zero),
                                   n_pairs_positive=int(n_exact - n_exact_zero - len(negs)), n_pairs_negative=len(negs),
                                   n_q_resolved_exactly=int(len(q_exact)),
                                   n_q_resolved_exactly_all_zero=int(len(set(map(int, q_exact)) - posq)),
                                   n_q_resolved_exactly_some_positive=int(len(posq)),
                                   n_q_identity_by_corner_test_not_evaluated=int(len(skipped_ident_q)),
                                   ms_per_pair_per_process=exact_rate*1e3,
                                   min_positive_exact_value=(None if ex_min is None else str(ex_min)),
                                   argmin=(None if ex_arg is None else dict(q=ex_arg[0], vertex=[S(x) for x in V[ex_arg[1]]])))
    if negs:
        m, q, v = negs[0]
        out.update(certified=False, counterexample=dict(q=q, vertex=[S(x) for x in V[m]], exact_margin=v))
        out['timing_seconds'] = T; out['seconds_total'] = time.time() - t00
        return out
    # consistency of the two routes on the zero co-players
    if exact_all:
        allzero_q = set(map(int, q_exact)) - posq
        ident_ns = set(map(int, np.where(ident & ~struct)[0]))
        out['exact_resolution']['all_zero_q_equal_corner_identities'] = (allzero_q == ident_ns)
        out['exact_resolution']['n_corner_identities_nonstructural'] = len(ident_ns)

    # 6. the smallest positive vertex margin, exactly -------------------------------------------
    t = time.time()
    # m0 = smallest lower end over all certified-positive (vertex, q) pairs; every pair whose lower end
    # is <= m0 + 1e-9 belongs to a q with minpos <= m0 + 1e-9: re-scan those q at every vertex
    m0 = float(minpos.min())
    candq = np.where(minpos <= m0 + 1e-9)[0]
    cand = []
    for m, v in enumerate(V):
        R = en.i_scan(*en.point_box(v), sign, candq, nthreads=NTH)
        for i in np.where((R['cls'] == 1) & (R['mlo'] <= m0 + 1e-9))[0]: cand.append((m, int(candq[i])))
    best = None
    for m, q in cand:
        v = ex.margin(V[m], ex.g2p(int(q)), sign)
        if best is None or v < best[0]: best = (v, int(q), int(m))
    if ex_min is not None and ex_min < best[0]: best = (ex_min, ex_arg[0], ex_arg[1])
    # rigour of the argmin: every pair outside cand has lower end > m0 + 1e-9 > (its exact value)
    assert float(best[0]) < m0 + 1e-9
    T['min_margin'] = time.time() - t
    out['min_vertex_margin'] = dict(
        min_over_all_pairs='0' if (struct.any() or n_exact_zero or ident.any()) else str(best[0]),
        attained_by='structural zeros and exact identities (margin identically 0 on the stratum)',
        min_positive_exact=str(best[0]), min_positive_float=float(best[0]), argmin_q=best[1],
        argmin_vertex=[S(x) for x in V[best[2]]], n_candidates_evaluated_exactly=int(len(cand)),
        interval_lower_bound=m0)
    log(f"[{name}] smallest positive vertex margin {best[0]} = {float(best[0]):.6g} at q={best[1]}")

    # 7. corner route summary -----------------------------------------------------------------
    out['corner_route'] = dict(n_q_pos_ok=int(Cs['pos_ok'].sum()), n_q_identity=int(ident.sum()),
                               n_q_identity_nonstructural=int((ident & ~struct).sum()),
                               corner_counts_total={k: int(v) for k, v in zip(
                                   ['pos', 'neg', 'zero_struct', 'zero_integral', 'E_c>E', 'E_c<E', 'unresolved'],
                                   cnt.sum(0))}, seconds=T['corner_scan'])

    # 8. exact multilinearity check of N_q on the box (lemma spot check) -------------------------
    t = time.time()
    rng = random.Random(seed)
    posl = [int(q) for q in np.where(allpos)[0]]; idl = [int(q) for q in np.where(ident & ~struct)[0]]
    qs = [best[1]] + rng.sample(posl, min(len(posl), n_lemma//2)) + rng.sample(idl, min(len(idl), n_lemma - n_lemma//2 - 1))
    jobs = []
    for q in qs:
        xs = [[a + (b - a)*F(rng.randint(1, 96), 97) for a, b in zip(lo, hi)] for _ in range(2)]
        jobs.append((q, sign, fcs, lo, hi, C, xs))
    with Pool(NPROC) as pool: LR = pool.map(_lemma_job, jobs)
    n_match = sum(r[0] for _, _, _, res in LR for r in res); n_tot = sum(len(res) for _, _, _, res in LR)
    out['lemma_check'] = dict(n_q=len(qs), n_interior_points=n_tot, n_multilinear_exact_match=n_match,
                              E_constant=all(r[2] for r in LR),
                              interior_nonneg=all(r[1] for _, _, _, res in LR for r in res),
                              seconds=time.time() - t)
    T['lemma_check'] = time.time() - t
    log(f"[{name}] multilinearity: {n_match}/{n_tot} exact matches over {len(qs)} co-players")

    # 9. float sanity scan at random interior points (not part of the proof) ---------------------
    t = time.time(); nr = np.random.default_rng(seed); fmin = np.inf
    for _ in range(n_float):
        x = np.array([float(v) for v in C])
        for j, a, b in zip(fcs, lo, hi): x[j] = float(a) + (float(b) - float(a))*nr.random()
        R = en.f_scan(x, nthreads=NTH); W = R['W']
        fmin = min(fmin, float((sign*(W[:, 2] - W[:, 1])).min()))
    out['float_sanity'] = dict(n_points=n_float, min_margin_float=fmin, seconds=time.time() - t)
    T['float_sanity'] = time.time() - t

    n_int = int(allpos.sum()); n_ex = int(len(q_exact)); n_st = int(struct.sum()); n_id = int(len(skipped_ident_q))
    out['summary_counts'] = dict(struct_zero=n_st, positive_by_interval_all_vertices=n_int,
                                 resolved_exactly=n_ex, identity_by_corner_test=n_id,
                                 total=n_st + n_int + n_ex + n_id)
    ok = (out['efficiency']['meets'] and spv_ok and eff_ok(spc, eff) and n_st + n_int + n_ex + n_id == NQ
          and Econst and Sconst)
    out['certified'] = bool(ok)
    out['timing_seconds'] = T; out['seconds_total'] = time.time() - t00
    return out


# ------------------------------------------------------------------ write-up
def write_results(out, path_base, statement, extra_lines=()):
    with open(path_base + '.json', 'w') as f: json.dump(jsonable(out), f, indent=1)
    L = []
    L.append(f"# {out['name']}: {'CERTIFIED' if out.get('certified') else 'NOT CERTIFIED'}")
    L.append(statement)
    L.append("")
    L.append(f"wedge {out['wedge']} (margin = {'+' if out['sign'] > 0 else '-'}(w_DC - w_CD)), required efficiency class {out['efficiency_required']}")
    L.append(f"stratum {out['stratum']}  (dimension {out['d']}; free coordinates p_s, s = {out['free_coordinates']})")
    L.append(f"centre  {ptstr(out['centre'])}")
    L.append("box (closed, rational bounds, pinned coordinates as in the stratum):")
    for k, (a, b) in out['box'].items(): L.append(f"   {k:>4} in [{a}, {b}]")
    L.append(f"volume {out['volume']} = {out['volume_float']:.4g};  {out['n_vertices']} vertices")
    e = out['efficiency']
    L.append("")
    L.append("EFFICIENCY (exact).  L1: the eps-exponent of every spanning tree of the self-play chain depends only")
    L.append("on the pattern, and the leading coefficients are positive on the open stratum, so the limit support")
    L.append(f"is the same at every point of the open stratum: class {e['class_of_open_stratum']} (required {e['required']}).")
    L.append(f"   exact self-play at the centre (w_CC,w_CD,w_DC,w_DD) = ({', '.join(e['selfplay_centre_exact'])})")
    L.append(f"   exact self-play efficient at all {out['n_vertices']} vertices: {e['selfplay_all_vertices_exact_efficient']}")
    if 'vertex_route' in out:
        v = out['vertex_route']; sc = out.get('summary_counts', {})
        L.append("")
        L.append("COMPETITIVENESS, vertex route (L2), all 65536 binary co-players x all vertices:")
        L.append(f"   structurally zero (w_CD = w_DC = 0 by exponents, on the whole open stratum): {v['n_struct_zero']}")
        L.append(f"   positive by interval at every vertex (outward-rounded enclosure > 0):        {v['n_positive_by_interval']}")
        if 'exact_resolution' in out:
            x = out['exact_resolution']
            L.append(f"   resolved exactly (some vertex enclosure contains 0; those vertices in mpq):   {x['n_q_resolved_exactly']}")
            L.append(f"        of which exactly 0 at every such vertex: {x['n_q_resolved_exactly_all_zero']}; positive somewhere: {x['n_q_resolved_exactly_some_positive']}")
            if x['n_q_identity_by_corner_test_not_evaluated']:
                L.append(f"   identities by the L3/L4 corner test (N_q = 0 identically; not evaluated): {x['n_q_identity_by_corner_test_not_evaluated']}")
            L.append(f"   exact (vertex, co-player) evaluations: {x['n_pairs']} (zero {x['n_pairs_zero']}, positive {x['n_pairs_positive']}, negative {x['n_pairs_negative']})")
            if 'all_zero_q_equal_corner_identities' in x:
                L.append(f"   the exactly-zero co-players coincide with the corner-test identities: {x['all_zero_q_equal_corner_identities']}")
        if sc: L.append(f"   total {sc['total']} (must be 65536)")
        L.append(f"   certified-negative vertex margins: {v['n_certified_negative_pairs']}")
        L.append(f"   tree exponent E(q) constant over centre, vertices and corners: {v['E_constant_over_box_and_corners']}; support constant: {v['support_constant_over_box']}")
    if 'min_vertex_margin' in out:
        mm = out['min_vertex_margin']
        L.append("")
        L.append(f"MINIMUM EXACT VERTEX MARGIN: {mm['min_over_all_pairs']} (attained by the {mm['attained_by']});")
        L.append(f"   smallest strictly positive vertex margin = {mm['min_positive_exact']} = {mm['min_positive_float']:.6g}")
        L.append(f"   at co-player {mm['argmin_q']} and vertex {ptstr(mm['argmin_vertex'])}")
    if 'corner_route' in out:
        c = out['corner_route']
        L.append("")
        L.append(f"CORNER ROUTE (L3/L4, independent): N_q >= 0 on the whole closed stratum cube for {c['n_q_pos_ok']} co-players;")
        L.append(f"   N_q identically 0 for {c['n_q_identity']} ({c['n_q_identity_nonstructural']} non-structural identities).")
    if 'centre_point' in out:
        p = out['centre_point']
        L.append("")
        L.append(f"CENTRE (exact point certificate; SI Table 2 entry): min margin {p['min_margin']} attained by {p['n_attaining_min']} co-players")
        L.append(f"   = {p['n_struct_zero']} structural + {p['n_exact_zero_nonstructural']} exact non-structural zeros; positive {p['n_pos']}; negative {p['n_neg']}")
    if 'lemma_check' in out:
        l = out['lemma_check']
        L.append("")
        L.append(f"LEMMA CHECK (exact): N_q = D* x margin equals its multilinear interpolation from the vertices at")
        L.append(f"   {l['n_multilinear_exact_match']}/{l['n_interior_points']} random rational interior points ({l['n_q']} co-players); E constant: {l['E_constant']}; interior margins >= 0: {l['interior_nonneg']}")
    if 'float_sanity' in out:
        s = out['float_sanity']
        L.append(f"FLOAT SANITY (not part of the proof): min float margin over {s['n_points']} random interior points x 65536 = {s['min_margin_float']:.3g}")
    for x in extra_lines: L.append(x)
    L.append("")
    L.append("TIMING (s, Martin's Mac, <= 3 threads/processes, nice 10): " +
             ', '.join(f"{k} {v:.1f}" for k, v in out.get('timing_seconds', {}).items()) +
             f"; total {out.get('seconds_total', 0):.1f}")
    txt = '\n'.join(L) + '\n'
    with open(path_base + '.txt', 'w') as f: f.write(txt)
    return txt


# ------------------------------------------------------------------ hybrid: ratio route + vertex route
def hybrid_box_certificate(name, wedge, c, lo, hi, log=print, max_vertex_seconds=780, us_per_eval=2.4):
    """PROVES (if out['certified']) that the closed box prod [lo_s, hi_s] (free coordinates of c's
    stratum, strictly inside it) consists of friendly rivals of the wedge.  Per co-player q:
      (a) w_unf = 0 structurally (w_fav > 0 or both 0): margin >= 0 on the whole open stratum;
      (b) ratio route: rho(c) = w_fav(c)/w_unf(c) >= U/L, interval lower bound (see certify_ratio_box);
      (c) otherwise the vertex route (L2): all 2^d vertex margins >= 0, by a positive interval
          enclosure or an exact rational value.
    Used where d is too large for (c) alone (d = 15 in S) but most co-players pass (b)."""
    t00 = time.time(); T = {}
    sign, eff = en.WEDGES[wedge]
    C = [F(x) for x in c]; pat = en.pattern_of(C); fcs = en.free_coords(pat); d = len(fcs); nv = 1 << d
    lo = [F(x) for x in lo]; hi = [F(x) for x in hi]
    Lr = F(1); Ur = F(1); vol = F(1)
    for j, a, b in zip(fcs, lo, hi):
        if not (0 < a <= C[j] <= b < 1): raise ValueError("box must contain c and lie inside the open stratum")
        Lr *= min(a/C[j], (1 - b)/(1 - C[j])); Ur *= max(b/C[j], (1 - a)/(1 - C[j])); vol *= b - a
    UL = Ur/Lr
    ecls = en.efficiency_class(pat); spc = ex.selfplay(C)
    fav, unf = (2, 1) if sign > 0 else (1, 2)
    t = time.time()
    SW = en.i_scanW(*en.point_box(C), nthreads=NTH)
    unf_on = en.outcome_present(SW['supp'], unf); fav_on = en.outcome_present(SW['supp'], fav)
    n_both0 = int((~unf_on & ~fav_on).sum()); n_unf0 = int((~unf_on & fav_on).sum())
    live = np.where(unf_on)[0]
    rho_lo = np.nextafter(SW['Wlo'][live, fav]/SW['Whi'][live, unf], -np.inf)
    ULf = float(UL); ULhi = float(np.nextafter(ULf, np.inf)) if F(ULf) < UL else ULf
    passr = live[rho_lo >= ULhi]; rest = live[rho_lo < ULhi]
    T['ratio_screen'] = time.time() - t
    est = len(rest)*nv*us_per_eval*1e-6/NTH
    log(f"[{name}] d={d} U/L={float(UL):.4f}: both-0 {n_both0}, unf-0 {n_unf0}, ratio pass {len(passr)}, vertex route {len(rest)} (est {est:.0f}s)")
    out = dict(name=name, wedge=wedge, stratum=en.pattern_str(pat), d=d, centre=[S(x) for x in C],
               box={f"p{j}": [S(a), S(b)] for j, a, b in zip(fcs, lo, hi)}, volume=str(vol), volume_float=float(vol),
               U_over_L=str(UL), efficiency=dict(class_of_open_stratum=ecls, selfplay_centre_exact=[S(x) for x in spc]),
               n_struct_zero=n_both0, n_unfav_struct_zero=n_unf0, n_ratio_interval=int(len(passr)),
               n_vertex_route=int(len(rest)), estimated_vertex_seconds=est)
    if est > max_vertex_seconds:
        out.update(certified=False, refused=True); return out
    t = time.time()
    vlo = np.empty((nv, 16)); vhi = np.empty((nv, 16))
    base_lo, base_hi = en.point_box(C)
    enc_lo = [en.enclose(a) for a in lo]; enc_hi = [en.enclose(b) for b in hi]
    for m in range(nv):
        vlo[m] = base_lo; vhi[m] = base_hi
        for t_, j in enumerate(fcs):
            e = enc_hi[t_] if (m >> t_) & 1 else enc_lo[t_]
            vlo[m, j], vhi[m, j] = e
    n_neg = 0; n_unres = 0; n_exact_zero = 0; minlo = np.inf; negq = []
    unres_all = []
    B = 256
    for a in range(0, len(rest), B):
        ql = rest[a:a + B]
        VS = en.vertex_scan(vlo, vhi, sign, ql, nthreads=NTH, cap=2_000_000)
        n_neg += int((VS['cnt'][:, 1] > 0).sum())
        if n_neg:
            negq = [int(ql[i]) for i in np.where(VS['cnt'][:, 1] > 0)[0][:5]]; break
        pm = VS['cnt'][:, 0] == nv
        if pm.any(): minlo = min(minlo, float(VS['mmin'][pm].min()))
        if VS['nunres'] > 2_000_000: raise RuntimeError("too many unresolved")
        for i, m in VS['unres']: unres_all.append((int(ql[i]), int(m)))
    T['vertex_route'] = time.time() - t
    out['n_vertex_certified_negative_q'] = n_neg
    if n_neg:
        out.update(certified=False, counterexample_q=negq, seconds_total=time.time() - t00, timing_seconds=T)
        return out
    t = time.time()
    def vpt(m):
        x = list(C)
        for t_, j in enumerate(fcs): x[j] = hi[t_] if (m >> t_) & 1 else lo[t_]
        return x
    vals = ce.exact_margins([vpt(m) for q, m in unres_all], [q for q, m in unres_all], sign, nproc=NPROC) if unres_all else []
    negx = [(q, m, v) for (q, m), v in zip(unres_all, vals) if v < 0]
    T['exact'] = time.time() - t
    out.update(n_vertex_pairs_exact=len(unres_all), n_vertex_pairs_exact_zero=sum(1 for v in vals if v == 0),
               n_vertex_pairs_exact_negative=len(negx), vertex_min_lower_bound_all_positive_q=minlo)
    ok = (not negx) and en.meets(ecls, eff) and eff_ok(spc, eff)
    out.update(certified=bool(ok), timing_seconds=T, seconds_total=time.time() - t00)
    if negx: out['counterexample'] = dict(q=negx[0][0], vertex=[S(x) for x in vpt(negx[0][1])], margin=str(negx[0][2]))
    return out
