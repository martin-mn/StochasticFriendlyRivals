#!/usr/bin/env python3
r"""EXACT CERTIFICATE for the two-parameter family of ms.tex SI 'The wedge W':

    p(a,b) = (1,1,1,0, 0,0,1,0, 1,a,1,b, 0,0,0,0),   a = p_9 = p_<DC,CD>,  b = p_11 = p_<DC,DD>,

claimed to be a friendly rival of W (T > S, below the switch line: efficient means w_CC = 1)
"precisely for 0 < a <= 1/2 and 0 < b <= 1/2".

THEOREM PROVED HERE (set equality over the closed square).  For (a,b) in [0,1]^2,
    p(a,b) is a friendly rival of W   <=>   (a,b) in (0,1/2] x (0,1/2].
The square [0,1]^2 is the union of 9 strata (a, b each 0, free or 1); each is settled exactly:
  * a = 0 or b = 0 (5 strata): NOT EFFICIENT, by L1 (the support of the self-play limit is fixed by
    integer tree exponents, constant on each stratum; its class is not C), with the exact self-play
    limit at one point of each as a check.
  * the open square (0,1)^2, stratum 111000101*1*0000 (class C, efficient by L1).  For each binary
    co-player q, N_q(a,b) = D*(a,b) (w_DC - w_CD) is BILINEAR in (a,b) with the least tree order E
    constant (structural lemma).  N_q and D* are computed EXACTLY at the grid {1/3,2/3}^2, the
    bilinear interpolant is formed exactly and checked exactly at (1/5,5/7) (and D* too), for all
    65536 co-players.  Then:
      - N_q >= 0 at the four vertices (0,0), (0,1/2), (1/2,0), (1/2,1/2) for EVERY q  =>  N_q >= 0 on
        the closed square [0,1/2]^2 (a bilinear function attains its extremes at vertices)  =>
        margin >= 0 on (0,1/2]^2 against every q: competitive, hence a friendly rival there;
      - co-player q_A = 9088 has N = 1 - 2a EXACTLY, so its margin (1-2a)/D*(a,b) is < 0 whenever
        a > 1/2 (D* > 0 on the open stratum), and q_B = 24576 has N = 2(1 - 2b), < 0 whenever b > 1/2.
        So no point of (0,1)^2 with a > 1/2 or b > 1/2 is competitive.
    Cross-checks (independent routes): (i) the vertex values at the boundary vertices (0,0), (0,1/2),
    (1/2,0) are compared, for every q, with the L3 face rule from the exact face-point computation
    (D*_f * margin_f if E_f = E, else 0); (ii) the interval-engine box certificate
    cert.certify.certify_box on [2^-20, 1/2]^2 (corner test, box intervals and vertices);
    (iii) the explicit margins against q_A and q_B are compared with direct exact evaluations at 24
    random rational points.
  * a = 1, 0 < b < 1 and 0 < a < 1, b = 1 (two 1-dimensional strata, class C) and (1,1) (binary,
    class C): NOT COMPETITIVE -- an explicit co-player strictly beats every point, by the exact
    affine interpolation on the edge (values at 1/3, 2/3, extrapolated to the ends, checked at 1/2)
    and by cert.certify.exclude_stratum (L3 + L4, negative direction); at (1,1) by the exact margin.
The tolerance false positive: at a = 1/2 + 2^-40, b = 1/4 the float census rule (min margin
>= -1e-9) accepts the point, but exactly 2496 co-players have a strictly negative margin; the
exact minimum is printed.  At a = 1/2, b = 1/4 every non-structural near-zero margin is exactly 0.

Run:   cd PRFull/check && nice -n 10 python3 families/cert_ab.py    (3 processes, about 3 min)
Writes families/cert_ab.json.
"""
import os, sys, json, time, random
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
from fractions import Fraction as F
import numpy as np
import famlib as fl
from cert import engine as en, exact as ex, certify as ce

SIGN, EFF = en.WEDGES['W']
H = F(1, 2)
PAT = '111000101*1*0000'
QA, QB = 9088, 24576

def P(a, b):
    return [F(1), F(1), F(1), F(0), F(0), F(0), F(1), F(0), F(1), F(a), F(1), F(b), F(0), F(0), F(0), F(0)]

def bilinear(v11, v21, v12, v22):
    """coefficients (al, be, ga, de) of al + be a + ga b + de a b from exact values at
    (1/3,1/3), (2/3,1/3), (1/3,2/3), (2/3,2/3)"""
    c0 = v11; cu = v21 - v11; cw = v12 - v11; cuw = v22 - v21 - v12 + v11      # in u = 3a-1, w = 3b-1
    return (c0 - cu - cw + cuw, 3*cu - 3*cuw, 3*cw - 3*cuw, 9*cuw)

def ev(c, a, b): return c[0] + c[1]*a + c[2]*b + c[3]*a*b

def poly_str(c):
    t = []
    for k, m in zip(c, ('', 'a', 'b', 'ab')):
        if k: t.append('%s%s' % (k, ('*' + m) if m else ''))
    return ' + '.join(t) if t else '0'

def main():
    T = fl.timer(); out = dict(family="p(a,b) = (1,1,1,0,0,0,1,0,1,a,1,b,0,0,0,0)", wedge='W', sign=SIGN,
                               theorem="for (a,b) in [0,1]^2: friendly rival of W <=> 0 < a <= 1/2 and 0 < b <= 1/2")
    assert en.pattern_str(en.pattern_of(P(F(1, 3), F(1, 3)))) == PAT
    # ------------------------------------------------------------ the nine strata: efficiency (L1)
    reps = {0: F(0), 1: F(1, 2), 2: F(1)}           # 0, free (representative 1/2), 1
    strata = {}
    for ia in (0, 1, 2):
        for ib in (0, 1, 2):
            p = P(reps[ia], reps[ib]); pat = en.pattern_str(en.pattern_of(p))
            cls = en.efficiency_class(pat); W = ex.selfplay(p)
            key = 'a=%s,b=%s' % ('0*1'[ia], '0*1'[ib])
            strata[key] = dict(stratum=pat, efficiency_class=cls, efficient=en.meets(cls, EFF),
                               selfplay_at=[fl.fstr(reps[ia]), fl.fstr(reps[ib])], selfplay=fl.vec_str(W))
            assert (W[0] == 1) == en.meets(cls, EFF)
    eff_ok = all(v['efficient'] == ('0' not in k.replace('a=', '').replace('b=', '').split(',')) for k, v in strata.items())
    print("nine strata (L1 efficiency, exact self-play at the representative point):")
    for k, v in strata.items(): print("   %-9s %s  class %s  efficient %-5s  self-play%s" % (k, v['stratum'], v['efficiency_class'], v['efficient'], tuple(v['selfplay'])))
    print("   efficient exactly on the strata with a != 0 and b != 0: %s" % eff_ok)
    out['strata'] = strata; out['efficiency_pattern_ok'] = eff_ok
    # ------------------------------------------------------------ open square: exact bilinear N_q
    t1 = time.time()
    grid = [P(F(1, 3), F(1, 3)), P(F(2, 3), F(1, 3)), P(F(1, 3), F(2, 3)), P(F(2, 3), F(2, 3)), P(F(1, 5), F(5, 7)),
            P(0, 0), P(0, H), P(H, 0), P(H, H)]
    R = fl.exact_scan(grid, SIGN)
    nEv = nbil = nneg = nL3 = 0; coef = {}; Dcoef = {}
    cats = dict(structural_zero=0, identity=0, not_identically_zero=0)
    for g in range(fl.NQ):
        v = R[g]
        if len(set(x['E'] for x in v[:5])) > 1 or len(set(x['supp'] for x in v[:5])) > 1: nEv += 1
        c = bilinear(*[x['N'] for x in v[:4]]); d = bilinear(*[x['D'] for x in v[:4]])
        if ev(c, F(1, 5), F(5, 7)) != v[4]['N'] or ev(d, F(1, 5), F(5, 7)) != v[4]['D']: nbil += 1
        coef[g] = c; Dcoef[g] = d
        vals = [ev(c, 0, 0), ev(c, 0, H), ev(c, H, 0), ev(c, H, H)]
        if min(vals) < 0: nneg += 1
        E = v[0]['E']
        for Nx, f in zip(vals[:3], v[5:8]):          # L3 face rule at the boundary vertices
            pred = f['N'] if f['E'] == E else F(0)
            if f['E'] < E or Nx != pred: nL3 += 1
        if vals[3] != v[8]['N']: nL3 += 1             # (1/2,1/2) is interior: direct value
        if fl.structural_zero(v[0]['supp']): cats['structural_zero'] += 1
        elif all(x == 0 for x in c): cats['identity'] += 1
        else: cats['not_identically_zero'] += 1
    comp_ok = (nEv == 0 and nbil == 0 and nneg == 0 and nL3 == 0)
    print("\nopen square, all 65536 co-players, exact (%.1f s):" % (time.time() - t1))
    print("   E and support constant on the 5 interior points: %s; N_q and D* bilinear (checked at (1/5,5/7)): %s"
          % (nEv == 0, nbil == 0))
    print("   N_q >= 0 at (0,0), (0,1/2), (1/2,0), (1/2,1/2) for every q: %s   (q with a negative vertex: %d)" % (nneg == 0, nneg))
    print("   L3 face rule reproduces the boundary vertex values for every q: %s   (mismatches %d)" % (nL3 == 0, nL3))
    print("   categories: %s" % cats)
    # the binding co-players
    cA, cB = coef[QA], coef[QB]
    bindA = [g for g, c in coef.items() if c[1] < 0 and c[0] == -c[1]*H and c[2] == 0 and c[3] == 0]
    bindB = [g for g, c in coef.items() if c[2] < 0 and c[0] == -c[2]*H and c[1] == 0 and c[3] == 0]
    vanA = [g for g, c in coef.items() if any(c) and ev(c, H, 0) == 0 and ev(c, H, 1) == 0]
    vanB = [g for g, c in coef.items() if any(c) and ev(c, 0, H) == 0 and ev(c, 1, H) == 0]
    sep_ok = (cA == (F(1), F(-2), F(0), F(0)) and cB == (F(2), F(0), F(-4), F(0)))
    print("   q_A = %d: N = %s,  D* = %s" % (QA, poly_str(cA), poly_str(Dcoef[QA])))
    print("   q_B = %d: N = %s,  D* = %s" % (QB, poly_str(cB), poly_str(Dcoef[QB])))
    print("   co-players with N = c(1-2a), c > 0: %d;  with N = c(1-2b): %d;  vanishing on the whole line a=1/2: %d, b=1/2: %d"
          % (len(bindA), len(bindB), len(vanA), len(vanB)))
    # D* > 0 on the open square: the lemma; check the D* polynomials of q_A, q_B at the closed square's corners
    DA = [ev(Dcoef[QA], x, y) for x in (0, 1) for y in (0, 1)]; DB = [ev(Dcoef[QB], x, y) for x in (0, 1) for y in (0, 1)]
    dpos = min(DA) >= 0 and min(DB) >= 0 and max(DA) > 0 and max(DB) > 0
    # direct exact margins at random rational points against q_A, q_B vs the formula
    rng = random.Random(20260927); nform = 0; npts = 24
    for _ in range(npts):
        a = F(rng.randint(1, 999), 1000); b = F(rng.randint(1, 999), 1000)
        for q, c in ((QA, cA), (QB, cB)):
            d = fl.point_data(P(a, b), q, SIGN)
            if d['m'] != ev(c, a, b)/ev(Dcoef[q], a, b) or d['N'] != ev(c, a, b): nform += 1
    print("   explicit margins (1-2a)/D*_A and 2(1-2b)/D*_B = direct exact margins at %d random points: %s;  D*_A, D*_B >= 0 at the corners: %s"
          % (npts, nform == 0, dpos))
    out['open_square'] = dict(stratum=PAT, n_q=fl.NQ, E_support_constant=(nEv == 0), bilinear_checked=(nbil == 0),
                              n_q_negative_at_a_vertex_of_half_square=nneg, L3_face_rule_mismatches=nL3,
                              competitive_on_half_square=comp_ok, categories=cats,
                              qA=dict(g=QA, N=poly_str(cA), Dstar=poly_str(Dcoef[QA]), margin="(1-2a)/D*"),
                              qB=dict(g=QB, N=poly_str(cB), Dstar=poly_str(Dcoef[QB]), margin="2(1-2b)/D*"),
                              n_coplayers_N_c_times_1_minus_2a=len(bindA), n_coplayers_N_c_times_1_minus_2b=len(bindB),
                              n_vanishing_on_line_a_half=len(vanA), n_vanishing_on_line_b_half=len(vanB),
                              formula_checks_failed=nform, separation_ok=(sep_ok and dpos and nform == 0),
                              seconds=round(time.time() - t1, 1))
    # ------------------------------------------------------------ route (ii): interval box certificate
    t2 = time.time()
    lo = [F(1, 2**20)]*2; hi = [H, H]
    rb = ce.certify_box(PAT, lo, hi, wedge='W')
    print("\ninterval-engine box certificate of [2^-20, 1/2]^2 (certify_box): certified=%s  by corners %s, by box %s, by vertices %s  (%.1f s)"
          % (rb['certified'], rb.get('n_by_corners'), rb.get('n_by_box_interval'), rb.get('n_vertex_q'), time.time() - t2))
    out['interval_box_certificate'] = dict(lo=fl.vec_str(lo), hi=fl.vec_str(hi), certified=bool(rb['certified']),
                                           n_by_corners=rb.get('n_by_corners'), n_by_box_interval=rb.get('n_by_box_interval'),
                                           n_vertex_q=rb.get('n_vertex_q'), n_vertex_exact=rb.get('n_vertex_exact'))
    # control: a box reaching past a = 1/2 must fail
    rc = ce.certify_box(PAT, [F(1, 4), F(1, 4)], [F(1, 2) + F(1, 64), F(1, 2)], wedge='W')
    print("   control [1/4, 1/2+1/64] x [1/4, 1/2]: certified=%s, counterexample q=%s" % (rc['certified'], rc.get('counterexample_q')))
    out['interval_box_control'] = dict(certified=bool(rc['certified']), counterexample_q=rc.get('counterexample_q'))
    # ------------------------------------------------------------ the pinned edges a = 1 and b = 1, and (1,1)
    t3 = time.time(); edges = {}
    for key, mk in (('a=1,b=*', lambda s: P(1, s)), ('a=*,b=1', lambda s: P(s, 1))):
        pat = en.pattern_str(en.pattern_of(mk(H)))
        Re = fl.exact_scan([mk(F(1, 3)), mk(F(2, 3)), mk(H)], SIGN)
        best = None
        for g in range(fl.NQ):
            x, y, z = Re[g]
            if x['E'] != y['E']: raise RuntimeError("E varies on an open edge")
            s = (y['N'] - x['N'])*3; N0 = x['N'] - s/3; N1 = N0 + s
            if z['N'] != N0 + s/2: raise RuntimeError("affinity fails on an edge")
            if N0 <= 0 and N1 <= 0 and z['N'] < 0:
                cand = (z['m'], g, N0, N1)
                if best is None or cand < best: best = cand
        exq = ce.exclude_stratum(pat, wedge='W', candidates=[best[1]] if best else None)
        edges[key] = dict(stratum=pat, efficiency_class=en.efficiency_class(pat),
                          excluding_q=best[1] if best else None, N_at_ends=[fl.fstr(best[2]), fl.fstr(best[3])] if best else None,
                          margin_at_half=fl.fstr(best[0]) if best else None,
                          exclude_stratum=dict(excluded=bool(exq['excluded']), q=exq['q']))
        print("edge %s (%s, class %s): co-player %s has N <= 0 at both ends (%s, %s) and margin %s < 0 at 1/2 -> beats every point; "
              "exclude_stratum: excluded=%s by q=%s"
              % (key, pat, edges[key]['efficiency_class'], best[1] if best else None, best[2] if best else '', best[3] if best else '',
                 best[0] if best else None, exq['excluded'], exq['q']))
    R11 = fl.exact_scan([P(1, 1)], SIGN)
    m11 = {g: R11[g][0]['m'] for g in range(fl.NQ)}
    mn11 = min(m11.values()); arg11 = min(g for g, m in m11.items() if m == mn11)
    nneg11 = sum(1 for m in m11.values() if m < 0)
    cls11 = en.efficiency_class(en.pattern_of(P(1, 1)))
    edges['a=1,b=1'] = dict(stratum=en.pattern_str(en.pattern_of(P(1, 1))), efficiency_class=cls11,
                            exact_min_margin=fl.fstr(mn11), argmin=arg11, n_neg=nneg11)
    print("corner (1,1): class %s, exact min margin %s at co-player %d (%d co-players negative)" % (cls11, mn11, arg11, nneg11))
    edges_ok = all(v.get('exclude_stratum', {}).get('excluded', True) and (v.get('excluding_q') is not None or k == 'a=1,b=1')
                   for k, v in edges.items()) and nneg11 > 0
    out['pinned_edges'] = edges; out['pinned_edges_ok'] = edges_ok; out['pinned_edges_seconds'] = round(time.time() - t3, 1)
    # ------------------------------------------------------------ the tolerance false positive
    t4 = time.time()
    sl = {}
    for name, a in (('a=1/2+2^-40', H + F(1, 2**40)), ('a=1/2', H)):
        p = P(a, F(1, 4))
        Rs = fl.exact_scan([p], SIGN)
        ms = {g: Rs[g][0]['m'] for g in range(fl.NQ)}
        neg = sorted(g for g, m in ms.items() if m < 0)
        struct = sum(1 for g in range(fl.NQ) if fl.structural_zero(Rs[g][0]['supp']))
        zero_ns = sum(1 for g, m in ms.items() if m == 0) - struct
        mn = min(ms.values()); arg = min(g for g, m in ms.items() if m == mn)
        # float engine and census rule
        S = en.f_scan(np.array([float(x) for x in p]))
        fm = SIGN*(S['W'][:, 2] - S['W'][:, 1])
        I = en.i_scan(*en.point_box(p), sign=SIGN)
        pred_ok = all(ms[g] == ev(coef[g], a, F(1, 4))/ev(Dcoef[g], a, F(1, 4)) for g in range(fl.NQ))
        sl[name] = dict(point=fl.vec_str(p), exact_min_margin=fl.fstr(mn), exact_min_margin_float=float(mn), argmin=arg,
                        n_negative_exact=len(neg), negative_coplayers=neg, n_structural_zero=struct, n_nonstructural_exact_zero=zero_ns,
                        float_min_margin=float(fm.min()), census_rule_accepts=bool(fm.min() >= -1e-9),
                        interval_counts=dict(pos=int((I['cls'] == 1).sum()), neg=int((I['cls'] == -1).sum()),
                                             struct=int((I['cls'] == 0).sum()), unresolved=int((I['cls'] == 2).sum())),
                        exact_equals_bilinear_formula=pred_ok)
        print("%s, b=1/4: exact min margin %s (= %.6e) at q=%d, %d co-players exactly negative, %d non-structural exact zeros; "
              "float min %.4e -> census rule (>= -1e-9) accepts: %s; interval scan %s; exact = bilinear formula for all q: %s"
              % (name, mn if len(str(mn)) < 60 else '(long rational)', float(mn), arg, len(neg), zero_ns, fm.min(), fm.min() >= -1e-9,
                 sl[name]['interval_counts'], pred_ok))
    sliver_ok = (sl['a=1/2+2^-40']['n_negative_exact'] > 0 and sl['a=1/2+2^-40']['census_rule_accepts']
                 and sl['a=1/2']['n_negative_exact'] == 0)
    out['tolerance_false_positive'] = sl; out['tolerance_false_positive_settled'] = sliver_ok
    out['tolerance_seconds'] = round(time.time() - t4, 1)
    ok = eff_ok and comp_ok and sep_ok and dpos and nform == 0 and rb['certified'] and not rc['certified'] and edges_ok
    out['certified'] = bool(ok); out['seconds'] = round(T(), 1)
    print("\nCERTIFIED (set equality: friendly rival of W <=> (a,b) in (0,1/2]^2, over the closed square): %s   (%.1f s)" % (ok, T()))
    with open(os.path.join(HERE, 'cert_ab.json'), 'w') as f: json.dump(out, f, indent=1)

if __name__ == '__main__':
    main()
