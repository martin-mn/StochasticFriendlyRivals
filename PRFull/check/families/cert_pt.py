#!/usr/bin/env python3
r"""EXACT CERTIFICATE for the one-parameter family of ms.tex (main text 'binary-invisible', Figure 3c,
SI 'The other wedges'):

    p(t) = (0,0,0,t, 0,0,0,0, 0,1,0,1, 1,0,1,0),   the stratum 000*000001011010 (free p_3 = p_<CC,DD>).

CLAIM PROVED.  For EVERY t in the open interval (0,1), p(t) is a friendly rival of the wedge N
(T > S, above the switch line):
  (i)  efficiency: its self-play limit is exactly (0, 1/2, 1/2, 0);
  (ii) competitiveness: w_DC - w_CD >= 0 against every one of the 65536 binary memory-two
       co-players (hence, by the box-vertex lemma, against every co-player of any memory).
At the two binary endpoints it is NOT efficient: self-play (1/6,1/3,1/3,1/6) at t = 0 and
(2/9,1/3,1/3,1/9) at t = 1 (exact).

PROOF ROUTES (all three are run and must agree):
  (1) efficiency by L1: the support of the self-play limit is decided by integer tree exponents, which
      depend on the pattern only, so it is constant on the open stratum; at the centre it is the
      class A support (w_CC = w_DD = 0).  Self-play is symmetric under exchanging the players (the
      pair chain commutes with the mirror map), so w_CD = w_DC and hence both are 1/2.  Checked
      exactly at 11 rational t as well.
  (2) competitiveness by the corner test (L3 + L4, positive direction), cert.certify:
      certify_open_stratum_all_q -- every stratum corner has N_q(c) >= 0 for every q.
  (3) competitiveness by exact affine interpolation, independent of L3: for fixed q, N_q(t) =
      D*(t) (w_DC - w_CD)(t) is affine in t on the open interval (each spanning-tree product uses the
      focal coordinate p_3 at most once, and the least tree order E is constant there).  N_q is
      computed EXACTLY at t = 1/3 and 2/3, checked exactly at t = 1/2 and t = 1/7 (affinity), and
      extrapolated exactly to t = 0 and t = 1: both values >= 0 for every q, so N_q >= 0 on [0,1]
      and the margin N_q/D* >= 0 on (0,1) (D* > 0 there).  As a by-product the extrapolated
      endpoint values are compared, for every q, with the L3 prediction from the exact binary
      corner (D*_c * margin_c if E_c = E, else 0).
Categories per co-player (constant on the open interval): structural zero (w_CD = w_DC = 0 by
exponents), identity (N_q identically 0 but not structural), positive (N_q > 0 on the whole open
interval: affine, >= 0 at both ends, not identically 0).

Run:   cd PRFull/check && nice -n 10 python3 families/cert_pt.py      (3 processes, about 1-2 min)
Writes families/cert_pt.json and prints a log (families/cert_pt.out when redirected).
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
from fractions import Fraction as F
import numpy as np
import famlib as fl
from cert import engine as en, exact as ex, certify as ce

PAT = '000*000001011010'
SIGN, EFF = en.WEDGES['N']

def P(t):
    return [F(0), F(0), F(0), F(t), F(0), F(0), F(0), F(0), F(0), F(1), F(0), F(1), F(1), F(0), F(1), F(0)]

def main():
    T = fl.timer(); out = dict(claim="p(t) = (0,0,0,t,0,0,0,0,0,1,0,1,1,0,1,0) is a friendly rival of wedge N "
                                     "for every t in (0,1); not efficient at t = 0 and t = 1",
                                stratum=PAT, wedge='N', sign=SIGN)
    assert en.pattern_str(en.pattern_of(P(F(1, 2)))) == PAT
    # ---------------------------------------------------------------- (1) efficiency
    cls = en.efficiency_class(PAT)
    ts = [F(1, 1000), F(1, 10), F(1, 7), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(6, 7), F(9, 10), F(999, 1000)]
    sp = {str(t): fl.vec_str(ex.selfplay(P(t))) for t in ts}
    sp_ok = all(v == ['0', '1/2', '1/2', '0'] for v in sp.values())
    ends = {}
    for t in (F(0), F(1)):
        W = ex.selfplay(P(t)); ends[str(t)] = dict(selfplay=fl.vec_str(W), efficiency_class=en.efficiency_class(en.pattern_of(P(t))),
                                                     efficient_N=(W[1] == F(1, 2) and W[2] == F(1, 2)))
    ends_ok = (ends['0']['selfplay'] == ['1/6', '1/3', '1/3', '1/6'] and ends['1']['selfplay'] == ['2/9', '1/3', '1/3', '1/9']
               and not ends['0']['efficient_N'] and not ends['1']['efficient_N'])
    out['efficiency'] = dict(L1_class_open_stratum=cls, required=EFF, selfplay_exact_at=sp, all_interior_exact_equal_0_half_half_0=sp_ok,
                             endpoints=ends, endpoints_not_efficient=ends_ok)
    print("(1) efficiency: L1 class of the open stratum = %s (required %s); exact self-play at %d interior t all (0,1/2,1/2,0): %s"
          % (cls, EFF, len(ts), sp_ok))
    for k, v in ends.items(): print("    t = %s: self-play %s, class %s, efficient: %s" % (k, tuple(v['selfplay']), v['efficiency_class'], v['efficient_N']))
    # ---------------------------------------------------------------- (2) corner test
    t2 = time.time()
    r2 = ce.certify_open_stratum_all_q(PAT, wedge='N', stop_on_fail=False)
    out['corner_test'] = dict(certified=bool(r2['certified']), n_q_checked=r2['n_q_checked'], n_q_failing=r2['n_q_failing'],
                              corner_counts=r2['corner_counts'], seconds=round(time.time() - t2, 3))
    print("(2) corner test (L3+L4): certified=%s  q checked %d, failing %d, corners %s  (%.2f s)"
          % (r2['certified'], r2['n_q_checked'], r2['n_q_failing'], r2['corner_counts'], time.time() - t2))
    # ---------------------------------------------------------------- (3) exact affine interpolation
    t3 = time.time()
    pts = [P(F(1, 3)), P(F(2, 3)), P(F(1, 2)), P(F(1, 7)), P(F(0)), P(F(1))]
    R = fl.exact_scan(pts, SIGN)
    n_E_var = n_aff_bad = n_neg = n_L3_bad = 0
    cat = dict(structural_zero=[], identity=[], positive=[])
    pos_zero_end = dict(t0=0, t1=0)
    min_pos_half = None; arg_pos_half = None
    for g in range(fl.NQ):
        a, b, c, d, c0, c1 = R[g]
        if not (a['E'] == b['E'] == c['E'] == d['E'] and a['supp'] == b['supp'] == c['supp'] == d['supp']): n_E_var += 1
        # N(t) = N0 + (N1 - N0) t, from the exact values at 1/3 and 2/3
        slope = (b['N'] - a['N'])*3; N0 = a['N'] - slope/3; N1 = N0 + slope
        if c['N'] != N0 + slope/2 or d['N'] != N0 + slope/7: n_aff_bad += 1
        # the same for D* (affine too)
        sD = (b['D'] - a['D'])*3; D0 = a['D'] - sD/3
        if c['D'] != D0 + sD/2 or d['D'] != D0 + sD/7: n_aff_bad += 1
        if N0 < 0 or N1 < 0: n_neg += 1
        # L3 prediction at the corners t = 0, 1
        E = a['E']
        for Nx, cc in ((N0, c0), (N1, c1)):
            pred = cc['N'] if cc['E'] == E else F(0)
            if cc['E'] < E or Nx != pred: n_L3_bad += 1
        if fl.structural_zero(a['supp']): cat['structural_zero'].append(g); assert N0 == 0 and N1 == 0
        elif N0 == 0 and N1 == 0: cat['identity'].append(g)
        else:
            cat['positive'].append(g)
            if N0 == 0: pos_zero_end['t0'] += 1
            if N1 == 0: pos_zero_end['t1'] += 1
            if min_pos_half is None or c['m'] < min_pos_half: min_pos_half, arg_pos_half = c['m'], g
    ok3 = (n_E_var == 0 and n_aff_bad == 0 and n_neg == 0 and n_L3_bad == 0)
    counts = {k: len(v) for k, v in cat.items()}
    out['exact_affine'] = dict(certified=ok3, n_q=fl.NQ, n_E_or_support_varies=n_E_var, n_affinity_failures=n_aff_bad,
                               n_negative_endpoint=n_neg, n_L3_mismatch=n_L3_bad, categories=counts,
                               positive_with_N_zero_at_t0=pos_zero_end['t0'], positive_with_N_zero_at_t1=pos_zero_end['t1'],
                               min_positive_margin_at_t_half=fl.fstr(min_pos_half), argmin_positive_at_t_half=arg_pos_half,
                               seconds=round(time.time() - t3, 1))
    print("(3) exact affine interpolation over all 65536 co-players (%.1f s): E/support constant %s, affine %s, "
          "N_q(0), N_q(1) >= 0 for all q: %s, L3 corner rule exact for all q: %s"
          % (time.time() - t3, n_E_var == 0, n_aff_bad == 0, n_neg == 0, n_L3_bad == 0))
    print("    categories on the whole open interval: %s  (zero margin: %d)" % (counts, counts['structural_zero'] + counts['identity']))
    print("    positive co-players whose N_q vanishes at t=0: %d, at t=1: %d (their margin tends to 0 at that end)"
          % (pos_zero_end['t0'], pos_zero_end['t1']))
    print("    smallest positive margin at t=1/2: %s = %.6g (co-player %d)" % (min_pos_half, float(min_pos_half), arg_pos_half))
    # the minimum-margin co-player lists are constant on (0,1)
    out['zero_margin_coplayers'] = dict(structural_zero=cat['structural_zero'], identity=cat['identity'])
    # ---------------------------------------------------------------- endpoints: competitiveness (for the record)
    for t in (F(0), F(1)):
        c = ce.certify_point(P(t), wedge='N')
        ends[str(t)].update(min_margin=fl.fstr(c['min_margin']), n_neg=c['n_neg'], n_zero=c['n_zero'], n_pos=c['n_pos'])
        print("    t = %s: min margin %s, zero at %d, negative at %d co-players (competitive but not efficient)"
              % (t, c['min_margin'], c['n_zero'], c['n_neg']))
    ok = (cls == EFF and sp_ok and ends_ok and r2['certified'] and ok3)
    out['certified'] = bool(ok); out['seconds'] = round(T(), 1)
    print("CERTIFIED: p(t) is a friendly rival of wedge N for every t in (0,1), and not efficient at t = 0, 1: %s  (%.1f s)" % (ok, T()))
    with open(os.path.join(HERE, 'cert_pt.json'), 'w') as f: json.dump(out, f, indent=1)

if __name__ == '__main__':
    main()
