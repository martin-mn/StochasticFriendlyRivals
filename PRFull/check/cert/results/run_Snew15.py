"""Referee P2 (f), replacement S row in the 15-dimensional stratum 1*************** (the maximal S stratum).
Centre c15 = the 11-dimensional replacement c' (run_Snew11.py) with its four pinned coordinates p4, p8, p11
moved from 1 to 7/8 and p12 from 0 to 1/8:
    c15 = (1,3/4,3/4,1/2, 7/8,3/4,3/4,1/2, 7/8,3/4,3/4,7/8, 1/8,1/2,1/2,1/2).
(1) ratio route (certify_ratio_box) on c15 +- 1/256; (2) hybrid ratio + vertex route (p2lib) on c15 +- 1/128,
if its estimated cost is below the local limit; (3) the exact point certificate at c15 (SI Table 2 entry);
(4) the float defensibility LP of M2f m2lib (the SI's method).
Run: cd PRFull/check && nice -n 10 python3 cert/results/run_Snew15.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
from run_Snew11 import defensibility
Q, Hf, Se, Ei = F(3, 4), F(1, 2), F(7, 8), F(1, 8)
CEN = [1, Q, Q, Hf, Se, Q, Q, Hf, Se, Q, Q, Se, Ei, Hf, Hf, Hf]
if __name__ == '__main__':
    t00 = time.time(); T = {}
    pat = en.pattern_str(en.pattern_of(CEN)); fcs = en.free_coords(pat); assert pat == '1' + '*'*15
    t = time.time(); R1 = ce.certify_ratio_box(CEN, F(1, 256), wedge='S', nthreads=NTH); T['ratio_1_256'] = time.time() - t
    t = time.time(); P = ce.certify_point(CEN, wedge='S', nthreads=NTH); T['centre_point'] = time.time() - t
    h2 = F(1, 128)
    t = time.time(); R2 = hybrid_box_certificate('S_new15_h128', 'S', CEN, [CEN[j] - h2 for j in fcs], [CEN[j] + h2 for j in fcs])
    T['hybrid_1_128'] = time.time() - t
    t = time.time(); dfn = defensibility(CEN, -1); T['defensibility'] = time.time() - t
    spc = ex.selfplay(CEN)
    out = dict(name='S_new15_box', centre=[S(x) for x in CEN], stratum=pat, d=15,
               efficiency=dict(class_of_open_stratum=R1['efficiency_class'], selfplay_centre_exact=[S(x) for x in spc]),
               ratio_box_1_256={k: (str(v) if isinstance(v, F) else v) for k, v in R1.items()},
               hybrid_box_1_128=R2,
               centre_point=dict(n_pos=P['n_pos'], n_struct_zero=P['n_struct_zero'], n_exact_zero=P['n_exact_zero'],
                                 n_neg=P['n_neg'], min_margin=S(P['min_margin']), n_zero=P['n_zero']),
               defensibility_float_lp=dfn, timing_seconds=T, seconds_total=time.time() - t00)
    out['certified_1_256'] = bool(R1['certified']); out['certified_1_128'] = bool(R2.get('certified'))
    with open(os.path.join(HERE, 'S_new15_box.json'), 'w') as f: json.dump(jsonable(out), f, indent=1)
    L = [f"# S_new15_box: c15 +- 1/256 {'CERTIFIED' if R1['certified'] else 'NOT certified'}; c15 +- 1/128 "
         f"{'CERTIFIED' if R2.get('certified') else ('refused (cost)' if R2.get('refused') else 'NOT certified')}",
         "STATEMENT.  Every p with p0 = 1 and |p_s - c15_s| <= h for s = 1..15 is a friendly rival of wedge S in the",
         "eps -> 0 limit (self-play (1,0,0,0), w_CD >= w_DC against every co-player), for the h certified below.",
         "", f"stratum {pat} (dimension 15);  centre c15 = {ptstr(out['centre'])}",
         f"EFFICIENCY (exact, L1): class {R1['efficiency_class']}; exact self-play at c15 ({', '.join(S(x) for x in spc)})",
         "", f"(1) ratio route, h = 1/256: U/L = {float(R1['U_over_L']):.6f}, rho_min = {R1['rho_min']} = {float(R1['rho_min']):.6f} (q={R1['rho_argmin']})",
         f"    structural zeros {R1['n_struct_zero']}, w_DC = 0 < w_CD structurally {R1['n_unfav_struct_zero']}, "
         f"ratio by interval {R1['n_ratio_interval']}, by exact {R1['n_ratio_exact']}, identities {R1['n_identities']}, failing {R1['n_fail']}",
         f"    volume (1/128)^15 = {R1['volume']:.3g}"]
    if R2.get('refused'):
        L.append(f"(2) hybrid, h = 1/128: refused, estimated vertex-route cost {R2['estimated_vertex_seconds']:.0f}s for {R2['n_vertex_route']} co-players")
    else:
        L += [f"(2) hybrid, h = 1/128: U/L = {float(F(R2['U_over_L'])):.4f}; structural zeros {R2['n_struct_zero']}, w_DC = 0 structurally {R2['n_unfav_struct_zero']},",
              f"    ratio route {R2['n_ratio_interval']}, vertex route (all 32768 vertices) {R2['n_vertex_route']} co-players: "
              f"certified-negative {R2.get('n_vertex_certified_negative_q')}, exact vertex pairs {R2.get('n_vertex_pairs_exact')} "
              f"(negative {R2.get('n_vertex_pairs_exact_negative')}), min lower bound {R2.get('vertex_min_lower_bound_all_positive_q')}",
              f"    certified: {R2.get('certified')}; volume (1/64)^15 = {R2['volume_float']:.3g}"]
    L += [f"(3) exact point certificate at c15: positive {P['n_pos']}, structural zero {P['n_struct_zero']}, exact non-structural zero {P['n_exact_zero']}, "
          f"negative {P['n_neg']}; minimum margin {P['min_margin']} [{P['n_zero']}]",
          f"(4) defensibility, float LP (M2f m2lib, eps = 0): min (w_CD - w_DC) over all co-player behaviours = {dfn:.3g}",
          "", "TIMING (s, <= 3 threads, nice 10): " + ', '.join(f"{k} {v:.1f}" for k, v in T.items()) + f"; total {out['seconds_total']:.1f}"]
    txt = '\n'.join(L) + '\n'; open(os.path.join(HERE, 'S_new15_box.txt'), 'w').write(txt); print(txt)
