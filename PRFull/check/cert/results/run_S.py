"""Referee P2 (e): positive-volume certificate in wedge S (T<S, below the switch line), 15-dimensional.
Self-contained replacement of the triage's scratch certS_ratio.py (which imported PartnersRivalsM2f).
Centre c = (32,27,15,25,29,28,9,6,21,14,27,15,6,12,30,22)/32 in the stratum 1*************** (p0 = 1, all
other 15 coordinates free), box c +- 1/1024.  RATIO ROUTE: on the stratum each A*_k(x) is a nonnegative
combination of monomials prod_s f_s(x_s), f_s in {x_s, 1-x_s, 1}, so on the box A*_k(x)/A*_k(c) lies in
[L, U] with L = prod_s (1 - h/m_s), U = prod_s (1 + h/m_s), m_s = min(c_s, 1-c_s); hence w_CD(x)/w_DC(x)
>= rho(c) L/U and margin >= 0 on the box whenever w_DC = 0 structurally or rho(c) = w_CD(c)/w_DC(c) >= U/L.
Three independent evaluations: (1) certify.certify_ratio_box (interval screen, exact where needed);
(2) a fully exact mpq scan of all 65536 co-players at c (the scratch script's computation, 3 processes);
(3) the exact point certificate at c.  Also a larger, proportional box h_s = t m_s.
Run: cd PRFull/check && nice -n 10 python3 cert/results/run_S.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
NUM = [32, 27, 15, 25, 29, 28, 9, 6, 21, 14, 27, 15, 6, 12, 30, 22]
C = [F(n, 32) for n in NUM]
H = F(1, 1024)

def _exact_freq_chunk(qs):
    r = []
    for q in qs:
        W = ex.freqs(C, ex.g2p(q)); r.append((q, str(W[1]), str(W[2])))
    return r

def full_exact_scan(c):
    with Pool(NPROC) as pool:
        res = [x for ch in pool.map(_exact_freq_chunk, [list(range(a, a + 2048)) for a in range(0, NQ, 2048)]) for x in ch]
    nz = nfav = neg = ties = 0; rho = None; arg = None
    for q, a, b in res:
        wcd, wdc = F(a), F(b)                    # S: favourable CD, unfavourable DC
        if wcd == 0 and wdc == 0: nz += 1; continue
        if wdc == 0: nfav += 1; continue
        r = wcd/wdc
        if r < 1: neg += 1
        if r == 1: ties += 1
        if rho is None or r < rho: rho, arg = r, q
    return dict(n_struct_zero=nz, n_unfav_zero=nfav, n_negative=neg, n_ties=ties, rho_min=rho, argmin=arg)

def UL_of(c, lo, hi, fcs):
    L = F(1); U = F(1)
    for j, a, b in zip(fcs, lo, hi):
        L *= min(a/c[j], (1 - b)/(1 - c[j])); U *= max(b/c[j], (1 - a)/(1 - c[j]))
    return U/L

if __name__ == '__main__':
    t00 = time.time(); T = {}
    pat = en.pattern_of(C); fcs = en.free_coords(pat); assert en.pattern_str(pat) == '1' + '*'*15
    t = time.time(); R = ce.certify_ratio_box(C, H, wedge='S', nthreads=NTH); T['ratio_box'] = time.time() - t
    t = time.time(); X = full_exact_scan(C); T['full_exact_scan'] = time.time() - t
    t = time.time(); P = ce.certify_point(C, wedge='S', nthreads=NTH); T['centre_point'] = time.time() - t
    spc = ex.selfplay(C)
    lo = [C[j] - H for j in fcs]; hi = [C[j] + H for j in fcs]; UL = UL_of(C, lo, hi, fcs)
    vol = (2*H)**15
    exact_ok = (X['n_negative'] == 0 and X['n_ties'] == 0 and X['rho_min'] >= UL)
    # proportional box: h_s = t m_s, largest t = k/4096 with U/L <= rho_min
    m = {j: min(C[j], 1 - C[j]) for j in fcs}; rho = X['rho_min']; k = 1
    while UL_of(C, [C[j] - F(k + 1, 4096)*m[j] for j in fcs], [C[j] + F(k + 1, 4096)*m[j] for j in fcs], fcs) <= rho: k += 1
    tp = F(k, 4096); lo2 = [C[j] - tp*m[j] for j in fcs]; hi2 = [C[j] + tp*m[j] for j in fcs]
    t = time.time(); R2 = ce.certify_ratio_box(C, lo=lo2, hi=hi2, wedge='S', nthreads=NTH); T['ratio_box_proportional'] = time.time() - t
    vol2 = F(1)
    for a, b in zip(lo2, hi2): vol2 *= b - a
    # exact spot check inside the box against the argmin co-player and random ones (not needed for the proof)
    rng = random.Random(11); nck = 0; nneg = 0
    for q in [X['argmin']] + rng.sample(range(NQ), 20):
        for _ in range(3):
            x = list(C)
            for j, a, b in zip(fcs, lo2, hi2): x[j] = a + (b - a)*F(rng.randint(0, 64), 64)
            nck += 1; nneg += ex.margin(x, ex.g2p(q), -1) < 0
    out = dict(name='S_box', wedge='S', stratum=en.pattern_str(pat), d=15, centre=[S(x) for x in C], h=S(H),
               box={f"p{j}": [S(a), S(b)] for j, a, b in zip(fcs, lo, hi)}, volume=str(vol), volume_float=float(vol),
               U_over_L=str(UL), U_over_L_float=float(UL),
               efficiency=dict(class_of_open_stratum=R['efficiency_class'], selfplay_centre_exact=[S(x) for x in spc]),
               ratio_box={k: (str(v) if isinstance(v, F) else v) for k, v in R.items()},
               full_exact_scan={k: (str(v) if isinstance(v, F) else v) for k, v in X.items()},
               full_exact_certified=exact_ok,
               centre_point=dict(n_pos=P['n_pos'], n_struct_zero=P['n_struct_zero'], n_exact_zero=P['n_exact_zero'],
                                 n_neg=P['n_neg'], min_margin=S(P['min_margin']), n_zero=P['n_zero']),
               proportional_box=dict(t=str(tp), box={f"p{j}": [S(a), S(b)] for j, a, b in zip(fcs, lo2, hi2)},
                                     volume=str(vol2), volume_float=float(vol2), volume_ratio_to_symmetric=float(vol2/vol),
                                     certified=R2['certified'], U_over_L=str(R2['U_over_L']), n_ratio_interval=R2['n_ratio_interval'],
                                     n_ratio_exact=R2['n_ratio_exact'], n_fail=R2['n_fail']),
               spot_check=f"{nck - nneg}/{nck} exact interior margins >= 0",
               timing_seconds=T, seconds_total=time.time() - t00)
    out['certified'] = bool(R['certified'] and exact_ok and R['efficiency_class'] == 'C' and spc[0] == 1)
    with open(os.path.join(HERE, 'S_box.json'), 'w') as f: json.dump(jsonable(out), f, indent=1)
    Lines = [f"# S_box: {'CERTIFIED' if out['certified'] else 'NOT CERTIFIED'}",
      "STATEMENT.  Every p with p0 = 1 and |p_s - c_s| <= 1/1024 for s = 1..15 (c below) is a friendly rival of",
      "wedge S in the eps -> 0 limit: self-play (1,0,0,0) and w_CD >= w_DC against every binary memory-two co-player",
      "(hence every memory-two co-player).  The box is 15-dimensional, the dimension of the facet p0 = 1.",
      f"   Also certified: the larger proportional box |p_s - c_s| <= {tp} * min(c_s, 1-c_s) (volume {float(vol2):.3g}, "
      f"{float(vol2/vol):.3g} x the symmetric one): {R2['certified']}.",
      "",
      f"stratum {out['stratum']} (dimension 15)",
      f"centre  {ptstr(out['centre'])}",
      f"box     c_s +- 1/1024, s = 1..15; volume (1/512)^15 = {float(vol):.3g}",
      "",
      f"EFFICIENCY (exact, L1): class {R['efficiency_class']} on the whole open stratum; exact self-play at c = ({', '.join(S(x) for x in spc)})",
      "",
      "COMPETITIVENESS (ratio route), all 65536 binary co-players:",
      f"   U/L = {UL} = {float(UL):.6f}",
      f"   (1) interval screen + exact:  both w_CD, w_DC structurally 0: {R['n_struct_zero']};  w_DC structurally 0 < w_CD: {R['n_unfav_struct_zero']};",
      f"       rho(c) >= U/L by interval lower bound: {R['n_ratio_interval']};  by exact rational: {R['n_ratio_exact']};  identities: {R['n_identities']};  failing: {R['n_fail']}",
      f"       rho_min = {R['rho_min']} = {float(R['rho_min']):.6f} at co-player {R['rho_argmin']}",
      f"   (2) fully exact mpq scan at c (no floating point): structural zeros {X['n_struct_zero']}, w_DC = 0 < w_CD {X['n_unfav_zero']},",
      f"       rho < 1: {X['n_negative']}, ties w_CD = w_DC > 0: {X['n_ties']}, rho_min = {X['rho_min']} at {X['argmin']}; rho_min >= U/L: {X['rho_min'] >= UL}",
      f"   (3) exact point certificate at c: positive {P['n_pos']}, structural zero {P['n_struct_zero']}, exact non-structural zero {P['n_exact_zero']}, negative {P['n_neg']};",
      f"       minimum margin {P['min_margin']} attained by {P['n_zero']} co-players",
      "   Number of co-players certified positive by interval: " + str(R['n_ratio_interval'] + R['n_unfav_struct_zero']) +
      f" (of which {R['n_unfav_struct_zero']} have w_DC = 0 structurally); resolved exactly: {R['n_ratio_exact']}; structurally zero: {R['n_struct_zero']}.",
      f"   Smallest ratio margin: rho_min - U/L = {float(R['rho_min'] - UL):.6f} (> 0).  No vertex enumeration (d = 15).",
      f"   Exact spot check (not needed): {out['spot_check']} (argmin co-player + 20 random, proportional box).",
      "",
      "TIMING (s, <= 3 threads/processes, nice 10): " + ', '.join(f"{k} {v:.1f}" for k, v in T.items()) + f"; total {out['seconds_total']:.1f}"]
    txt = '\n'.join(Lines) + '\n'; open(os.path.join(HERE, 'S_box.txt'), 'w').write(txt); print(txt)
