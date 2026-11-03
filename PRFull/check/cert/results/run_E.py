"""Referee P2 (c): positive-volume certificate in wedge E (T<S, above the switch line).
E is the mirror image of N under C <-> D:  pbar_s = 1 - p_{15-s}, which commutes with the error clamp
(1 - (eps + (1-2eps)p) = eps + (1-2eps)(1-p)), maps the pair state s = 4x+y to 15-s, the co-player q to
qbar, CC <-> DD and CD <-> DC; so w_CD(p,q) = w_DC(pbar,qbar), the N margin +(w_DC - w_CD) at (p,q)
equals the E margin -(w_DC - w_CD) at (pbar,qbar), and self-play (0,1/2,1/2,0) is preserved.  Hence the
N box certificate (N_box.txt) proves the mirrored box below.  This script ALSO certifies the mirrored
box directly (vertex route; identities by the L3/L4 corner test instead of 4.5M exact evaluations) and
checks the mirror correspondence exactly.
Run: cd PRFull/check && nice -n 10 python3 cert/results/run_E.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
import run_N
if __name__ == '__main__':
    t0 = time.time()
    patN, cN, h = run_N.PAT, run_N.CEN, run_N.H
    patE = mirror_pattern(patN); cE = mirror_point(cN)
    fN = en.free_coords(patN); fE = en.free_coords(patE)
    loN = {j: cN[j] - h for j in fN}; hiN = {j: cN[j] + h for j in fN}
    loE = [1 - hiN[15 - s] for s in fE]; hiE = [1 - loN[15 - s] for s in fE]
    table_E = [1, 0, 1, F(2,3), F(2,3), F(1,3), 0, F(1,3), 1, 1, F(2,3), F(2,3), F(2,3), 1, 0, F(1,3)]
    assert patE == '101***0*11***10*' and cE == table_E, (patE, cE)
    # --- mirror correspondence, exact and interval
    rng = random.Random(7); qs = rng.sample(range(NQ), 300)
    exact_ok = all(ex.margin(cN, ex.g2p(q), +1) == ex.margin(cE, ex.g2p(mirror_genome(q)), -1) for q in qs)
    vN = [cN[j] if j not in fN else (hiN[j] if rng.random() < .5 else loN[j]) for j in range(16)]
    vE = mirror_point(vN)
    exact_v_ok = all(ex.margin(vN, ex.g2p(q), +1) == ex.margin(vE, ex.g2p(mirror_genome(q)), -1) for q in qs[:100])
    SN = en.i_scan(*en.point_box(cN), +1, nthreads=NTH); SE = en.i_scan(*en.point_box(cE), -1, nthreads=NTH)
    mg = np.array([mirror_genome(q) for q in range(NQ)])
    cls_ok = bool((SN['cls'] == SE['cls'][mg]).all())
    t_mirror = time.time() - t0
    out = box_certificate('E_box', 'E', patE, cE, loE, hiE, exact_all=False)
    out['timing_seconds']['mirror_checks'] = t_mirror
    out['mirror'] = dict(of='N_box', pattern_N=patN, centre_N=[S(x) for x in cN],
                         exact_centre_margins_equal_300q=exact_ok, exact_vertex_margins_equal_100q=exact_v_ok,
                         interval_class_at_centre_equal_all_q=cls_ok)
    st = ("STATEMENT.  Every p in the closed 8-dimensional box below is a friendly rival of wedge E in the\n"
          "eps -> 0 limit: self-play (0,1/2,1/2,0) and w_CD >= w_DC against every co-player.  The box is the\n"
          "C<->D mirror image of the N box (N_box.txt), i.e. the SI Table 2 E witness +- 1/12, and is proved\n"
          "twice: (i) by the mirror symmetry from the N certificate, (ii) directly, below.")
    extra = ["", "MIRROR CHECKS: exact margins N(c,q) = E(cbar,qbar) for 300 random q at the centres: "
             f"{exact_ok}; at a random vertex for 100 q: {exact_v_ok}; interval class at the centres equal for all 65536 q "
             f"under q -> qbar: {cls_ok}"]
    print(write_results(out, os.path.join(HERE, 'E_box'), st, extra))
