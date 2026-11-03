"""Referee P2 (a), larger box: greedy face-by-face enlargement of the W box around the SI Table 2 witness
(stratum 1*1*001*1*1**0*0).  The symmetric boxes +-5/64 and wider fail (a certified-negative vertex
against co-player 32768), so the enlargement is asymmetric.  Screen: the L3/L4 corner test disposes of
50472 co-players once (N_q >= 0 on the whole stratum); a trial box is accepted if every vertex margin of
the remaining co-players has a strictly positive or structural enclosure.  The final box is then
certified from scratch by the full vertex route with exact resolution (box_certificate, exact_all=True).
Run: cd PRFull/check && nice -n 10 python3 cert/results/run_W_large.py"""
import sys, os, itertools; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
PAT = '1*1*001*1*1**0*0'
CEN = [1, F(1,2), 1, F(1,2), 0, 0, 1, F(1,2), 1, F(1,4), 1, F(1,4), F(1,4), 0, F(1,4), 0]
def verts(lo, hi, fcs):
    nv = 1 << len(fcs); vlo = np.empty((nv, 16)); vhi = np.empty((nv, 16)); b0, b1 = en.point_box(CEN)
    for m in range(nv):
        vlo[m] = b0; vhi[m] = b1
        for t, j in enumerate(fcs):
            vlo[m, j], vhi[m, j] = en.enclose(hi[t] if (m >> t) & 1 else lo[t])
    return vlo, vhi
if __name__ == '__main__':
    t00 = time.time(); sign = +1
    fcs = en.free_coords(PAT); d = len(fcs)
    Cs = en.corner_scan(PAT, sign, None, nthreads=NTH); rest = np.where(~Cs['pos_ok'])[0]
    def ok(lo, hi):
        vlo, vhi = verts(lo, hi, fcs)
        VS = en.vertex_scan(vlo, vhi, sign, rest, nthreads=NTH)
        return bool((VS['cnt'][:, 1] == 0).all() and (VS['cnt'][:, 3] == 0).all())
    h = F(1, 16); lo = [CEN[j] - h for j in fcs]; hi = [CEN[j] + h for j in fcs]
    assert ok(lo, hi)
    trials = 0; log = []
    for step in [F(1, 32), F(1, 64), F(1, 128), F(1, 256)]:
        changed = True
        while changed:
            changed = False
            for t in range(d):
                for side in (0, 1):
                    l2, h2 = list(lo), list(hi)
                    if side == 0: l2[t] -= step
                    else: h2[t] += step
                    if not (0 < l2[t] and h2[t] < 1): continue
                    trials += 1
                    if ok(l2, h2):
                        lo, hi = l2, h2; changed = True; log.append((fcs[t], 'lo' if side == 0 else 'hi', str(step)))
    t_search = time.time() - t00
    vol0 = F(1, 8)**d; vol = F(1)
    for a, b in zip(lo, hi): vol *= b - a
    print(f"greedy: {trials} trials, {len(log)} accepted, {t_search:.0f}s; volume {float(vol):.4g} = {float(vol/vol0):.3g} x the +-1/16 box", flush=True)
    out = box_certificate('W_box_large', 'W', PAT, CEN, lo, hi, exact_all=True)
    out['greedy'] = dict(trials=trials, accepted=log, seconds=t_search, volume_ratio_to_pm_1_16=float(vol/vol0))
    out['timing_seconds']['greedy_search'] = t_search
    st = ("STATEMENT.  Every p in the closed 7-dimensional box below (pinned coordinates as in the stratum) is a\n"
          "friendly rival of wedge W in the eps -> 0 limit (self-play (1,0,0,0); w_DC >= w_CD against every co-player).\n"
          f"The box contains the +-1/16 box of W_box.txt and has {float(vol/vol0):.3g} times its volume; it was found by a\n"
          "greedy face-by-face enlargement (steps 1/32 .. 1/256) and then certified from scratch.")
    print(write_results(out, os.path.join(HERE, 'W_box_large'), st,
                        ["", f"GREEDY SEARCH: {trials} trial boxes ({t_search:.0f}s), {len(log)} accepted face moves: " +
                         ', '.join(f"p{a} {b} x{n}" for (a, b), n in sorted(__import__('collections').Counter((a, b) for a, b, _ in log).items()))]))
