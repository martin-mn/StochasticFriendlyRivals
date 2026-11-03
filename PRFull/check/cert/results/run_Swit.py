"""Referee P2 (f): is the SI Table 2 S witness c = (1,1/2,1/2,1/2, 1,1/2,1/2,1/2, 1,1/2,1/2,1, 0,1/2,1/2,1/2)
(stratum 1***1***1**10***, dimension 11) an interior or a boundary point of the friendly-rival set?

Criterion (exact).  For a co-player q with N_q(c) = 0 there are two cases.  If N_q is identically 0 on the
stratum (all corner values 0, L3), q imposes nothing near c.  If N_q is not identically 0, then c is a
boundary point: N_q is multilinear, hence harmonic (each d^2/dx_s^2 N_q = 0), so by the strong minimum
principle it cannot have a local minimum 0 at the interior point c without vanishing identically; so every
neighbourhood of c in the open stratum contains x with N_q(x) < 0, i.e. margin(x, q) < 0 (D* > 0), and
x is not competitive.  The script finds every zero-margin co-player at c exactly, splits them by the corner
identity test, and exhibits explicit exact non-competitive points c - delta e_s for delta down to 2^-40.
Run: cd PRFull/check && nice -n 10 python3 cert/results/run_Swit.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
PAT = '1***1***1**10***'
h_ = F(1, 2)
CW = [1, h_, h_, h_, 1, h_, h_, h_, 1, h_, h_, 1, 0, h_, h_, h_]

def _dir_job(args):
    j, sgn, delta, qs = args
    x = list(CW); x[j] = x[j] + sgn*F(delta)
    vals = [(ex.margin(x, ex.g2p(q), -1), q) for q in qs]
    neg = [(v, q) for v, q in vals if v < 0]
    return j, sgn, len(neg), (min(neg) if neg else None)

if __name__ == '__main__':
    t00 = time.time(); T = {}; sign = -1
    assert en.pattern_str(en.pattern_of(CW)) == PAT
    fcs = en.free_coords(PAT)
    # 1. exact margins at c -----------------------------------------------------------------
    t = time.time()
    Sc = en.i_scan(*en.point_box(CW), sign, nthreads=NTH)
    un = np.where(Sc['cls'] == 2)[0]
    vals = ce.exact_margins([CW]*len(un), un, sign, nproc=NPROC)
    Z = [int(q) for q, v in zip(un, vals) if v == 0]
    n_neg = int((Sc['cls'] == -1).sum()) + sum(1 for v in vals if v < 0)
    n_struct = int((Sc['cls'] == 0).sum()); n_pos = NQ - n_struct - len(Z) - n_neg
    selfp = ex.selfplay(CW); ecls = en.efficiency_class(PAT)
    T['exact_at_c'] = time.time() - t
    # 2. identity test on the zero co-players (2^11 corners each) ------------------------------
    t = time.time()
    Cz = en.corner_scan(PAT, sign, np.array(Z), nthreads=NTH); cnt = Cz['cnt']
    ident = (cnt[:, en.C_POS] == 0) & (cnt[:, en.C_NEG] == 0) & (cnt[:, en.C_UNR] == 0) & (cnt[:, en.C_LOW] == 0)
    B = [q for q, i in zip(Z, ident) if not i]
    T['identity_test'] = time.time() - t
    # 3. explicit exterior points: c +- 1/100 e_s against the non-identity zero co-players --------
    t = time.time()
    jobs = [(j, sg, F(1, 100), B) for j in fcs for sg in (-1, +1)]
    with Pool(NPROC) as pool: D = pool.map(_dir_job, jobs)
    T['directions'] = time.time() - t
    negD = [r for r in D if r[3] is not None]
    if not negD: raise SystemExit("no explicit exterior point at delta = 1/100 among the non-identity zeros")
    best = min(negD, key=lambda r: r[3][0])
    j0, s0, _, (v0, q0) = best
    shrink = []
    for k in [7, 10, 20, 40]:
        x = list(CW); x[j0] += s0*F(1, 2**k); shrink.append((f"2^-{k}", str(ex.margin(x, ex.g2p(q0), -1))))
    # 4. the critic's six points c - 0.01 e_s, s in {1,2,5,6,9,10}: exact point certificates ------
    t = time.time(); crit = []
    for j in [1, 2, 5, 6, 9, 10]:
        x = list(CW); x[j] -= F(1, 100)
        P = ce.certify_point(x, wedge='S', nthreads=NTH)
        crit.append(dict(s=j, n_neg=P['n_neg'], min_margin=str(P['min_margin']), min_margin_float=float(P['min_margin']), argmin=P['argmin']))
    T['critic_points'] = time.time() - t
    boundary = len(B) > 0 and n_neg == 0
    out = dict(name='S_witness_boundary', point=[S(x) for x in CW], stratum=PAT, d=len(fcs),
               efficiency_class=ecls, selfplay_exact=[S(x) for x in selfp],
               at_c=dict(n_pos=n_pos, n_struct_zero=n_struct, n_exact_zero_nonstructural=len(Z), n_neg=n_neg),
               zero_coplayers_identities=int(ident.sum()), zero_coplayers_not_identities=len(B),
               non_identity_zero_coplayers_first50=B[:50],
               directions_delta_1_100=[dict(s=j, dir=sg, n_q_negative=n, most_negative=(None if m is None else [str(m[0]), m[1]])) for j, sg, n, m in D],
               explicit=dict(s=j0, dir=s0, q=q0, margin_at_1_100=str(v0), margin_at_smaller_delta=shrink),
               critic_points=crit, boundary_point=boundary, timing_seconds=T, seconds_total=time.time() - t00)
    with open(os.path.join(HERE, 'S_witness_boundary.json'), 'w') as f: json.dump(jsonable(out), f, indent=1)
    L = [f"# S witness of SI Table 2: {'BOUNDARY POINT (proved)' if boundary else 'not shown to be a boundary point'}",
         f"point   {ptstr(out['point'])}", f"stratum {PAT} (dimension {len(fcs)}), efficiency class {ecls}, exact self-play ({', '.join(out['selfplay_exact'])})",
         "",
         f"exact margins -(w_DC - w_CD) at c over the 65536 binary co-players: positive {n_pos}, structural zero {n_struct},",
         f"   exact non-structural zero {len(Z)}, negative {n_neg}  (so c itself is a friendly rival; minimum 0 [{n_struct + len(Z)}])",
         f"of the {len(Z)} non-structural zeros, {int(ident.sum())} are identities (N_q = 0 on the whole stratum) and {len(B)} are NOT:",
         f"   for each of these {len(B)} co-players N_q(c) = 0 while N_q is not identically 0; N_q is multilinear, hence",
         "   harmonic, so by the strong minimum principle c is not a local minimum of N_q: every neighbourhood of c in",
         "   the open stratum contains points with margin < 0.  c is a BOUNDARY point of the friendly-rival set.",
         "",
         "explicit exterior points c + dir * delta * e_s, delta = 1/100 (number of the non-identity zero co-players with",
         "negative exact margin, and the most negative):"]
    for j, sg, n, m in D:
        L.append(f"   s={j:2d} {'+' if sg > 0 else '-'}:  {n:5d}   " + ('' if m is None else f"{float(m[0]):+.4e} (q={m[1]})"))
    L += [f"most negative: s={j0}, dir {'+' if s0 > 0 else '-'}, q={q0}: exact margin {v0} = {float(v0):.4e};",
          "   the same direction at smaller delta (exact): " + ', '.join(f"{a}: {float(F(b)):.3e}" for a, b in shrink),
          "",
          "the critic's points c - 0.01 e_s (exact point certificates over all 65536 co-players):"]
    dmap = {(j, sg): (n, m) for j, sg, n, m in D}
    for c in crit:
        n_dir, m_dir = dmap[(c['s'], -1)]
        c['all_negatives_among_the_zero_coplayers'] = (n_dir == c['n_neg'])
        c['exact_global_min'] = str(m_dir[0]) if n_dir == c['n_neg'] else None
        L.append(f"   s={c['s']:2d}: {c['n_neg']} co-players beat it (all among the {len(B)} zero co-players: {n_dir == c['n_neg']}); "
                 f"exact minimum margin {c['exact_global_min']} = {float(F(c['exact_global_min'])):.4e}")
    with open(os.path.join(HERE, 'S_witness_boundary.json'), 'w') as f: json.dump(jsonable(out), f, indent=1)
    L += ["", "TIMING (s, <= 3 threads/processes, nice 10): " + ', '.join(f"{k} {v:.1f}" for k, v in T.items()) + f"; total {out['seconds_total']:.1f}"]
    txt = '\n'.join(L) + '\n'; open(os.path.join(HERE, 'S_witness_boundary.txt'), 'w').write(txt); print(txt)
