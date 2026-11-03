"""Referee P2 (f), replacement S row, same stratum as the old SI Table 2 S witness.
The old witness (1,1/2,1/2,1/2, 1,1/2,1/2,1/2, 1,1/2,1/2,1, 0,1/2,1/2,1/2) is a boundary point
(S_witness_boundary.txt): 6728 co-players tie with it and none of the ties is an identity; lowering p1, p2,
p5, p6, p9 or p10 makes it lose.  Raising those six coordinates from 1/2 to 3/4 gives
    c' = (1,3/4,3/4,1/2, 1,3/4,3/4,1/2, 1,3/4,3/4,1, 0,1/2,1/2,1/2)   (stratum 1***1***1**10***, dimension 11),
which is tie-free (every zero margin is structural) and the centre of a certified box c' +- 3/16.
Vertex route with exact resolution (box_certificate, exact_all=True), 2^11 vertices x 65536 co-players.
Run: cd PRFull/check && nice -n 10 python3 cert/results/run_Snew11.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p2lib import *
Q, Hf = F(3, 4), F(1, 2)
CEN = [1, Q, Q, Hf, 1, Q, Q, Hf, 1, Q, Q, 1, 0, Hf, Hf, Hf]
H = F(3, 16)
def defensibility(p, sense):
    """float LP of PartnersRivalsM2f/check/m2lib.py (read-only import): min over all co-player behaviours
    of sense*(W_DC - W_CD) at eps = 0 (the SI's defensibility test, double precision, tolerance 1e-9)"""
    sys.path.insert(0, '/Users/martin/Documents/CL1/PartnersRivalsM2f/check')
    import m2lib
    v, _ = m2lib.occupation_lp(np.array([float(x) for x in p]), sense)
    return float(v)
if __name__ == '__main__':
    pat = en.pattern_str(en.pattern_of(CEN)); fcs = en.free_coords(pat); assert pat == '1***1***1**10***'
    out = box_certificate('S_new11_box', 'S', pat, CEN, [CEN[j] - H for j in fcs], [CEN[j] + H for j in fcs],
                          exact_all=True, n_float=100)
    t = time.time()
    old = [1, Hf, Hf, Hf, 1, Hf, Hf, Hf, 1, Hf, Hf, 1, 0, Hf, Hf, Hf]
    dnew = defensibility(CEN, -1); dold = defensibility(old, -1)
    out['defensibility_float_lp'] = dict(new_centre=dnew, old_witness=dold, seconds=time.time() - t)
    st = ("STATEMENT.  Every p in the closed 11-dimensional box below (p0 = p4 = p8 = p11 = 1, p12 = 0; the six\n"
          "coordinates p1,p2,p5,p6,p9,p10 in [9/16,15/16]; p3,p7,p13,p14,p15 in [5/16,11/16]) is a friendly rival\n"
          "of wedge S in the eps -> 0 limit: self-play (1,0,0,0), w_CD >= w_DC against every co-player.  Its centre\n"
          "c' replaces the boundary point of SI Table 2 (same stratum, six entries 1/2 -> 3/4).")
    extra = ["", f"DEFENSIBILITY (float LP of M2f m2lib.occupation_lp, eps = 0, the SI's method; not exact): "
             f"min over all co-player behaviours of (w_CD - w_DC) = {dnew:.3g} at c' (defensible iff >= -1e-9); "
             f"old witness {dold:.3g}"]
    print(write_results(out, os.path.join(HERE, 'S_new11_box'), st, extra))
