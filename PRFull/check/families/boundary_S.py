#!/usr/bin/env python3
r"""The SI Table 2 witness of wedge S lies on the BOUNDARY of the friendly-rival set (within its
11-dimensional stratum 1***1***1**10***), unlike the W and N (and E) witnesses.

Why (exact): at the S witness p* (every free coordinate 1/2) 6728 co-players have margin exactly 0
without being structural zeros or identities (cert/witnesses/S_SI2.json, category pointwise_zero).
For such q, N_q = D* (w_CD - w_DC) is multilinear on the stratum, vanishes at the interior point p*
and is not identically 0.  A multilinear function that is >= 0 near an interior zero has a local
minimum there; along each coordinate line it is affine, so it is constant along every coordinate
line through the point, and inductively identically 0 near p* -- a contradiction.  Hence every
neighbourhood of p* contains points where the margin against q is < 0: p* is not an interior point
of the friendly-rival set.  (For the W, N and E witnesses every non-structural zero is an identity,
so no such obstruction exists; W and N carry certified boxes, cert.certify.certify_box.)

This script exhibits it explicitly: for every pointwise-zero co-player it looks for an exact
rational point p* +- delta e_j (delta = 1/64, one free coordinate moved) with a strictly negative
exact margin, and reports how many are caught that way; for the rest it tries the 2-coordinate
moves p* +- delta e_i +- delta e_j.  All values are exact rationals.

Run:   cd PRFull/check && nice -n 10 python3 families/boundary_S.py     (3 processes, about 1 min)
"""
import os, sys, json, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
from fractions import Fraction as F
import famlib as fl
from cert import engine as en

WIT = os.path.join(os.path.dirname(HERE), 'cert', 'witnesses', 'S_SI2.json')

def main():
    T = fl.timer()
    J = json.load(open(WIT)); p = [F(x) for x in J['strategy']]; sign = J['sign_T_minus_S']
    pwz = J['competitiveness']['argmin']['pointwise_zero']; fr = J['free_coordinates']
    delta = F(1, 64)
    moves1 = [(j, s) for j in fr for s in (+1, -1)]
    pts1 = []
    for j, s in moves1:
        x = list(p); x[j] += s*delta; pts1.append(x)
    R = fl.exact_scan(pts1, sign, qlist=pwz)
    caught = {}; rest = []
    for g in pwz:
        best = min((R[g][i]['m'], i) for i in range(len(moves1)))
        if best[0] < 0: caught[g] = (moves1[best[1]], best[0])
        else: rest.append(g)
    print("S witness %s, stratum %s: %d pointwise-zero co-players" % (J['strategy'], J['stratum'], len(pwz)))
    print("  one-coordinate moves by +-1/64: %d of them have an exact strictly negative margin at a moved point" % len(caught))
    caught2 = {}
    if rest:
        moves2 = [((i, si), (j, sj)) for i, j in itertools.combinations(fr, 2) for si in (1, -1) for sj in (1, -1)]
        pts2 = []
        for (i, si), (j, sj) in moves2:
            x = list(p); x[i] += si*delta; x[j] += sj*delta; pts2.append(x)
        R2 = fl.exact_scan(pts2, sign, qlist=rest)
        for g in rest:
            best = min((R2[g][k]['m'], k) for k in range(len(moves2)))
            if best[0] < 0: caught2[g] = (moves2[best[1]], best[0])
        print("  two-coordinate moves: %d more" % len(caught2))
    g0 = min(caught) if caught else None
    ex0 = None
    if g0 is not None:
        (j, s), v = caught[g0]; x = list(p); x[j] += s*delta
        ex0 = dict(coplayer=g0, moved_coordinate=j, direction=s, delta=str(delta), point=fl.vec_str(x), exact_margin=str(v),
                   stratum_of_point=en.pattern_str(en.pattern_of(x)), efficiency_class=en.efficiency_class(en.pattern_of(x)))
        print("  example: q = %d, p_%d %s 1/64 -> exact margin %s (%.4g); the moved point is in the same stratum %s, class %s"
              % (g0, j, '+' if s > 0 else '-', v, float(v), ex0['stratum_of_point'], ex0['efficiency_class']))
    out = dict(witness='S_SI2', n_pointwise_zero=len(pwz), delta=str(delta), n_negative_after_one_coordinate_move=len(caught),
               n_negative_after_two_coordinate_move=len(caught2), n_not_exhibited=len(pwz) - len(caught) - len(caught2),
               example=ex0, conclusion="the S witness of SI Table 2 is a boundary point of the friendly-rival set in its stratum",
               seconds=round(T(), 1))
    with open(os.path.join(HERE, 'boundary_S.json'), 'w') as f: json.dump(out, f, indent=1)
    print("  not exhibited by these moves: %d (the argument above still applies to them)   (%.1f s)" % (out['n_not_exhibited'], T()))

if __name__ == '__main__':
    main()
