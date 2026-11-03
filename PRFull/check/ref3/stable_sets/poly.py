"""Exact 2-D polygon utilities over Fractions.  A half-plane is (A, B, C) meaning A*u + B*v + C <= 0.
All polygons are intersected with the box |u|, |v| <= BOX, so they are compact and are the convex hull
of their vertices, which are feasible pairwise intersections of the boundary lines."""
from fractions import Fraction as F
from itertools import combinations

BOX = F(1000)
BOXHP = [(F(1), F(0), -BOX), (F(-1), F(0), -BOX), (F(0), F(1), -BOX), (F(0), F(-1), -BOX)]


def feasible(hps, pt):
    u, v = pt
    return all(A * u + B * v + C <= 0 for (A, B, C) in hps)


def vertices(hps):
    H = [h for h in hps if h[0] != 0 or h[1] != 0]
    consts = [h for h in hps if h[0] == 0 and h[1] == 0]
    if any(C > 0 for (_, _, C) in consts):
        return []
    H = H + BOXHP
    pts = set()
    for (a1, b1, c1), (a2, b2, c2) in combinations(H, 2):
        det = a1 * b2 - a2 * b1
        if det == 0:
            continue
        u = (-c1 * b2 + c2 * b1) / det
        v = (-a1 * c2 + a2 * c1) / det
        if feasible(H, (u, v)):
            pts.add((u, v))
    return sorted(pts)


def affine_dim(pts):
    if not pts:
        return -1
    p0 = pts[0]
    vecs = [(p[0] - p0[0], p[1] - p0[1]) for p in pts[1:]]
    vecs = [w for w in vecs if w != (0, 0)]
    if not vecs:
        return 0
    w0 = vecs[0]
    for w in vecs[1:]:
        if w0[0] * w[1] - w0[1] * w[0] != 0:
            return 2
    return 1


def compare(T, Cl, D=()):
    """Polygons T (truth) and Cl (claim), both restricted to the closed domain D and the box.
    Returns (ok, list of witnesses): a vertex of one polygon that is not in the other."""
    TD = list(T) + list(D)
    CD = list(Cl) + list(D)
    VT = vertices(TD)
    VC = vertices(CD)
    bad = []
    for pt in VT:
        if not feasible(CD, pt):
            bad.append(('true-not-claimed', pt))
    for pt in VC:
        if not feasible(TD, pt):
            bad.append(('claimed-not-true', pt))
    return (not bad), bad, VT, VC
