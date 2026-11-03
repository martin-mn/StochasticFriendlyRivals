"""
The 512-point sunflower and its VORONOI tessellation of the k = 4 disk.

    python3 vor.py            self-test: the cells tile the disk, and how

This is the geometry module for the sunflower block, the counterpart of
`disk.py` (the map) and `square.py` (the raster).  It holds no data and no
process; it turns a dict {ipt: value} into 512 filled polygons.

WHY A SUNFLOWER, STATED CORRECTLY.  The obvious thing to say against the
16 x 8 polar grid of `maps.sector_patches` is that its cells are uneven -- fat
at the rim, slivers at the centre.  THAT IS FALSE, and the self-test below
prints the refutation: the radial bins are r = sqrt((j-1/2)/8), which is an
equal-area rule, so all 128 polar cells have area exactly pi/128 and their
area ratio is 1.00 against the sunflower's 2.34.  On equal area the polar grid
WINS.

What is actually wrong with it is SHAPE, and the self-test measures that too,
as the isoperimetric quotient 4 pi A / P^2 (1 for a disc, 0.907 for a regular
hexagon, small for anything elongated).  A polar cell is elongated radially at
the centre and tangentially at the rim -- IQ 0.43 in the inner ring and 0.39
in the outer against a sunflower median near 0.83 -- and sixteen of them meet
at a single point at the origin, which the eye reads as a pinwheel that is not
in the data.  The sunflower is near-isotropic everywhere, has 512 cells to the
polar grid's 128, and reaches further.  CL/BinM2Ev kit ca6 built it for the
same plane as a Fermat (sunflower) spiral,

    rho_i = rho_max sqrt((i - 1/2)/M),   theta_i = 2 pi i / phi^2 ,

whose Voronoi cells are near-hexagonal -- equal area per point from the sqrt,
and phi^2 the most irrational rotation number, so no spoke forms at any M.
Martin asked for that tessellation here, and the points
are ca6's own, copied verbatim into ../Ref/binm2ev_ca6_points.txt rather than
regenerated, so that the two campaigns are the SAME 512 GAMES and their maps
overlay cell for cell.

TWO THINGS TO KNOW WHEN READING A SUNFLOWER PANEL BESIDE A POLAR ONE.

  * The sunflower reaches FURTHER.  Its cut is r = 38.16 (rho = 0.9946, 99.0%
    of the disk by area); the polar grid stops at r = 15.49 (rho = 0.968,
    93.7%).  The outer annulus 0.968 < rho < 0.995 is sampled only here.
  * Near-equal area ON THE DISK is not equal area in the plane.  The map
    compresses
    ferociously outward -- the outermost cells span tens of units of (u,v) --
    so a sunflower cell is a statement about a DIRECTION at large |(u,v)|,
    which is exactly what the disk exists to show.

The clipping is ca6's method (`CL/BinM2Ev/cannon/ca6/p_disk.py`, function
`cells`), reimplemented here rather than imported so that BinM2Ev stays
read-only from this project: mirror every point in the unit circle to supply
the missing neighbours that would otherwise leave the outer cells unbounded --
exact for a circular boundary and it moves no real point -- then
Sutherland-Hodgman against the circle and push the cut vertices back onto it
so the union has no white seam at the rim.
"""
import math
import os

import numpy as np
from scipy.spatial import Voronoi

import disk

PTS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   os.pardir, "Ref", "binm2ev_ca6_points.txt")
_CACHE = {}


def points(path=PTS):
    """(ipt, u, v) for the 512 games, in file order"""
    if "p" not in _CACHE:
        d = np.loadtxt(path)
        _CACHE["p"] = (d[:, 0].astype(int), d[:, 1].copy(), d[:, 2].copy())
    return _CACHE["p"]


def uv_of(ipt):
    i, u, v = points()
    k = int(np.where(i == ipt)[0][0])
    return float(u[k]), float(v[k])


def index_of(tag):
    """V_pNNN_QQW -> NNN; None for every other tag"""
    if not tag.startswith("V_p"):
        return None
    return int(tag.split("_")[1][1:])


def xy():
    """the 512 points as they sit on the unit disk"""
    _, u, v = points()
    x, y = disk.to_disk(u, v)
    return np.column_stack([x, y])


def clip_disk(poly, nseg=64):
    """Sutherland-Hodgman against the unit circle, approximated from OUTSIDE by
    an nseg-gon, then the cut vertices pushed onto the circle so the cells
    leave no white seam at the rim."""
    if len(poly) < 3:
        return poly
    for k in range(nseg):
        t = 2 * math.pi * k / nseg
        a, b, c = math.cos(t), math.sin(t), 1.0
        out, n = [], len(poly)
        for i in range(n):
            p, q = poly[i], poly[(i + 1) % n]
            sp = a * p[0] + b * p[1] - c
            sq = a * q[0] + b * q[1] - c
            if sp <= 0:
                out.append(p)
            if (sp > 0) != (sq > 0):
                out.append(p + (sp / (sp - sq)) * (q - p))
        poly = np.array(out)
        if len(poly) < 3:
            return poly
    rr = np.hypot(*poly.T)
    hot = rr > 1.0 - 1e-9
    poly[hot] /= rr[hot][:, None]
    return poly


def cells():
    """the 512 Voronoi polygons, clipped to the unit disk, in ipt order"""
    if "c" in _CACHE:
        return _CACHE["c"]
    p = xy()
    r2 = (p ** 2).sum(1)
    mir = p / r2[:, None]                     # inversion in the unit circle
    far = 3.0 * np.column_stack([np.cos(np.linspace(0, 2 * math.pi, 33)[:-1]),
                                 np.sin(np.linspace(0, 2 * math.pi, 33)[:-1])])
    vor = Voronoi(np.vstack([p, mir, far]))
    out = []
    for i in range(len(p)):
        reg = vor.regions[vor.point_region[i]]
        q = vor.vertices[reg] if reg and -1 not in reg else np.empty((0, 2))
        out.append(clip_disk(np.asarray(q, float)))
    _CACHE["c"] = out
    return out


def area(poly):
    if len(poly) < 3:
        return 0.0
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def iq(poly):
    """isoperimetric quotient 4 pi A / P^2: 1 for a disc, 0.9069 for a regular
    hexagon, -> 0 for a sliver.  This, not area, is what separates the
    sunflower from the polar grid."""
    if len(poly) < 3:
        return 0.0
    per = np.hypot(*(poly - np.roll(poly, 1, axis=0)).T).sum()
    return 4 * math.pi * area(poly) / per ** 2 if per > 0 else 0.0


def polar_cells(nang=16, nrad=8):
    """`maps.sector_patches`'s own 16 x 8 grid as polygons, for the comparison
    the docstring makes -- so the claim is measured, not asserted."""
    rb = np.sqrt(np.arange(nrad + 1) / float(nrad))
    out = []
    for ia in range(nang):
        t0, t1 = 2 * math.pi * ia / nang, 2 * math.pi * (ia + 1) / nang
        t = np.linspace(t0, t1, 25)
        for j in range(nrad):
            arc0 = np.column_stack([rb[j] * np.cos(t), rb[j] * np.sin(t)])
            arc1 = np.column_stack([rb[j + 1] * np.cos(t[::-1]),
                                    rb[j + 1] * np.sin(t[::-1])])
            out.append(np.vstack([arc0, arc1]))
    return out


def patches(ax, vals, cmap, norm, lw=0.0, ec="none", zorder=1):
    """fill each game's Voronoi cell with its value.  `vals` is {ipt: value};
    a point not in `vals` is left unfilled, which is how a partial campaign
    shows itself rather than being silently interpolated over."""
    from matplotlib.collections import PolyCollection
    C = cells()
    pol, val = [], []
    for ipt, v in vals.items():
        q = C[int(ipt)]
        if len(q) >= 3:
            pol.append(q); val.append(v)
    pc = PolyCollection(pol, array=np.asarray(val, float), cmap=cmap,
                        norm=norm, linewidths=lw, edgecolors=ec,
                        zorder=zorder)
    ax.add_collection(pc)
    return pc


def patches_rgb(ax, cols, miss="0.85", zorder=1):
    """explicit per-cell colours, the counterpart of `mechmaps._sector_rgb`.
    A point with no entry in `cols` is drawn in `miss` rather than skipped, so
    a hole in the campaign shows as grey instead of as white paper."""
    from matplotlib.collections import PolyCollection
    C = cells()
    pol = [q for q in C if len(q) >= 3]
    fc = [cols.get(i, miss) for i, q in enumerate(C) if len(q) >= 3]
    pc = PolyCollection(pol, facecolors=fc, linewidths=0.0,
                        edgecolors="none", zorder=zorder)
    ax.add_collection(pc)
    return pc


def draw_frame(ax, **kw):
    """the sunflower lives ON the disk, so the frame is disk.draw_frame"""
    return disk.draw_frame(ax, **kw)


def selftest():
    i, u, v = points()
    p = xy()
    C = cells()
    A = np.array([area(q) for q in C])
    r = np.hypot(u, v)
    rho = np.hypot(*p.T)
    print("%d sunflower points" % len(i))
    print("  |(u,v)|  %.4f to %.4f, median %.3f" % (r.min(), r.max(),
                                                    np.median(r)))
    print("  rho      %.4f to %.4f  (the polar grid stops at 0.9682)"
          % (rho.min(), rho.max()))
    print("  cells with fewer than 3 vertices: %d" % int((A == 0).sum()))
    nv = np.array([len(q) for q in C])
    print("  vertices per cell: median %d, mean %.2f  (hexagonal is 6)"
          % (int(np.median(nv)), nv.mean()))
    print("  the cells tile the WHOLE disk, not just out to rho_max: total"
          " area %.5f = %.3f%% of pi" % (A.sum(), 100 * A.sum() / math.pi))
    print("     (the 0.06% short is the 64-gon used to clip the rim, not a"
          " gap between cells)")

    Q = np.array([iq(q) for q in C])
    P = polar_cells()
    PA = np.array([area(q) for q in P])
    PQ = np.array([iq(q) for q in P])
    print("\n  sunflower against the 16x8 polar grid it is meant to improve on:")
    print("     %-22s %-26s %s" % ("", "sunflower (512)", "polar (128)"))
    print("     %-22s %-26s %s"
          % ("area, max/min", "%.2f" % (A.max() / A.min()),
             "%.2f   <- the polar grid WINS here" % (PA.max() / PA.min())))
    print("     %-22s %-26s %s"
          % ("isoperimetric Q, min", "%.3f" % Q.min(), "%.3f" % PQ.min()))
    print("     %-22s %-26s %s"
          % ("isoperimetric Q, med", "%.3f" % np.median(Q),
             "%.3f   <- and loses here" % np.median(PQ)))
    print("     %-22s %-26s %s"
          % ("cells meeting at o", "3 (generic Voronoi)", "16 (a pinwheel)"))
    print("     %-22s %-26s %s"
          % ("outermost rho", "%.4f" % rho.max(), "0.9682"))
    print("  a regular hexagon has Q = %.4f; a disc has 1" % (math.pi * math.sqrt(3) / 6))

    assert (A > 0).all(), "a cell came out empty"
    assert abs(A.sum() - math.pi) < 0.01, "the cells do not tile the disk"
    assert np.median(Q) > 0.75, "the cells are not near-isotropic"
    assert np.median(Q) > np.median(PQ), "no shape gain over the polar grid"
    print("  OK")


if __name__ == "__main__":
    selftest()
