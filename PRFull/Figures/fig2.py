#!/usr/bin/env python3
"""
Figure 2 of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity": the memory-one cube [0,1]^4.

    python3 fig2.py      ->  Fig2.pdf / Fig2.png (here), cache/fig2_partner_disk.npz

Top row: three disks of games in the frame of Figure 3 of PartnersRivals1/Revised2 -- the drawing constants,
conics, guide circles, rim names and colour scales are those of PartnersRivals1/FinalFigures/atomcounts.py and
PartnersRivals1/Figures/disklib.py, imported read-only.
  a  1**  the dimension of the efficient set: 3 below the switch line (Eff_CC, in the facet p1 = 1), 2 above it
          (the alternators, the two sheets (p2, p3) = (0, 1), (1, 0)).  PartnersRivalsM1f/Text/friendly.md section 1.
  b  *1*  the dimension of the stable set: 3 on {u < 1} u {v < 0} u {(1, 0)}, 2 on {u >= 1, v >= 0} minus (1, 0);
          never empty.  PartnersRivalsM1f/Text/nash.md.  The point (1, 0) is drawn as a dot of the dimension-3 colour.
  c  11*  the fraction of the efficient set that is a partner (relative 3-volume below the switch line, relative
          2-area of the two sheets above it): the graded partner measure of PartnersRivalsM1f/Text/partners.md
          section 3, evaluated from the closed form partners_volume.frac_of_efficient of
          PartnersRivalsM1f/check/working (the module fig_partners.py imports for its panel b) on a grid of the disk.
Bottom row: the competitive sets in the two facets that carry them (rivals.md Theorem, friendly.md section 2), and
the friendly rivals, one panel under each disk.
  d  Riv+ = {p4 = 0, p2 + p3 <= 1, p2 < 1}, T > S, in the facet p4 = 0 (axes p1, p2, p3), as fig_rivals.py panel b.
  e  Riv- = {p1 = 1, p2 + p3 >= 1, p3 > 0}, T < S, in the facet p1 = 1, the 3-cube of (p2, p3, p4), drawn like d
     with the free coordinate p4 along the first axis (as fig_rivals.py panel c) so that the two prisms, their
     deleted edge p2 = 1, p3 = 0 and the segment (1, s, 1-s, 0) sit in the same places.
  f  the friendly rivals of wedge S in the same prism, FR_S = Riv- n Eff_CC
     = {p1 = 1, p2 > 0, p3 > 0, p2 + p3 >= 1, (p4 > 0 or max(p2, p3) = 1)} (fig_friendly.py panel b).  FR_S has the
     same volume 1/2 as Riv-; what efficiency removes is the face p4 = 0 except its two edges p2 = 1 and p3 = 1, and
     the edge p2 = 0 (TFT--1011).  Those removed pieces are drawn grey (competitive, not efficient); the deleted edge
     p2 = 1, p3 = 0 stays dotted crimson.
Colours of the bottom row: the atom colours of Figure 4 (atomshares.WINCOL: 001 competitive, 111 friendly rival).
Receipts: this script checks panel c's closed form against deterministic quadrature and reproduces fig_partners.py's
square-grid numbers; fig2_check.py re-derives panel c, the carriers of panel b and the vertex census of d, e from the
raw limit tests of PartnersRivalsM1f/check/m1lib.py (cache/fig2_check.txt); fig2_check_f.py checks panel f piece by
piece (interior, faces, edges of the prism) against the same raw tests at two wedge-S games (cache/fig2_check_f.txt).
Nothing outside PRFull is written (no bytecode either).
Notation: the comments here write p1, p2, p3, p4 for the paper's p_CC, p_CD, p_DC, p_DD, which are what the figure
itself shows (axes of d, e, facet titles, key).
"""
import os
import sys
import warnings

sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FF = "/Users/martin/Documents/CL1/PartnersRivals1/FinalFigures"
M1F = "/Users/martin/Documents/CL1/PartnersRivalsM1f"
sys.path.insert(0, FF)
import atomcounts as C                     # noqa: E402  Figure 3's frame; also puts PartnersRivals1/Figures on the path
import disklib                             # noqa: E402
sys.path.insert(0, os.path.join(M1F, "check", "working"))
import partners_volume as PV               # noqa: E402  the closed form of the graded partner measure

import matplotlib                          # noqa: E402
import matplotlib.pyplot as plt            # noqa: E402
from matplotlib.cm import ScalarMappable   # noqa: E402
from matplotlib.collections import PolyCollection   # noqa: E402
from matplotlib.colors import Normalize    # noqa: E402
from matplotlib.lines import Line2D        # noqa: E402
from matplotlib.patches import Circle, Rectangle    # noqa: E402
from mpl_toolkits.mplot3d import proj3d   # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection   # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ====================================================================== the plane of games, exact polygons
RB, TOL = 2000.0, 1e-5                     # as m1atoms.py: the plane clipped at |u|,|v| <= RB, conic tolerance
BOX = [(-RB, -RB), (RB, -RB), (RB, RB), (-RB, RB)]


def clip(poly, a, b, c, keep):
    """the part of a convex polygon where keep * (a u + b v - c) >= 0"""
    out, n = [], len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        sp = (a * p[0] + b * p[1] - c) * keep
        sq = (a * q[0] + b * q[1] - c) * keep
        if sp >= 0:
            out.append(p)
        if (sp > 0) != (sq > 0) and sp != sq:
            t = sp / (sp - sq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def todisk(poly):
    pts = []
    for i in range(len(poly)):
        pts.extend(disklib.arc(poly[i], poly[(i + 1) % len(poly)], C.phi, TOL))
    return np.asarray(pts, float)


def area(q):
    return 0.5 * abs(np.dot(q[:, 0], np.roll(q[:, 1], -1)) - np.dot(q[:, 1], np.roll(q[:, 0], -1)))


def faces(regions):
    F = [todisk(p) for p, _ in regions]
    DA = np.array([area(q) for q in F])
    assert 0 <= (np.pi - DA.sum()) / np.pi < 1e-4, "the mapped faces do not tile the disk"
    return F, DA, np.array([k for _, k in regions], float)


BELOW, ABOVE = clip(BOX, 1, 1, 1, -1), clip(BOX, 1, 1, 1, +1)                 # u + v <= 1, u + v >= 1
EFF = faces([(BELOW, 3), (ABOVE, 2)])
QUAD = clip(clip(BOX, 1, 0, 1, +1), 0, 1, 0, +1)                              # u >= 1, v >= 0
STAB = faces([(clip(BOX, 1, 0, 1, -1), 3), (clip(clip(BOX, 1, 0, 1, +1), 0, 1, 0, -1), 3), (QUAD, 2)])
for tag, (F, DA, K) in (("a efficient", EFF), ("b stable", STAB)):
    print("%-12s disk share:  dim 3 %.4f %%   dim 2 %.4f %%"
          % (tag, 100 * DA[K == 3].sum() / DA.sum(), 100 * DA[K == 2].sum() / DA.sum()))
SATB = area(todisk(clip(clip(BOX, 1, 0, 1, -1), 0, 1, 0, -1))) / np.pi       # u <= 1, v <= 0 (below the line)
SATA = area(todisk(clip(clip(ABOVE, 1, 0, 0, +1), 0, 1, 0, +1))) / np.pi     # u >= 0, v >= 0, above the line
print("c            disk share where every efficient strategy is a partner (u <= 1, v <= 0): %.4f %%" % (100 * SATB))
print("             disk share where the whole live sheet is (u, v >= 0 above the line, value 1/2): %.4f %%"
      % (100 * SATA))

# ====================================================================== panel c: the graded partner measure
NPIX = 1200
g = (np.arange(NPIX) + 0.5) / NPIX * 2.0 - 1.0
X, Y = np.meshgrid(g, g)                          # rows are y
R2 = X * X + Y * Y
IN = R2 < 1.0
sc = C.K / np.sqrt(1.0 - R2[IN])
FRAC = np.full(X.shape, np.nan)
FRAC[IN] = PV.frac_of_efficient(X[IN] * sc, Y[IN] * sc)
fin = FRAC[IN]
assert np.isfinite(fin).all() and fin.min() > 0.0 and fin.max() <= 1.0 + 1e-12
print("c            disk grid %d^2: min %.3e, max %.6f, disk mean %.4f, share == 1 %.3f %%, share == 1/2 %.3f %%"
      % (NPIX, fin.min(), fin.max(), fin.mean(), 100 * np.mean(fin == 1.0), 100 * np.mean(fin == 0.5)))
os.makedirs(os.path.join(HERE, "cache"), exist_ok=True)
np.savez_compressed(os.path.join(HERE, "cache", "fig2_partner_disk.npz"), x=g, frac=FRAC, K=np.array(C.K))

# the numbers fig_partners.py prints for its panel b, recomputed on its own grid of the square |u|, |v| <= 8
NG, SQ = 1400, 8.0
gq = (np.arange(NG) + 0.5) / NG * (2 * SQ) - SQ
UU, VV = np.meshgrid(gq, gq, indexing="ij")
FQ = PV.frac_of_efficient(UU, VV)
print("check        square grid %d^2 (fig_partners.py's): min %.5f, max %.5f, mean %.5f" % (NG, FQ.min(), FQ.max(), FQ.mean()))
for u0, v0 in ((-3.0, 0.7), (2.5, -4.0), (6.0, -3.0), (-4.0, 6.0), (0.3, 0.3)):   # fig_partners.py's spot checks
    ref = PV.vol_below_quad(u0, v0, 1500) if u0 + v0 < 1 else PV.area_above_quad(u0, v0, 200000) / 2.0
    got = float(PV.frac_of_efficient(np.array(u0), np.array(v0)))
    assert abs(ref - got) < 2e-5, (u0, v0, ref, got)
    print("check        (u,v) = (%5.1f,%5.1f): closed form %.6f, deterministic quadrature %.6f" % (u0, v0, got, ref))

# ====================================================================== layout (inches)
PW, LIM, EXTT, EXTB = C.PW, C.LIM, C.EXTT, 0.03          # no sub-line under the disks, so the bottom extension shrinks
PH = PW * (2 * LIM + EXTT + EXTB) / (2 * LIM)
ML, MR, MB, MT = C.ML, C.MR, 0.20, C.MT
COLW, COLGAP = C.COLW, C.COLGAP
W = ML + 3 * COLW + 2 * COLGAP + MR
KH, H3, ROWG = 1.00, 6.00, 0.30                          # key strip, 3D row, gap between the rows
H = MB + KH + H3 + ROWG + PH + MT
FS = C.FS
fig = plt.figure(figsize=(W, H))
YTOP = MB + KH + H3 + ROWG


def disk_axes(col):
    x = ML + col * (COLW + COLGAP)
    return x, fig.add_axes([x / W, YTOP / H, PW / W, PH / H])


def furniture(ax, letter, head):
    """atomcounts.panel's furniture, verbatim in effect: guides, the three conics, ticks, u/v, letter, head, rim names"""
    ax.plot([-1, 1], [0, 0], color=C.GUIDE_C, lw=C.GUIDE_LW, zorder=5, alpha=C.GUIDE_A)
    ax.plot([0, 0], [-1, 1], color=C.GUIDE_C, lw=C.GUIDE_LW, zorder=5, alpha=C.GUIDE_A)
    for s_, col, dash in C.CONICS:
        ax.plot(s_[:, 0], s_[:, 1], color=col, lw=1.5, zorder=7, linestyle=dash)
    for r in C.GUIDE:
        rr = float(np.hypot(*C.phi([[r, 0.0]])[0]))
        ax.plot(rr * np.cos(C.th), rr * np.sin(C.th), "-", color=C.GUIDE_C, lw=C.GUIDE_LW, zorder=5, alpha=C.GUIDE_A)
        for sg in (-1, 1):
            ax.plot([sg * rr, sg * rr], [-0.016, 0.016], color="0.2", lw=1.4, zorder=8, path_effects=C.PEL)
            ax.plot([-0.016, 0.016], [sg * rr, sg * rr], color="0.2", lw=1.4, zorder=8, path_effects=C.PEL)
            lb = ("%g" % r) if sg > 0 else "$-$%g" % r
            ax.text(sg * rr, -0.048, lb, fontsize=FS * 0.7, color="0.2", ha="center", va="top", zorder=9,
                    path_effects=C.PEL)
            ax.text(-0.048, sg * rr, lb, fontsize=FS * 0.7, color="0.2", ha="right", va="center", zorder=9,
                    path_effects=C.PEL)
    ax.text(1.055, 0.0, "$u$", fontsize=FS * 1.4, ha="left", va="center")
    ax.text(0.0, 1.055, "$v$", fontsize=FS * 1.4, ha="center", va="bottom")
    ax.text(-LIM + 0.01, LIM + 0.055, letter, fontsize=FS * 1.7, fontweight="bold", ha="left", va="bottom")
    ax.text(0.0, LIM + 0.055, head, fontsize=FS * 1.05, color="0.12", ha="center", va="bottom")
    for a0, nm in disklib.NAMES:
        _, _, ro = disklib.arc_text(ax, nm, a0, 1.0, C.NM_FS, C.DPU, disklib.NM_GAP, disklib.NM_TRK)
        assert ro < LIM, "%s reaches r = %.4f" % (nm, ro)
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-(LIM + EXTB), LIM + EXTT)
    ax.set_aspect("equal")
    ax.axis("off")


def bar(x, cm, norm, ticks, lab, disc, label):
    """atomcounts.panel's colour bar, to the right of the disk whose left edge is at x inches"""
    xb = x + PW + C.CBGAP
    yc = YTOP + PH * (LIM + EXTB) / (2 * LIM + EXTT + EXTB)
    cax = fig.add_axes([xb / W, (yc - 0.5 * C.CBH) / H, C.CBW / W, C.CBH / H])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cm), cax=cax, ticks=ticks, spacing="uniform", drawedges=disc)
    cb.ax.set_yticklabels(lab)
    cb.ax.minorticks_off()
    for lb in (cb.ax.get_yticklabels()[0], cb.ax.get_yticklabels()[-1]):
        lb.set_fontweight("bold")
    cb.set_label(label, fontsize=FS * 0.95, labelpad=10)
    cb.ax.tick_params(labelsize=FS * (0.72 if disc else 0.82), length=5, width=1.0)
    cb.outline.set_linewidth(1.0)
    if disc:
        cb.dividers.set_color("white")
        cb.dividers.set_linewidth(0.6)


# ---- a, b: dimension fields, exact faces, Figure 3's discrete scale (atomcounts.atom_scale)
for col, (F, DA, K), letter, head in ((0, EFF, "a", r"$\mathbf{1{*}{*}}$   efficient"),
                                      (1, STAB, "b", r"$\mathbf{{*}1{*}}$   stable")):
    x, ax = disk_axes(col)
    cm, norm, ticks, lab, disc = C.atom_scale(K, DA)
    c = cm(norm(K))
    ax.add_collection(PolyCollection(F, facecolors=c, edgecolors=c, linewidths=0.30, antialiased=True, zorder=1))
    if letter == "b":                                  # the isolated point (u, v) = (1, 0), where the dimension is 3
        pc = C.phi([[1.0, 0.0]])[0]
        ax.plot([pc[0]], [pc[1]], "o", ms=6.5, mfc=cm(norm(3.0)), mec="0.15", mew=0.9, zorder=10)
    furniture(ax, letter, head)
    bar(x, cm, norm, ticks, lab, disc, "Dimension of the set")

# ---- c: the graded partner measure, a raster of the disk clipped to the unit circle
x, ax = disk_axes(2)
CMC, NMC = plt.get_cmap("magma"), Normalize(0.0, 1.0)
im = ax.imshow(FRAC, origin="lower", extent=(-1, 1, -1, 1), cmap=CMC, norm=NMC, interpolation="nearest", zorder=1)
im.set_clip_path(Circle((0, 0), 1.0, transform=ax.transData))
furniture(ax, "c", r"$\mathbf{11{*}}$   efficient and stable: partners")
bar(x, CMC, NMC, [0, 0.25, 0.5, 0.75, 1.0], ["0", "1/4", "1/2", "3/4", "1"], False, "Fraction of the efficient set")

# ====================================================================== d, e: the competitive sets in their facets
COMP = "#f0a48c"                                   # atomshares.WINCOL["001"]: competitive
COMP_E = "#b3232b"                                 # eulerlib.COL_C: the competitive outline of Figure 1
FRC = "#7a4fa8"                                    # atomshares.WINCOL["111"]: efficient, stable and competitive
EXCL = "crimson"
CUBE = [[(0, 0, 0), (1, 0, 0)], [(0, 0, 0), (0, 1, 0)], [(0, 0, 0), (0, 0, 1)], [(1, 1, 1), (0, 1, 1)],
        [(1, 1, 1), (1, 0, 1)], [(1, 1, 1), (1, 1, 0)], [(1, 0, 0), (1, 1, 0)], [(1, 0, 0), (1, 0, 1)],
        [(0, 1, 0), (1, 1, 0)], [(0, 1, 0), (0, 1, 1)], [(0, 0, 1), (1, 0, 1)], [(0, 0, 1), (0, 1, 1)]]
NOTEFF = "0.62"                                    # competitive but not efficient (panel f)
AXEDGES = [((0, 0, 0), (1, 0, 0)), ((1, 0, 0), (1, 1, 0)), ((1, 1, 0), (1, 1, 1))]   # x, y, z


def cube_axes(i, axnames, azim):
    """a 3D axes showing the unit cube of the facet: grey frame, three dark edges named as axes (no ticks)"""
    x0 = ML + i * (COLW + COLGAP)                   # the columns of the disks above
    ax = fig.add_axes([(x0 + 0.04 * COLW) / W, (MB + KH) / H, 0.92 * COLW / W, (H3 - 0.45) / H], projection="3d")
    ax.computed_zorder = False                    # draw in the zorder given below, not by depth
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.set_box_aspect((1, 1, 1), zoom=ZOOM)
    ax.view_init(elev=ELEV, azim=azim)
    ax.set_axis_off()
    ax.add_collection3d(Line3DCollection(CUBE, colors="0.72", linewidths=0.8, zorder=1))
    ax.add_collection3d(Line3DCollection(AXEDGES, colors="0.15", linewidths=1.3, zorder=1.5))
    fig.canvas.draw()
    M = ax.get_proj()

    def disp(p):
        xs, ys, _ = proj3d.proj_transform(p[0], p[1], p[2], M)
        return np.array([xs, ys]), ax.transData.transform((xs, ys))
    _, cen = disp((0.5, 0.5, 0.5))
    for (P, Q), nm in zip(AXEDGES, axnames):
        mid = tuple(0.5 * (np.array(P) + np.array(Q)))
        xy, d = disp(mid)
        u = (d - cen) / np.hypot(*(d - cen))
        LAB.append(ax.annotate(nm, xy=tuple(xy), xycoords="data", xytext=tuple(26 * u), textcoords="offset points",
                               fontsize=FS * 1.3, ha="center", va="center", color="0.1"))
    return x0, ax


def head3d(x0, letter, head):
    yt = (MB + KH + H3 - 0.05) / H
    fig.text((x0 + 0.02) / W, yt, letter, fontsize=FS * 1.7, fontweight="bold", ha="left", va="top")
    fig.text((x0 + 0.5 * COLW) / W, yt - 0.05 / H, head, fontsize=FS * 1.05, color="0.12", ha="center", va="top")


def prism(tri):
    """faces of the prism [0,1] x tri, the free coordinate first: the caps at 0 and 1, then the three walls"""
    caps = [[(t, a, b) for a, b in tri] for t in (0.0, 1.0)]
    walls = []
    for i in range(3):
        (a1, b1), (a2, b2) = tri[i], tri[(i + 1) % 3]
        walls.append([(0.0, a1, b1), (0.0, a2, b2), (1.0, a2, b2), (1.0, a1, b1)])
    return caps, walls


def seg(ax, P, Q, z, **kw):
    ax.plot([P[0], Q[0]], [P[1], Q[1]], [P[2], Q[2]], zorder=z, solid_capstyle="round", **kw)


def excluded(ax, P, Q):
    seg(ax, P, Q, 8, color="white", lw=3.6)
    seg(ax, P, Q, 9, color=EXCL, lw=3.0, ls=(0, (1.6, 1.4)))


def vertices(ax, V):
    """filled/open circles at the vertices, and their labels placed in points around the projected vertex"""
    fig.canvas.draw()
    M = ax.get_proj()
    for (xyz, name, kind, dx, dy, ha) in V:
        fc = {"comp": "k", "fr": FRC, "out": "white"}[kind]
        ax.scatter([xyz[0]], [xyz[1]], [xyz[2]], s=85, c=[fc], edgecolors="k", linewidths=1.3, zorder=20,
                   depthshade=False)
        xs, ys, _ = proj3d.proj_transform(xyz[0], xyz[1], xyz[2], M)
        LAB.append(ax.annotate(name, xy=(xs, ys), xycoords="data", xytext=(dx, dy), textcoords="offset points",
                               fontsize=FS * 1.0, ha=ha, va="center", color="0.08", zorder=21))


ZOOM, ELEV = 1.00, 20
LAB = []

# ---- d: Riv+ in the facet p4 = 0; axes (p1, p2, p3); the prism [0,1] x {p2 + p3 <= 1}
x0, ax = cube_axes(0, ("$p_{CC}$", "$p_{CD}$", "$p_{DC}$"), -56)
caps, walls = prism([(0, 0), (1, 0), (0, 1)])
ax.add_collection3d(Poly3DCollection(caps + walls, facecolors=COMP, alpha=0.40, edgecolors=COMP_E, linewidths=1.0,
                                     zorder=3))
excluded(ax, (0, 1, 0), (1, 1, 0))                                      # (p1, 1, 0, 0): not competitive
seg(ax, (1, 0, 1), (1, 1, 0), 10, color="k", lw=2.6)                    # (1, s, 1-s, 0)
vertices(ax, [((0, 0, 0), "0000 ALLD", "comp", -10, -6, "right"),
              ((1, 0, 0), "1000 Grim", "comp", 11, -5, "left"),
              ((0, 0, 1), "0010", "comp", -10, 0, "right"),
              ((1, 0, 1), "1010 TFT", "comp", 0, 15, "center"),
              ((0, 1, 0), "0100", "out", 9, 10, "left"),
              ((1, 1, 0), "1100", "out", 10, -2, "left")])
head3d(x0, "d", r"$\mathbf{{*}{*}1}$   competitive, $T>S$  facet $p_{DD}=0$")

# ---- e: Riv- in the facet p1 = 1; axes (p4, p2, p3); the prism [0,1] x {p2 + p3 >= 1}, drawn like d
TRI_M = [(1, 0), (0, 1), (1, 1)]
x0, ax = cube_axes(1, ("$p_{DD}$", "$p_{CD}$", "$p_{DC}$"), -56)
caps, walls = prism(TRI_M)
ax.add_collection3d(Poly3DCollection(caps + walls, facecolors=COMP, alpha=0.40, edgecolors=COMP_E, linewidths=1.0,
                                     zorder=3))
excluded(ax, (0, 1, 0), (1, 1, 0))                                      # (1, 1, 0, p4): not competitive
seg(ax, (0, 0, 1), (0, 1, 0), 10, color="k", lw=2.6)                    # (1, s, 1-s, 0)
vertices(ax, [((0, 0, 1), "1010 TFT", "comp", -10, 6, "right"),
              ((1, 0, 1), "1011", "comp", 0, 15, "center"),
              ((0, 1, 1), "1110", "comp", 10, 8, "left"),
              ((1, 1, 1), "1111 ALLC", "comp", 10, 6, "left"),
              ((0, 1, 0), "1100", "out", 9, 10, "left"),
              ((1, 1, 0), "1101", "out", 10, -2, "left")])
head3d(x0, "e", r"$\mathbf{{*}{*}1}$   competitive, $T<S$  facet $p_{CC}=1$")

# ---- f: the friendly rivals of wedge S, FR_S, in the same prism; the pieces efficiency removes in grey
x0, ax = cube_axes(2, ("$p_{DD}$", "$p_{CD}$", "$p_{DC}$"), -56)
caps, walls = prism(TRI_M)
cap0, cap1 = caps
ax.add_collection3d(Poly3DCollection([cap1] + walls, facecolors=FRC, alpha=0.30, edgecolors=COMP_E, linewidths=1.0,
                                     zorder=3))
ax.add_collection3d(Poly3DCollection([cap0], facecolors=NOTEFF, alpha=0.45, edgecolors=COMP_E, linewidths=1.0,
                                     zorder=4))                          # p4 = 0: competitive, not efficient
seg(ax, (0, 1, 0), (0, 1, 1), 7, color=FRC, lw=3.2)                     # (1, 1, p3, 0): friendly rivals
seg(ax, (0, 0, 1), (0, 1, 1), 7, color=FRC, lw=3.2)                     # (1, p2, 1, 0): friendly rivals
seg(ax, (0, 0, 1), (1, 0, 1), 7, color=NOTEFF, lw=4.0)                  # (1, 0, 1, p4): competitive, not efficient
excluded(ax, (0, 1, 0), (1, 1, 0))                                      # (1, 1, 0, p4): not competitive
vertices(ax, [((0, 0, 1), "1010 TFT", "comp", -10, 6, "right"),
              ((1, 0, 1), "1011", "comp", 0, 15, "center"),
              ((0, 1, 1), "1110", "fr", 10, 8, "left"),
              ((1, 1, 1), "1111 ALLC", "fr", 0, 16, "center"),
              ((0, 1, 0), "1100", "out", 9, 10, "left"),
              ((1, 1, 0), "1101", "out", 10, -2, "left")])
head3d(x0, "f", r"$\mathbf{111}$   friendly rivals, wedge $S$  facet $p_{CC}=1$")

# ---- the key of d, e, f
axk = fig.add_axes([ML / W, MB / H, (W - ML - MR) / W, KH / H])
axk.set_xlim(0, W - ML - MR)
axk.set_ylim(0, KH)
axk.axis("off")
ROWS = [[("patch", COMP, "competitive"), ("patch", FRC, "friendly rivals"), ("patch", NOTEFF, "competitive, not efficient"),
         ("dash", EXCL, r"$p_{CD}=1$, $p_{DC}=0$: not competitive")],
        [("dot", "k", "binary, competitive"), ("dot", FRC, "binary friendly rival"),
         ("dot", "white", "binary, not competitive"), ("line", "k", r"$(1,s,1-s,0)$: competitive, not efficient")]]
KFS = FS * 1.0
R = fig.canvas.get_renderer()
for r_, items in enumerate(ROWS):
    yk = KH * (0.72 if r_ == 0 else 0.26)
    xs, arts = 0.0, []
    for kind, colr, text in items:
        if kind == "patch":
            arts.append(axk.add_patch(Rectangle((xs, yk - 0.13), 0.40, 0.26, fc=colr,
                                                alpha={COMP: 0.55, NOTEFF: 0.55}.get(colr, 0.40), ec=COMP_E, lw=0.9)))
        elif kind == "dot":
            arts += axk.plot([xs + 0.20], [yk], "o", ms=9, mfc=colr, mec="k", mew=1.2)
        elif kind == "dash":
            arts += axk.plot([xs, xs + 0.40], [yk, yk], color=EXCL, lw=3.0, ls=(0, (1.6, 1.4)))
        else:
            arts += axk.plot([xs, xs + 0.40], [yk, yk], color=colr, lw=2.6)
        t = axk.text(xs + 0.52, yk, text, fontsize=KFS, ha="left", va="center", color="0.12")
        arts.append(t)
        xs += 0.52 + t.get_window_extent(R).width / fig.dpi + 0.42
    shift = 0.5 * ((W - ML - MR) - (xs - 0.42))             # centre each row of the key
    for art in arts:
        if isinstance(art, Rectangle):
            art.set_x(art.get_x() + shift)
        elif isinstance(art, Line2D):
            art.set_xdata(np.asarray(art.get_xdata()) + shift)
        else:
            xx, yy = art.get_position()
            art.set_position((xx + shift, yy))

# ====================================================================== checks and output
fig.canvas.draw()
R = fig.canvas.get_renderer()
OFF = []
for t_ in fig.texts + [t for a_ in fig.axes if a_.name != "3d" for t in a_.texts] + LAB:
    if t_.get_text() and t_.get_visible():
        bb = t_.get_window_extent(R)
        if not (0 < bb.x0 and bb.x1 < W * fig.dpi and 0 < bb.y0 and bb.y1 < H * fig.dpi):
            OFF.append(t_.get_text())
BB = [t_.get_window_extent(R) for t_ in LAB + fig.texts]
for i_ in range(len(BB)):
    for j_ in range(i_ + 1, len(BB)):
        if BB[i_].overlaps(BB[j_]):
            OFF.append("overlap %r / %r" % ((LAB + fig.texts)[i_].get_text(), (LAB + fig.texts)[j_].get_text()))
out = os.path.join(HERE, "Fig2.pdf")
fig.savefig(out)
fig.savefig(out.replace(".pdf", ".png"), dpi=200)
plt.close(fig)
print("Fig2: %.2f x %.2f in, pdf %.2f MB" % (W, H, os.path.getsize(out) / 1e6))
assert not OFF, "text runs off the sheet: %r" % OFF
