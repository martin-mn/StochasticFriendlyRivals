#!/usr/bin/env python3
"""
Figure 1 of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity": friendly rivals over the plane of games, for the
stochastic memory-one and memory-two strategies -- the analogue of Figure 1b, c and its key in
PartnersRivals1/Revised2/ms.tex.

    python3 fig1.py            ->  Fig1.pdf, Fig1.png (dpi 200), here

a, the disk of games for Memory-1, stochastic, [0,1]^4; b, the same for Memory-2, stochastic, [0,1]^16.  The disk is
the reference's: (u,v) -> (u,v)/sqrt(16+u^2+v^2), rim names, tick labels, faint circles, the switch line u+v=1 crimson
dashed, T=S white dashed, the donation game u+v=0, v>0 cyan dash-dotted -- drawn by PartnersRivals1/FinalFigures/
atomshares.furniture, imported read-only.  Each open wedge W, S, N, E is coloured by the dimension of its friendly-
rival set, on one discrete magma scale shared by both disks and the key (colours by rank, as the reference's discrete
bars); grey where the set is empty (the reference's grey 0.84).  The friendly-rival set is constant on each open wedge
(reference SI "Where the classes live in the plane of games"), so a disk is exactly four regions.  Below, the key: one
row per locus, the five boundary loci included, which have no area on the disks.  All numbers come from fig1data.py,
which takes the table values from fig1_reference_values.json and recomputes what it can (python3 fig1data.py
prints the provenance; python3 fig1data.py --verify-sources checks the values against the manuscripts, which only the
author's machine has).
"""
import os
import sys

sys.dont_write_bytecode = True                                   # import the reference modules without writing there
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivals1/FinalFigures")

import numpy as np                                               # noqa: E402
import matplotlib                                                # noqa: E402
matplotlib.use("Agg")
import matplotlib.patheffects as pe                              # noqa: E402
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.cm import ScalarMappable                         # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, to_rgb   # noqa: E402
from matplotlib.patches import Polygon, Rectangle                # noqa: E402

import atomshares as A                                           # noqa: E402  furniture, sizes, grey (read-only)
import fig1data as D                                             # noqa: E402

D.check(verbose=True)
K = A.K
LOCI = D.LOCI
WEDGES = ["W", "S", "N", "E"]

# ------------------------------------------------------------------ the shared discrete scale
DIMS = sorted({D.M1[k]["dim"] for k in LOCI if D.M1[k]["dim"] is not None} | {D.M2[k]["dim"] for k in LOCI})
assert DIMS == [2, 3, 7, 8, 14, 15], DIMS
CMAP = plt.get_cmap("magma")
COL = {d: CMAP(t) for d, t in zip(DIMS, np.linspace(0.0, 1.0, len(DIMS)))}
GREY = A.GREY


def colour(dim):
    return GREY if dim is None else COL[dim]


def ink(c):
    r, g, b = to_rgb(c)
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.45 else "0.08"


# ------------------------------------------------------------------ the four open wedges as disk polygons
def phi(p):
    return A.phi(p)


def ray(d):
    """the image of the ray (1,0) + s d, s from 0 to far out, as a disk polyline starting at the apex"""
    s = np.concatenate([np.linspace(0.0, 20.0, 800, endpoint=False), np.geomspace(20.0, 2.0e6, 800)])
    return phi(np.array([1.0, 0.0])[None, :] + s[:, None] * np.asarray(d, float)[None, :])


def rim(a0, a1, n=400):
    t = np.deg2rad(np.linspace(a0, a1, n))
    return np.column_stack([np.cos(t), np.sin(t)])


SW_UP, SW_DN, TS_UP, TS_DN = (-1, 1), (1, -1), (1, 1), (-1, -1)   # switch line: to 135 / 315 degrees; T=S: to 45 / 225
WPOLY = {"W": np.vstack([ray(SW_UP), rim(135, 225), ray(TS_DN)[::-1]]),
         "N": np.vstack([ray(TS_UP), rim(45, 135), ray(SW_UP)[::-1]]),
         "E": np.vstack([ray(SW_DN), rim(-45, 45), ray(TS_UP)[::-1]]),
         "S": np.vstack([ray(TS_DN), rim(225, 315), ray(SW_DN)[::-1]])}


def wedge_of(d):
    """the wedge of the game shown at disk point d"""
    d = np.asarray(d, float)
    u, v = K * d / np.sqrt(1.0 - (d ** 2).sum())
    return ("W" if u - v < 1 else "S") if u + v < 1 else ("N" if u - v < 1 else "E")


LABXY = {"W": (-0.60, -0.26), "S": (0.28, -0.62), "N": (0.28, 0.62), "E": (0.68, 0.19)}
for k, xy in LABXY.items():
    assert wedge_of(xy) == k, (k, xy)
for k, P in WPOLY.items():                                       # each polygon's centroid lies in its own wedge
    assert wedge_of(P.mean(0) * 0.9) == k, k

# ------------------------------------------------------------------ layout (inches), at the reference's native scale
PW, LIM, EXTT, EXTB = A.PW, A.LIM, A.EXTT, 0.04
PH = PW * (2 * LIM + EXTT + EXTB) / (2 * LIM)
FS = A.FS
ML, MR, MT, MB = 0.30, 0.20, 0.12, 0.14
COLGAP = 0.55
FIG_W = ML + 2 * PW + COLGAP + MR
CBS = 0.80                                                       # the colour-bar strip
RH, HH, SEP = 0.36, 0.34, 0.10                                   # key: row height, header row height, gap after the wedges
GAPK = 0.18                                                      # between the colour-bar strip and the key
CW = [0.50, 2.30, 1.40, 1.32, 1.08, 1.32, 1.08]                  # icon, locus, efficiency, M1 binary, M1 dim, M2 binary, M2 dim
KW = sum(CW)
KH = 2 * HH + len(LOCI) * RH + SEP
FIG_H = MB + KH + GAPK + CBS + PH + MT
fig = plt.figure(figsize=(FIG_W, FIG_H))

# ------------------------------------------------------------------ a, b: the disks
TITLE = {"a": r"Memory-1, stochastic, $[0,1]^4$", "b": r"Memory-2, stochastic, $[0,1]^{16}$"}
STROKE = [pe.withStroke(linewidth=3.2, foreground="white")]
ydisk = MB + KH + GAPK + CBS
for i, (letter, DATA) in enumerate((("a", D.M1), ("b", D.M2))):
    x = ML + i * (PW + COLGAP)
    ax = fig.add_axes([x / FIG_W, ydisk / FIG_H, PW / FIG_W, PH / FIG_H])
    for k in WEDGES:
        c = colour(DATA[k]["dim"])
        ax.add_patch(Polygon(WPOLY[k], closed=True, fc=c, ec=c, lw=0.4, zorder=1))
    A.furniture(ax, letter, TITLE[letter])
    ax.set_ylim(-(LIM + EXTB), LIM + EXTT)
    for k, (lx, ly) in LABXY.items():
        ax.text(lx, ly, k, fontsize=FS * 1.45, fontweight="bold", color="0.05", ha="center", va="center", zorder=12,
                path_effects=STROKE)

# ------------------------------------------------------------------ the colour bar: one discrete scale, and grey for empty
CBH, CBW = A.CBH, 3.10
yc = MB + KH + GAPK + 0.52 * CBS
LAB = "dimension of the friendly-rival set"
lab_w = 3.05                                                     # measured below and asserted
SWG, EMP = 0.40, 0.72                                            # gap to the grey swatch, width of its label
tot = lab_w + 0.18 + CBW + SWG + CBH + 0.08 + EMP
xb = (FIG_W - tot) / 2 + lab_w + 0.18
cax = fig.add_axes([xb / FIG_W, (yc - CBH / 2) / FIG_H, CBW / FIG_W, CBH / FIG_H])
cmap = ListedColormap([COL[d] for d in DIMS])
norm = BoundaryNorm(np.arange(len(DIMS) + 1), len(DIMS))
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", drawedges=True,
                  ticks=np.arange(len(DIMS)) + 0.5)
cb.ax.set_xticklabels([str(d) for d in DIMS])
cb.ax.minorticks_off()
cb.ax.tick_params(labelsize=FS * 0.82, length=5, width=1.0)
cb.outline.set_linewidth(1.0)
cb.dividers.set_color("white")
cb.dividers.set_linewidth(0.8)
axL = fig.add_axes([0, 0, 1, 1])
axL.set_xlim(0, FIG_W); axL.set_ylim(0, FIG_H); axL.axis("off")
tl = axL.text(xb - 0.18, yc, LAB, fontsize=FS * 0.95, color="0.12", ha="right", va="center")
xs = xb + CBW + SWG
axL.add_patch(Rectangle((xs, yc - CBH / 2), CBH, CBH, fc=GREY, ec="0.25", lw=0.8))
axL.text(xs + CBH + 0.08, yc, "empty", fontsize=FS * 0.95, color="0.12", ha="left", va="center")

# ------------------------------------------------------------------ the key
X0 = (FIG_W - KW) / 2
XC = X0 + np.concatenate([[0.0], np.cumsum(CW)])                 # column edges
ytop = MB + KH
BFS, HFS = FS * 0.86, FS * 0.80
EFFTXT = {"CC": r"$w_{CC}=1$", "ALT": "alternation", "Z": r"$w_{DD}=0$"}
LOCTXT = {"W": r"$u+v<1$,  $T>S$", "S": r"$u+v<1$,  $T<S$", "N": r"$u+v>1$,  $T>S$", "E": r"$u+v>1$,  $T<S$",
          "rayP": r"switch line,  $T>S$", "rayM": r"switch line,  $T<S$",
          "tsLO": r"$T=S$,  below the switch line", "tsHI": r"$T=S$,  above the switch line",
          "pt": r"$(1,0)$,  where $R=S=T$"}

# header: the two groups, then the columns
yg = ytop - HH
for j0, text in ((3, TITLE["a"]), (5, TITLE["b"])):
    xm = 0.5 * (XC[j0] + XC[j0 + 2])
    axL.text(xm, yg + 0.46 * HH, text, fontsize=FS * 0.84, color="0.12", ha="center", va="center")
    axL.plot([XC[j0] + 0.06, XC[j0 + 2] - 0.06], [yg + 0.06, yg + 0.06], color="0.45", lw=0.8)
yh = ytop - 2 * HH
HEAD = ["locus", "", "efficiency requires", "binary, of 16", "dimension", "binary, of 65536", "dimension"]
for j, h in enumerate(HEAD):
    if not h:
        continue
    if j == 0:
        axL.text(XC[1], yh + HH / 2, h, fontsize=HFS, color="0.25", ha="left", va="center")
    else:
        axL.text(0.5 * (XC[j] + XC[j + 1]), yh + HH / 2, h, fontsize=HFS, color="0.25", ha="center", va="center")


def icon(k, cx, cy, h):
    """the plane of games in miniature, the locus k highlighted: the switch line crimson, T=S dashed"""
    axL.add_patch(Rectangle((cx - h, cy - h), 2 * h, 2 * h, fc="0.95", ec="0.75", lw=0.5, zorder=3))
    corner = {"W": [(-1, 1), (-1, -1)], "N": [(1, 1), (-1, 1)], "E": [(1, -1), (1, 1)], "S": [(-1, -1), (1, -1)]}
    if k in corner:
        pts = [(cx, cy)] + [(cx + h * a, cy + h * b) for a, b in corner[k]]
        axL.add_patch(Polygon(pts, closed=True, fc="0.40", ec="none", zorder=4))
    axL.plot([cx - h, cx + h], [cy + h, cy - h], color="crimson", lw=1.1, zorder=5, solid_capstyle="butt")
    axL.plot([cx - h, cx + h], [cy - h, cy + h], color="0.35", lw=1.0, ls=(0, (2.2, 1.6)), zorder=5)
    end = {"rayP": (-1, 1), "rayM": (1, -1), "tsLO": (-1, -1), "tsHI": (1, 1)}
    if k in end:
        a, b = end[k]
        axL.plot([cx, cx + h * a], [cy, cy + h * b], color="0.05", lw=3.0, zorder=6, solid_capstyle="butt")
    if k == "pt":
        axL.plot([cx], [cy], "o", color="0.05", ms=5.5, zorder=6)


y = yh
for i, k in enumerate(LOCI):
    if k == "rayP":
        y -= SEP
    y -= RH
    cells = [None, None, EFFTXT[D.EFF[k]], str(D.M1[k]["bin"]), D.M1[k]["dim"], str(D.M2[k]["bin"]), D.M2[k]["dim"]]
    for j in range(7):
        x0, w = XC[j], CW[j]
        if j in (4, 6):
            c = colour(cells[j])
            axL.add_patch(Rectangle((x0, y), w, RH, fc=c, ec="white", lw=0.8, zorder=2))
            axL.text(x0 + w / 2, y + RH / 2, "empty" if cells[j] is None else str(cells[j]), fontsize=BFS, ha="center",
                     va="center", color=ink(c), fontweight="normal" if cells[j] is None else "bold", zorder=3)
        else:
            axL.add_patch(Rectangle((x0, y), w, RH, fc="white", ec="0.85", lw=0.5, zorder=2))
            if j >= 2:
                axL.text(x0 + w / 2, y + RH / 2, cells[j], fontsize=BFS, ha="center", va="center", color="0.12", zorder=3)
    icon(k, XC[0] + CW[0] / 2, y + RH / 2, 0.135)
    if k in WEDGES:
        axL.text(XC[1] + 0.02, y + RH / 2, k, fontsize=BFS * 1.05, fontweight="bold", color="0.05", ha="left", va="center")
        axL.text(XC[1] + 0.36, y + RH / 2, LOCTXT[k], fontsize=BFS, color="0.12", ha="left", va="center")
    else:
        axL.text(XC[1] + 0.02, y + RH / 2, LOCTXT[k], fontsize=BFS, color="0.12", ha="left", va="center")
for ya, yb in ((yh - 4 * RH, yh), (y, yh - 4 * RH - SEP)):     # frames around the two blocks
    axL.add_patch(Rectangle((X0, ya), KW, yb - ya, fc="none", ec="0.6", lw=0.7, zorder=7))

# ------------------------------------------------------------------ checks: nothing off the sheet, no text collisions in the key
fig.canvas.draw()
R = fig.canvas.get_renderer()
bbw = tl.get_window_extent(R).width / fig.dpi
assert bbw <= lab_w + 0.02, "the colour-bar label is %.2f in, not <= %.2f" % (bbw, lab_w)
for ax_ in fig.axes:
    for t_ in list(ax_.texts) + (list(ax_.get_xticklabels()) if ax_.axison else []):
        if t_.get_text() and t_.get_visible():
            bb = t_.get_window_extent(R)
            assert 0 < bb.x0 and bb.x1 < FIG_W * fig.dpi and 0 < bb.y0 and bb.y1 < FIG_H * fig.dpi, "off the sheet: %r" % t_.get_text()
keytexts = [t_ for t_ in axL.texts if t_.get_text()]
for i_ in range(len(keytexts)):
    for j_ in range(i_ + 1, len(keytexts)):
        a_, b_ = keytexts[i_].get_window_extent(R), keytexts[j_].get_window_extent(R)
        assert not a_.overlaps(b_), "texts overlap: %r / %r" % (keytexts[i_].get_text(), keytexts[j_].get_text())
# every key text inside its cell
for t_ in keytexts:
    bb = t_.get_window_extent(R)
    x0_, x1_ = bb.x0 / fig.dpi, bb.x1 / fig.dpi
    if bb.y1 / fig.dpi < ytop:                                   # a key text
        assert not [e for e in XC[1:-1] if x0_ + 0.01 < e < x1_ - 0.01 and t_.get_text() not in (TITLE["a"], TITLE["b"])], \
            "key text crosses a column edge: %r" % t_.get_text()

fig.savefig(os.path.join(HERE, "Fig1.pdf"))
fig.savefig(os.path.join(HERE, "Fig1.png"), dpi=200)
print("wrote Fig1.pdf and Fig1.png (%.2f x %.2f in); dimension colours: %s"
      % (FIG_W, FIG_H, ", ".join("%d %s" % (d, matplotlib.colors.to_hex(COL[d])) for d in DIMS)))
