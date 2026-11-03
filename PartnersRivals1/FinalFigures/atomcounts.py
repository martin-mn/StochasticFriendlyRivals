#!/usr/bin/env python3
"""
Shared machinery of Figure 2 and SI Figure 5: the number of strategies in each of the seven atoms of the
three properties -- efficient, nash, competitive -- at each game, on the disk of games (K = 4), seven
panels in the atom order 000, 100, 010, 001 (row 1) and 110, 011, 111 (row 2); the eighth cell is empty
(101 has no member at any game; until 2026-09-18 a note in the cell said so, and the ramp colour bars were
labelled "equal-area scale").

    draw(FACES, DA, ATOMS, ATOM_TOTAL, out)              the seven-atom page (two rows of four)
    sheet(nrow, ncol, labw) / panel(...) / rowlabel(...) / finish(...)   the pieces, for sheets with other rows
                                                          (fig2.py: memory one and two, four rows; fig3.py: the sets)

FACES are the faces of an exact eps -> 0 arrangement as polygons in disk coordinates, DA their areas,
ATOMS (nfaces x 8) the number of strategies of each atom on each face (code = 4 efficient + 2 nash +
competitive), ATOM_TOTAL (8) the number of different strategies that belong to the atom at some game.  Each
panel has its own colour scale over the values the atom takes: discrete bands when the atom takes at most
twelve values (one band if it takes one), otherwise an equal-area ramp; faces where the atom has no member
are grey, marked 0 on the bar.  Under each disk: the number of different strategies over the disk.  The
script asserts that no text leaves the sheet and that no two heads overlap.  Style and furniture are those
of ../Figures/fig2.py; fig2.py (memory two) and sifig5.py (memory one) are the two callers.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "Figures"))   # disklib

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap

import disklib
from disklib import NAMES, arc_text

K = 4.0
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 15, "axes.linewidth": 1.2})


def phi(p):
    return disklib.phi(p, K)


def grp(x):
    return "%d" % round(x)


# --------------------------------------------------------------- colour bars
NB = 256


def ladder(lo, hi):
    t, d = [], int(np.floor(np.log10(lo)))
    while 10.0 ** d <= hi * 10:
        for m in (1, 2, 5):
            x = m * 10.0 ** d
            if lo * 1.25 < x < hi / 1.25:
                t.append(x)
        d += 1
    return [lo] + t + [hi]


def thin(ticks, norm, nbands, mingap=0.052):
    pos = [(float(norm(t)) / float(nbands), t) for t in ticks]
    keep = [pos[0], pos[-1]]
    for q, t in sorted(pos[1:-1], key=lambda z: -z[0]):
        if all(abs(q - k[0]) >= mingap for k in keep):
            keep.append((q, t))
    return [t for _, t in sorted(keep)]


def equalised(vals, weights, cmap="magma", name=""):
    """Figure 2a's scale: the colormap spread over the field's own area-weighted distribution."""
    v = np.asarray(vals, dtype=float); w = np.asarray(weights, dtype=float)
    o = np.argsort(v, kind="stable")
    cum = np.cumsum(w[o]); cum /= cum[-1]
    bnd = np.interp(np.linspace(0.0, 1.0, NB + 1), cum, v[o])
    bnd[0], bnd[-1] = v.min(), v.max()
    bnd = np.unique(bnd); nb = len(bnd) - 1
    norm = BoundaryNorm(bnd, nb)
    cm = ListedColormap(plt.get_cmap(cmap)(np.linspace(0.0, 1.0, nb)))
    ticks = thin(ladder(float(v.min()), float(v.max())), norm, nb)
    print("bar %s: %d equal-area bands over %g to %g" % (name, nb, v.min(), v.max()))
    return cm, norm, ticks


# ------------------------------------------------------------------- layout
MODE = "atoms"
PORTRAIT = False
TXT = 2.0 if PORTRAIT else 1.45  # text scale relative to fig2.py's rule (many columns print smaller)
PW = 5.0 if MODE == "classes" else 4.3
LIM = 1.175
EXTT, EXTB = 0.22, 0.34
PH = PW * (2 * LIM + EXTT + EXTB) / (2 * LIM)
CBW, CBGAP, CBLAB = 0.22, 0.26, 1.30
COLW = PW + CBGAP + CBW + CBLAB
COLGAP, ROWGAP = 0.30, 0.25
ML, MR, MB, MT = 0.30, 0.20, 0.25, 0.25
NCOL, NROW = 4, 2
FIG_W = ML + NCOL * COLW + (NCOL - 1) * COLGAP + MR
FIG_H = MB + NROW * PH + (NROW - 1) * ROWGAP + MT
DPU = (2 * LIM / PW) / 72.0
SC = PW / 13.46
NM_FS = disklib.NM_FS * SC * (1.40 if PORTRAIT else 1.25)
FS = max(9.0, 13.0 * SC * 1.9) * TXT
CBH = 0.80 * PW
RB = 2000.0

th = np.linspace(0, 2 * np.pi, 4001)
GUIDE = (1, 2, 4, 10)
CONICS = [(s, col, dash)
          for (a, b, c), col, dash, arms in (((1, 1, 1), "crimson", (0, (7, 4)), (-1, 1)),
                                             ((1, -1, 1), "white", (0, (7, 4)), (-1, 1)),
                                             ((1, 1, 0), "#12d7e8", (0, (8, 3, 1.5, 3)), (1,)))
          for s in disklib.conic(a, b, c, phi, RB, arms)]
PEL = [pe.withStroke(linewidth=2.4, foreground="white")]
GUIDE_C, GUIDE_LW, GUIDE_A = "white", 0.6, 0.28

ATOM_ORDER = ["000", "100", "010", "001", "110", "011", "111"]
ATOM_NAME = {"000": "none of the three", "100": "efficient only", "010": "stable only", "001": "competitive only",
             "110": "efficient and stable only", "011": "stable and competitive only", "111": "efficient, stable and competitive"}
GREY = "0.84"


def cellxy(row, col):
    return (ML + col * (COLW + COLGAP), MB + (NROW - 1 - row) * (PH + ROWGAP))


def atom_scale(vals, weights):
    """ramp or bands over the positive values; zero is the under-colour GREY"""
    pos = vals > 0
    lv = np.unique(vals[pos])
    if len(lv) == 1:                                                  # one band (memory one: 111 is always 2)
        cm = ListedColormap([plt.get_cmap("magma")(0.75)])
        norm = BoundaryNorm([lv[0] - 0.5, lv[0] + 0.5], 1)
        ticks, lab, disc = [float(lv[0])], [grp(lv[0])], True
    elif len(lv) <= 12:
        cm, norm, ticks, lab, rec = disklib.discrete_bands(vals[pos], weights[pos], cmap="magma", vmin=lv.min(),
                                                            vmax=lv.max(), group=grp, by="rank", share=False)
        ticks, lab, disc = [float(x) for x in lv], [grp(x) for x in lv], True
        cm = ListedColormap(cm.colors)
    else:
        cm, norm, ticks = equalised(vals[pos], weights[pos], name="atom")
        lab, disc = [grp(t) for t in ticks], False
    cm.set_under(GREY)
    return cm, norm, ticks, lab, disc


def sheet(nrow, ncol=NCOL, labw=0.0):
    """a sheet of nrow x ncol cells of Figure 2's size, with a column of width labw (inches) on the left for row
    labels; returns (fig, FIG_W, FIG_H, cellxy) with cellxy(row, col) the lower-left corner of a disk in inches"""
    W = ML + labw + ncol * COLW + (ncol - 1) * COLGAP + MR
    H = MB + nrow * PH + (nrow - 1) * ROWGAP + MT
    fig = plt.figure(figsize=(W, H))

    def xy(row, col):
        return (ML + labw + col * (COLW + COLGAP), MB + (nrow - 1 - row) * (PH + ROWGAP))
    return fig, W, H, xy


def panel(fig, W, H, x, y, FACES, DA, vals, letter, head, sub, barlabel="Number of strategies"):
    """one disk of Figure 2 at (x, y) inches on a sheet of W x H inches: the faces FACES coloured by vals with the
    scale of atom_scale (grey where zero), the furniture, the head and the sub-line, and the colour bar to its right"""
    vals = np.asarray(vals, float)
    ax = fig.add_axes([x / W, y / H, PW / W, PH / H])
    zero = vals == 0
    cm, norm, ticks, lab, disc = atom_scale(vals, DA)
    c = cm(norm(vals))
    ax.add_collection(PolyCollection(FACES, facecolors=c, edgecolors=c, linewidths=0.30, antialiased=True, zorder=1))
    ax.plot([-1, 1], [0, 0], color=GUIDE_C, lw=GUIDE_LW, zorder=5, alpha=GUIDE_A)
    ax.plot([0, 0], [-1, 1], color=GUIDE_C, lw=GUIDE_LW, zorder=5, alpha=GUIDE_A)
    for s_, col, dash in CONICS:
        ax.plot(s_[:, 0], s_[:, 1], color=col, lw=1.5, zorder=7, linestyle=dash)
    for r in GUIDE:
        rr = float(np.hypot(*phi([[r, 0.0]])[0]))
        ax.plot(rr * np.cos(th), rr * np.sin(th), "-", color=GUIDE_C, lw=GUIDE_LW, zorder=5, alpha=GUIDE_A)
        for sg in (-1, 1):
            ax.plot([sg * rr, sg * rr], [-0.016, 0.016], color="0.2", lw=1.4, zorder=8, path_effects=PEL)
            ax.plot([-0.016, 0.016], [sg * rr, sg * rr], color="0.2", lw=1.4, zorder=8, path_effects=PEL)
            lb = ("%g" % r) if sg > 0 else "$-$%g" % r
            ax.text(sg * rr, -0.048, lb, fontsize=FS * 0.7, color="0.2", ha="center", va="top", zorder=9, path_effects=PEL)
            ax.text(-0.048, sg * rr, lb, fontsize=FS * 0.7, color="0.2", ha="right", va="center", zorder=9, path_effects=PEL)
    ax.text(1.055, 0.0, "$u$", fontsize=FS * 1.4, ha="left", va="center")
    ax.text(0.0, 1.055, "$v$", fontsize=FS * 1.4, ha="center", va="bottom")
    ax.text(-LIM + 0.01, LIM + 0.055, letter, fontsize=FS * 1.7, fontweight="bold", ha="left", va="bottom")
    ax.text(0.0, LIM + 0.055, head, fontsize=FS * 1.05, color="0.12", ha="center", va="bottom")
    for a0, nm in NAMES:
        _, _, ro = arc_text(ax, nm, a0, 1.0, NM_FS, DPU, disklib.NM_GAP, disklib.NM_TRK)
        assert ro < LIM, "%s reaches r = %.4f" % (nm, ro)
    ax.text(0.0, -LIM - 0.06, sub, fontsize=FS * 0.82, color="0.12", ha="center", va="top", zorder=12, linespacing=1.35)
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-(LIM + EXTB), LIM + EXTT)
    ax.set_aspect("equal")
    ax.axis("off")
    # the colour bar
    xb = x + PW + CBGAP
    yc = y + PH * (LIM + EXTB) / (2 * LIM + EXTT + EXTB)
    cax = fig.add_axes([xb / W, (yc - 0.5 * CBH) / H, CBW / W, CBH / H])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cm), cax=cax, ticks=ticks, spacing="uniform", drawedges=disc,
                      extend="min" if zero.any() else "neither", extendfrac=0.06)
    cb.ax.set_yticklabels(lab)
    cb.ax.minorticks_off()
    for lb in (cb.ax.get_yticklabels()[0], cb.ax.get_yticklabels()[-1]):
        lb.set_fontweight("bold")
    cb.set_label(barlabel, fontsize=FS * 0.95, labelpad=10)
    cb.ax.tick_params(labelsize=FS * (0.72 if disc else 0.82), length=5, width=1.0)
    if disc:
        cb.outline.set_linewidth(1.0)
        cb.dividers.set_color("white")
        cb.dividers.set_linewidth(0.6)
    if zero.any():
        cb.ax.text(1.45, -0.045, "0", transform=cb.ax.transAxes, fontsize=FS * 0.72, ha="left", va="center", color="0.25")
    return dict(min=vals[~zero].min() if (~zero).any() else np.nan, max=vals.max(), levels=len(np.unique(vals[~zero])),
                zshare=DA[zero].sum() / DA.sum())


def rowlabel(fig, W, H, x, y, text):
    """a rotated label centred on the disk whose lower-left corner is at (x, y) inches, printed 0.45 in to its left"""
    fig.text((x - 0.45) / W, (y + PH * (LIM + EXTB) / (2 * LIM + EXTT + EXTB)) / H, text, rotation=90,
             ha="center", va="center", fontsize=FS * 1.25, color="0.12")


def finish(fig, W, H, out):
    """check that no text leaves the sheet and no two heads overlap, then write out (.pdf) and its .png"""
    fig.canvas.draw()
    R = fig.canvas.get_renderer()
    for t_ in fig.texts:
        bb = t_.get_window_extent(R)
        assert 0 < bb.x0 and bb.x1 < W * fig.dpi and 0 < bb.y0 and bb.y1 < H * fig.dpi, "text runs off the sheet: %r" % t_.get_text()
    for ax_ in fig.axes:
        for t_ in list(ax_.texts) + (ax_.get_yticklabels() if ax_.axison else []):
            if t_.get_text() and t_.get_visible():
                bb = t_.get_window_extent(R)
                assert 0 < bb.x0 and bb.x1 < W * fig.dpi and 0 < bb.y0 and bb.y1 < H * fig.dpi, \
                    "text runs off the sheet: %r" % t_.get_text()
        heads = [t_ for t_ in ax_.texts if t_.get_text() and t_.get_va() == "bottom" and t_.get_ha() in ("center", "left")]
        for i_ in range(len(heads)):
            for j_ in range(i_ + 1, len(heads)):
                a_, b_ = heads[i_].get_window_extent(R), heads[j_].get_window_extent(R)
                assert not (a_.overlaps(b_)), "texts overlap: %r / %r" % (heads[i_].get_text(), heads[j_].get_text())
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=150)
    plt.close(fig)
    print("%s: %.2f x %.2f in, pdf %.1f MB" % (os.path.basename(out), W, H, os.path.getsize(out) / 1e6))


def draw(FACES, DA, ATOMS, ATOM_TOTAL, out):
    """draw the seven atom panels for the faces FACES (disk polygons) with areas DA and atom counts ATOMS
    (nfaces x 8), ATOM_TOTAL[c] = number of different strategies of atom c over the disk; write out (.pdf) and
    its .png.  Two rows of four; the eighth cell stays empty (101 has no member at any game)."""
    ATOMS = np.asarray(ATOMS, float)
    DA = np.asarray(DA, float)
    fig, W, H, xy = sheet(NROW, NCOL)
    for i, code in enumerate(ATOM_ORDER):
        row, col = divmod(i, NCOL)
        c = int(code, 2)
        x, y = xy(row, col)
        st = panel(fig, W, H, x, y, FACES, DA, ATOMS[:, c], "abcdefg"[i], r"$\mathbf{%s}$   %s" % (code, ATOM_NAME[code]),
                   "%d different strategies over the disk" % ATOM_TOTAL[c])      # the range per game is on the colour bar
        print("atom %s: %d strategies on an open set; %d-%d per game where present; %d levels; empty on %.1f%% of the disk area"
              % (code, ATOM_TOTAL[c], st["min"], st["max"], st["levels"], 100 * st["zshare"]))
    finish(fig, W, H, out)
