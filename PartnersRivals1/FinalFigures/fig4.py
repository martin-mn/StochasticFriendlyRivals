#!/usr/bin/env python3
"""
Figure 4 (Martin's design of 2026-09-18, final in his plan of 2026-09-19): one summary of the four Wright-Fisher runs.

    python3 fig4.py             ->  Fig4.pdf / Fig4.png  (here; drafted as figwf.py -> FigWF.*)

Twelve disks, one per (run, quantity):

    row 1   memory one, N = 100     (DiskM1WF run F2:  beta = 3,   mu = 1e-4, eps = 1e-4, 10 x 1e8 generations)
    row 2   memory one, N = 1000    (DiskM1WF run F1:  beta = 100, mu = 1e-2, eps = 1e-4, 10 x 1e7)
    row 3   memory two, N = 100     (DiskM2WF run F2:  the same as row 1 over the 65536 memory-two strategies)
    row 4   memory two, N = 1000    (DiskM2WF run F1:  the same as row 2)

    column 1   efficiency: the mean per-round payoff of the population over the second half of the run, divided
               by the joint optimum Emax = max(R, (S + T)/2) of the game (the "ef" of the runs; the realised
               efficiency of the Methods); viridis from 0 to 1, black where it is negative
    column 2   the most abundant atom: the atom whose strategies hold the largest share of the population
               (panel h of the pages of SI Figures 2 and 3)
    column 3   the most enriched atom: the atom with the largest share divided by its share of the strategy
               space (panel h of the retired enrichment pages in other/: Fig4_old, SIFig3_old, M1WF_enrichment_*)

The 512 games are the same Voronoi cells of the Fermat sunflower in all twelve panels; frame, furniture, atom
colours and swatch legend are those of Figures 3 and 4 (atomshares.py, imported).  Data: the memory-two runs
through DiskM2WF/Opt (figF4.collect for the shares, wfload for the efficiency), the memory-one runs through
DiskM1WF/Opt/m1load (its packs Data/mw1_e4.tsv, mw2_e4.tsv); all read-only.  Both classifications are the exact
eps -> 0 atoms of this paper.  Panel letters a-l row by row.

The sheet is drawn by draw(rows, stem), one row per (label, E, wa, we); sifig9.py draws the error-rate series (SI Figures 9 and 10)
with the same function, so the sheets share their frame; sifig8.py (the cube, 24 disks) reuses the colours and furniture.
The refactor of 2026-09-22 that exposed draw() left Fig4.png pixel-identical.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
import wfdata                                                      # noqa: E402  the four runs, read one way (also puts DiskM2WF/Opt on sys.path)
import atomshares as A                                             # noqa: E402  frame, furniture, colours

ORDER, WINCOL, GREY = A.ORDER, A.WINCOL, A.GREY
ECMAP = plt.get_cmap("viridis").copy()
ECMAP.set_under("black")                                           # a negative mean payoff
ENORM = Normalize(vmin=0.0, vmax=1.0)

ROWS = [(wfdata.LABEL[k],) + k for k in wfdata.ROWS if k[1] == "F1"]      # Figure 4: the large population (Martin, 2026-09-23)
SIROWS = [(wfdata.LABEL[k],) + k for k in wfdata.ROWS if k[1] == "F2"]    # SI Figure 4: the small population
COLS = ["efficiency", "most abundant atom", "most enriched atom"]


def winners(S, N, E, D, name):
    """per cell: the most abundant and the most enriched atom (indices into ORDER), with a printed summary"""
    share = np.where(N > 0, S, -np.inf)
    enr = np.where(N > 0, S / np.maximum(N, 1), -np.inf)           # times the constant |space|, which does not move the argmax
    wa, we = share.argmax(1), enr.argmax(1)
    print("%s (N = %d, beta = %d, mu = %s): efficiency %.3f .. %.3f, median %.3f, negative at %d games, below 0.9 at %d"
          % (name, D["N"], D["beta"], D["muplain"], E.min(), E.max(), np.median(E), (E < 0).sum(), (E < 0.9).sum()))
    print("   most abundant: " + ", ".join("%s %d" % (ORDER[i], (wa == i).sum()) for i in range(7) if (wa == i).any()))
    print("   most enriched: " + ", ".join("%s %d" % (ORDER[i], (we == i).sum()) for i in range(7) if (we == i).any()))
    return wa, we


def run(space, which):
    """per cell (in ipt order 0..511): efficiency, the most abundant atom, the most enriched atom (indices into ORDER)"""
    S, N, E, UV, D = wfdata.load(space, which)
    wa, we = winners(S, N, E, D, "%s %s" % (space, which))
    return E, wa, we, UV


# ------------------------------------------------------------------- layout (inches; the panel of Figures 3 and 4)
PW, LIM, EXTT = A.PW, A.LIM, A.EXTT
EXTB = 0.08                                                        # no swatch row under the disks here: the legends sit at the bottom
PH = PW * (2 * LIM + EXTT + EXTB) / (2 * LIM)
ML, MR, MT, MB = 0.95, 0.15, 0.15, 0.20
COLGAP, ROWGAP = 0.45, 0.15
LEG = 0.80                                                         # the legend strip at the bottom


def draw(rows, stem):
    """the sheet: rows = [(label, E, wa, we)] top to bottom, three panels each; writes HERE/<stem>.pdf and .png"""
    NR = len(rows)
    FIG_W = ML + 3 * PW + 2 * COLGAP + MR
    FIG_H = MB + LEG + NR * PH + (NR - 1) * ROWGAP + MT
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    def cell(row, col):
        x = ML + col * (PW + COLGAP)
        y = MB + LEG + (NR - 1 - row) * (PH + ROWGAP)
        return fig.add_axes([x / FIG_W, y / FIG_H, PW / FIG_W, PH / FIG_H]), x, y

    pol, keep = wfdata.cells()
    LET = "abcdefghijklmnopqrstuvwxyz"[:3 * NR]
    letters = iter(LET)
    for r, (label, E, wa, we) in enumerate(rows):
        for c, title in enumerate(COLS):
            ax, x, y = cell(r, c)
            if c == 0:
                fc = [ECMAP(ENORM(E[i])) if E[i] >= 0 else ECMAP.get_under() for i in keep]
            else:
                win = wa if c == 1 else we
                fc = [WINCOL[ORDER[win[i]]] for i in keep]
            ax.add_collection(PolyCollection(pol, facecolors=fc, edgecolors=fc, linewidths=0.3, antialiased=True, zorder=1))
            A.furniture(ax, next(letters), title if r == 0 else "")
            ax.set_ylim(-(LIM + EXTB), LIM + EXTT)                 # furniture set Figures 3-4's bottom extension; trim it
            if c == 0:
                fig.text((ML - 0.62) / FIG_W, (y + PH * (EXTB + LIM) / (2 * LIM + EXTT + EXTB)) / FIG_H, label, rotation=90,
                         ha="center", va="center", fontsize=A.FS * 1.25, color="0.12")

    # --------------------------------------------------------------- the legends: the efficiency bar under column 1, the atoms under columns 2-3
    CBH = A.CBH
    YLEG = MB + 0.5 * LEG                                          # centre line of the legend strip
    cax = fig.add_axes([(ML + 0.18 * PW) / FIG_W, (YLEG - 0.5 * CBH) / FIG_H, 0.64 * PW / FIG_W, CBH / FIG_H])
    cb = fig.colorbar(ScalarMappable(norm=ENORM, cmap=ECMAP), cax=cax, orientation="horizontal", extend="min", extendfrac=0.05)
    cb.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]); cb.ax.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"])
    cb.ax.minorticks_off(); cb.ax.tick_params(labelsize=A.FS * 0.82, length=5, width=1.0); cb.outline.set_linewidth(1.0)
    cax.text(0.5, 1.45, "efficiency (%s); black: below 0" % ", ".join(LET[0::3]), transform=cax.transAxes, ha="center", va="bottom",
             fontsize=A.FS * 0.95, color="0.12")
    axL = fig.add_axes([0, 0, 1, 1]); axL.set_xlim(0, FIG_W); axL.set_ylim(0, FIG_H); axL.set_aspect("equal"); axL.axis("off")
    SW, LFS = CBH, A.LFS
    LW = A.LW * PW / (2 * LIM)                                     # the width of a bold code, inches
    PITCH = SW + 0.03 + LW + 0.06
    xc = ML + PW + COLGAP + 0.5 * (2 * PW + COLGAP)                # centre of columns 2-3
    x0 = xc - 0.5 * (6 * PITCH + SW + 0.03 + LW)
    for i, code in enumerate(ORDER):
        xi = x0 + i * PITCH
        axL.add_patch(Rectangle((xi, YLEG - 0.5 * SW), SW, SW, fc=WINCOL[code], ec="0.25", lw=0.6, zorder=12))
        axL.text(xi + SW + 0.03, YLEG, code, fontsize=LFS, fontweight="bold", color="0.12", ha="left", va="center", zorder=12)
    axL.text(xc, YLEG + 0.5 * SW + 0.12, "the atom (%s)" % ", ".join(l for i, l in enumerate(LET) if i % 3), ha="center", va="bottom",
             fontsize=A.FS * 0.95, color="0.12")

    # --------------------------------------------------------------- checks: nothing off the sheet
    fig.canvas.draw()
    R = fig.canvas.get_renderer()
    for ax_ in fig.axes:
        for t_ in list(ax_.texts) + (ax_.get_xticklabels() if ax_.axison else []):
            if t_.get_text() and t_.get_visible():
                bb = t_.get_window_extent(R)
                assert 0 < bb.x0 and bb.x1 < FIG_W * fig.dpi and 0 < bb.y0 and bb.y1 < FIG_H * fig.dpi, "text runs off the sheet: %r" % t_.get_text()
    for t_ in fig.texts:
        bb = t_.get_window_extent(R)
        assert 0 < bb.x0 and bb.x1 < FIG_W * fig.dpi and 0 < bb.y0 and bb.y1 < FIG_H * fig.dpi, "text runs off the sheet: %r" % t_.get_text()

    out = os.path.join(HERE, stem + ".pdf")
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=110)
    plt.close(fig)
    print("wrote FinalFigures/%s.pdf and .png (%.2f x %.2f in)" % (stem, FIG_W, FIG_H))


if __name__ == "__main__":
    draw([(label,) + run(space, which)[:3] for label, space, which in ROWS], "Fig4")
    draw([(label,) + run(space, which)[:3] for label, space, which in SIROWS], os.path.join("SI", "SIFig4"))
