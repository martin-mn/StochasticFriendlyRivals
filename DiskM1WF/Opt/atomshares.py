#!/usr/bin/env python3
"""
atomshares.py -- the machinery of this project's four pages: the seven atoms of the three properties
(efficient, nash, competitive) over the disk of games in a Wright-Fisher run of the 16 BINARY
MEMORY-ONE strategies.  The memory-one counterpart of PartnersRivals/FinalFigures/atomshares.py,
which draws the same pages for the 65536 binary memory-two strategies, and built FROM it by
Opt/mkshares.py as an asserted list of eleven substitutions -- so the layout, the colours, the panel
letters, the legend of h and the checks are that module's, byte for byte, and this docstring plus
seven small edits are the whole of the difference.

    page F1   f1.py    run F1:  N = 1000,  beta = 100,  mu = 1e-2,  eps = 1e-4   share of the population
    page F2   f2.py    run F1                                                    enrichment
    page F3   f3.py    run F2:  N = 100,   beta = 3,    mu = 1e-4,  eps = 1e-4   share of the population
    page F4   f4.py    run F2                                                    enrichment

The PAGES are F1..F4 (renamed 2026-09-18 from M1Fig3, M1Fig4, M1SIFig2, M1SIFig3, which are the
PartnersRivals figures they correspond to).  The two RUNS also happen to be called F1 and F2, after
DiskM2WF; run F1 draws pages F1 and F2, run F2 draws pages F3 and F4.  Figures/parameters.tex
tabulates the lot.

The runs are DiskM1WF's own: dw1's process and dw1's run length at dw1's 512 sunflower games, with
the strategy space cut from 65536 to 16 (Cannon/mw, Opt/mksf.py).  The atom of each strategy at each
game is exact in the limit eps -> 0 (Opt/m1disk.py, through FinalFigures/m1atoms.py), and that module
also verified that the three properties of a memory-one strategy against the other fifteen are the
same as against all 65536 memory-two strategies -- so these panels are about memory one, not about a
truncation of memory two.

    row 1:   a  000     b  100     c  010     d  001
    row 2:   e  110     f  011     g  111     h  the most abundant atom at each game
    bottom:  the colour bar of a-g under columns 1-2, level with the swatch legend of h; no notes on the sheet

mode="share" is the share of the population held by the atom (M1 Figures 3 and SI 2; viridis,
linear); mode="enrichment" is that share divided by the atom's share of the 16 strategies, log10 on
RdBu with white at 1 (red depleted, blue enriched).  THE ENRICHMENT RANGE IS NOT THE MEMORY-TWO
MODULE'S.  An atom of n of the 16 strategies has enrichment at most 16/n -- 16 for a one-member atom,
and nothing above 1.21 decades ever -- so the memory-two range 10^-4 ... 10^4 would leave every
enriched cell within a few percent of white.  The map runs from 10^-4, the floor both runs reach,
to 16, the ceiling, with white at 1 and each half taking half the bar; both memory-one pages use it.

Cells where the atom has no strategy at that game are grey, as in SI Figure 5 (the memory-one
Figure 2).  Panel h colours each game by the atom with the largest share (largest enrichment in the
enrichment mode), with a swatch legend under the disk.  Frame and furniture are Figure 2's; the
projection is the same disk (K = 4).  Reads DiskM2WF/Opt and PartnersRivals/Figures read-only.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import FuncNorm, Normalize, TwoSlopeNorm
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath

HERE = "/Users/martin/Documents/CL1/DiskM1WF/Figures"           # where the pages are written
MOD = os.path.dirname(os.path.abspath(__file__))               # and where this module and m1load live
FIGDIR = "/Users/martin/Documents/CL1/PartnersRivals/Figures"  # disklib: the projection and the furniture
OPT = "/Users/martin/Documents/CL1/DiskM2WF/Opt"               # vor: the 512 sunflower cells
os.makedirs(os.path.join(HERE, "other"), exist_ok=True)
sys.path.insert(0, MOD)
sys.path.insert(0, FIGDIR)
sys.path.insert(0, OPT)
import disklib                                                     # noqa: E402
from disklib import NAMES, arc_text                                # noqa: E402
import vor                                                         # noqa: E402
# m1load carries collect(), DATASETS and K under the same names as DiskM2WF's figF4 and
# finalpanels, so every use of them below is the memory-two module's line unchanged.
import m1load as figF4                                             # noqa: E402
import m1load as finalpanels                                       # noqa: E402

matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 15, "axes.linewidth": 1.2})
K = 4.0
assert finalpanels.K == K


def phi(p):
    return disklib.phi(p, K)


ORDER = ["000", "100", "010", "001", "110", "011", "111"]
NAME = {"000": "none of the three", "100": "efficient only", "010": "nash only", "001": "competitive only",
        "110": "efficient and nash only", "011": "nash and competitive only", "111": "efficient, nash and competitive"}
GREY = "0.84"                                    # no strategy in the atom at that game (Figure 2's grey)
WINCOL = {"000": "0.36", "100": "#8fb8e8", "010": "#dba622", "001": "#f0a48c",     # panel h: the most abundant atom
          "110": "#2f5fb3", "011": "#c43d3b", "111": "#7a4fa8"}
CMAP = plt.get_cmap("viridis").copy()
CMAP.set_under("#270031")                        # logit scale: present, but holding less than FLOOR
FLOOR, CEIL = 1e-3, 0.99
# Enrichment: the log10 range of the diverging map, white at 1.  NOT the memory-two
# module's -4 .. 4: an atom of n of the 16 strategies has enrichment at most 16/n, so
# 16 (1.204 decades) is the ceiling and the upper half of that range is unreachable.
# 10^-4 is the floor both runs reach -- a strategy eliminated inside a generation is in
# the population only because mutation put it there.  Each half takes half the bar.
LO, HI = -4.0, np.log10(16.0)
ECMAP = plt.get_cmap("RdBu")
NSTRAT = 16.0                                    # the 16 binary memory-one strategies

# ------------------------------------------------------------------- layout (inches; Figure 2's rules)
PW = 4.3
LIM = 1.175
EXTT, EXTB = 0.22, 0.34
PH = PW * (2 * LIM + EXTT + EXTB) / (2 * LIM)
COLGAP, ROWGAP = 0.55, 0.25
ML, MR, MB, MT = 0.30, 0.20, 0.25, 0.25
NCOL, NROW = 4, 2
FOOT = 0.0                                       # no footer strip: the colour bar sits in the row-2 cells' own bottom margin
FIG_W = ML + NCOL * PW + (NCOL - 1) * COLGAP + MR
FIG_H = MB + FOOT + NROW * PH + (NROW - 1) * ROWGAP + MT
DPU = (2 * LIM / PW) / 72.0
SC = PW / 13.46
NM_FS = disklib.NM_FS * SC * 1.25
FS = max(9.0, 13.0 * SC * 1.9) * 1.45
RB = 2000.0

th = np.linspace(0, 2 * np.pi, 4001)
GUIDE = (1, 2, 4, 10)
CONICS = [(s, col, dash)
          for (a, b, c), col, dash, arms in (((1, 1, 1), "crimson", (0, (7, 4)), (-1, 1)),
                                             ((1, -1, 1), "white", (0, (7, 4)), (-1, 1)),
                                             ((1, 1, 0), "#12d7e8", (0, (8, 3, 1.5, 3)), (1,)))
          for s in disklib.conic(a, b, c, phi, RB, arms)]
PEL = [pe.withStroke(linewidth=2.4, foreground="white")]
HALO = [pe.withStroke(linewidth=2.7, foreground="0.35")]   # dark rim on the white T = S dashes: visible on the light cells too
GUIDE_C, GUIDE_LW, GUIDE_A = "white", 0.6, 0.28
CBH = 0.22                                       # the colour bar's height, inches
SW = CBH * 2 * LIM / PW                          # panel h key: swatch side = the bar's height, in data units
LFS = FS * 0.82                                  # its codes at the size of the bar's tick labels
LW = TextPath((0, 0), "000", size=LFS, prop=FontProperties(family="DejaVu Sans", weight="bold")).get_extents().width * DPU
PITCH = SW + 0.02 + LW + 0.035                   # swatch, gap, code, gap to the next swatch
X0 = -0.5 * (6 * PITCH + SW + 0.02 + LW)         # the row of seven centred under the disk
YL = -LIM - 0.04 - 0.5 * SW - 0.01               # its centre line; the colour bar sits level with it


def cellxy(row, col):
    return (ML + col * (PW + COLGAP), MB + FOOT + (NROW - 1 - row) * (PH + ROWGAP))


def cell(fig, row, col):
    x, y = cellxy(row, col)
    return fig.add_axes([x / FIG_W, y / FIG_H, PW / FIG_W, PH / FIG_H])


def furniture(ax, letter, title):
    ax.plot([-1, 1], [0, 0], color=GUIDE_C, lw=GUIDE_LW, zorder=5, alpha=GUIDE_A)
    ax.plot([0, 0], [-1, 1], color=GUIDE_C, lw=GUIDE_LW, zorder=5, alpha=GUIDE_A)
    for s, col, dash in CONICS:
        ax.plot(s[:, 0], s[:, 1], color=col, lw=1.5, zorder=7, linestyle=dash, **({"path_effects": HALO} if col == "white" else {}))
    for r in GUIDE:
        rr = float(np.hypot(*phi([[r, 0.0]])[0]))
        ax.plot(rr * np.cos(th), rr * np.sin(th), "-", color=GUIDE_C, lw=GUIDE_LW, zorder=5, alpha=GUIDE_A)
        for sg in (-1, 1):
            ax.plot([sg * rr, sg * rr], [-0.016, 0.016], color="0.2", lw=1.4, zorder=8, path_effects=PEL)
            ax.plot([-0.016, 0.016], [sg * rr, sg * rr], color="0.2", lw=1.4, zorder=8, path_effects=PEL)
            lab = ("%g" % r) if sg > 0 else "$-$%g" % r
            ax.text(sg * rr, -0.048, lab, fontsize=FS * 0.7, color="0.2", ha="center", va="top", zorder=9, path_effects=PEL)
            ax.text(-0.048, sg * rr, lab, fontsize=FS * 0.7, color="0.2", ha="right", va="center", zorder=9, path_effects=PEL)
    ax.text(1.055, 0.0, "$u$", fontsize=FS * 1.4, ha="left", va="center")
    ax.text(0.0, 1.055, "$v$", fontsize=FS * 1.4, ha="center", va="bottom")
    ax.text(-LIM + 0.01, LIM + 0.055, letter, fontsize=FS * 1.7, fontweight="bold", ha="left", va="bottom")
    ax.text(0.0, LIM + 0.055, title, fontsize=FS * 1.05, color="0.12", ha="center", va="bottom")
    for a0, nm in NAMES:
        _, _, ro = arc_text(ax, nm, a0, 1.0, NM_FS, DPU, disklib.NM_GAP, disklib.NM_TRK)
        assert ro < LIM, "%s reaches r = %.4f" % (nm, ro)
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-(LIM + EXTB), LIM + EXTT)
    ax.set_aspect("equal")
    ax.axis("off")


def logit(x):
    x = np.clip(x, 1e-12, 1 - 1e-12)
    return np.log(x / (1 - x))


def expit(y):
    return 1.0 / (1.0 + np.exp(-y))


def draw(which, stem, scale="linear", mode="share"):
    """Draw the figure of run `which` (F1 or F2) to HERE/<stem>.pdf/.png (HERE/other/<stem>_logit.* for scale="logit").

    mode="share": the share of the population in each atom (M1 Figure 3, M1 SI Figure 2);
    mode="enrichment": the share divided by the atom's share of the 16 strategies, log10, on the
    diverging map with white at 1 and the range 10^-4 .. 16 (M1 Figure 4, M1 SI Figure 3).
    """
    assert scale in ("linear", "logit") and mode in ("share", "enrichment"), (scale, mode)
    enrich = mode == "enrichment"
    assert not (enrich and scale != "linear"), "the enrichment figure has one scale"
    linear = scale == "linear"
    D = finalpanels.DATASETS[which]
    pts, S, N, UV = figF4.collect(which)         # per cell: the eight shares and the eight set sizes
    cells = vor.cells()
    assert S.shape == (512, 8) and not S[:, 5].any() and not N[:, 5].any()
    assert np.allclose(S.sum(1), 1.0)
    if enrich:
        norm = TwoSlopeNorm(vmin=LO, vcenter=0.0, vmax=HI)
        def colour(v):
            return ECMAP(norm(v))
    elif linear:
        norm = Normalize(vmin=0.0, vmax=1.0)
        def colour(v):
            return CMAP(norm(v))
    else:
        norm = FuncNorm((logit, expit), vmin=FLOOR, vmax=CEIL)
        def colour(v):
            return CMAP(norm(min(v, CEIL))) if v >= FLOOR else CMAP.get_under()

    keep = [k for k, i in enumerate(pts) if len(cells[int(i)]) >= 3]
    pol = [cells[int(pts[k])] for k in keep]
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    # -------------------------------------------------------------- the seven panels
    IDX = [int(c, 2) for c in ORDER]
    GAP = N[:, IDX] == 0
    if enrich:                                   # log10 of share / (|atom| / 16), undefined where the atom is empty
        V = np.where(GAP, np.nan, np.log10(np.maximum(S[:, IDX] * NSTRAT / np.maximum(N[:, IDX], 1), 1e-300)))
        nclip = int(((V[~GAP] < LO) | (V[~GAP] > HI)).sum())   # LO is no longer -HI
        V = np.clip(V, LO, HI)
        print("%s run (N = %d, beta = %d, mu = %s): enrichment per atom, 512 games%s"
              % (which, D["N"], D["beta"], D["muplain"], "" if nclip == 0 else ";  %d cells beyond the range, drawn at its ends" % nclip))
        print("%-4s %8s %8s %8s %8s  %6s  %6s  %6s" % ("atom", "geomean", "median", "min", "max", "empty", "r>1", "winner"))
    else:
        V = S[:, IDX]
        print("%s run (N = %d, beta = %d, mu = %s): share of the population per atom, 512 games" % (which, D["N"], D["beta"], D["muplain"]))
        print("%-4s %8s %8s %8s %8s  %6s  %6s  %6s" % ("atom", "mean", "median", "min", "max", "empty", "<0.1%", "winner"))
    arg = np.nanargmax(V, 1)
    for i, code in enumerate(ORDER):
        row, col = divmod(i, NCOL)
        vals, gap = V[:, i], GAP[:, i]
        ax = cell(fig, row, col)
        fc = [GREY if gap[k] else colour(vals[k]) for k in keep]
        ax.add_collection(PolyCollection(pol, facecolors=fc, edgecolors=fc, linewidths=0.3, antialiased=True, zorder=1))
        furniture(ax, "abcdefg"[i], r"$\mathbf{%s}$   %s" % (code, NAME[code]))
        ok = vals[~gap]
        if enrich:
            print("%-4s %8.3f %8.3f %8.3f %8.3f  %6d  %6d  %6d" % (code, 10 ** ok.mean(), 10 ** np.median(ok), 10 ** ok.min(), 10 ** ok.max(),
                                                                  gap.sum(), (ok > 0).sum(), (arg == i).sum()))
        else:
            print("%-4s %8.4f %8.4f %8.5f %8.4f  %6d  %6d  %6d" % (code, vals.mean(), np.median(vals), ok.min(), vals.max(),
                                                                  gap.sum(), (ok < FLOOR).sum(), (arg == i).sum()))

    # -------------------------------------------------------------- the eighth cell: the most abundant (enriched) atom
    SUPER = "enriched" if enrich else "abundant"
    ax = cell(fig, 1, 3)
    fc = [WINCOL[ORDER[arg[k]]] for k in keep]
    ax.add_collection(PolyCollection(pol, facecolors=fc, edgecolors=fc, linewidths=0.3, antialiased=True, zorder=1))
    furniture(ax, "h", "most %s atom" % SUPER)
    legend = []
    for i, code in enumerate(ORDER):
        xi = X0 + i * PITCH
        ax.add_patch(Rectangle((xi, YL - 0.5 * SW), SW, SW, fc=WINCOL[code], ec="0.25", lw=0.6, zorder=12, clip_on=False))
        legend.append(ax.text(xi + SW + 0.02, YL, code, fontsize=LFS, fontweight="bold", color="0.12", ha="left", va="center", zorder=12))
    u_, v_ = UV[:, 0], UV[:, 1]
    quad = np.where((u_ < 0) & (v_ > 0), "PD", np.where((u_ > 0) & (v_ > 0), "SD", np.where((u_ < 0) & (v_ < 0), "SH", "HG")))
    nwin = len(set(arg.tolist()))
    print("panel h, the most %s atom: " % SUPER + ";  ".join("%s in %d (%s)" % (c, (arg == i).sum(), " ".join("%s %d" % (q, ((arg == i) & (quad == q)).sum())
          for q in ["PD", "SD", "SH", "HG"] if ((arg == i) & (quad == q)).any())) for i, c in enumerate(ORDER) if (arg == i).any()))
    srt = np.sort(np.where(np.isnan(V), -np.inf, V), 1)
    margin = srt[:, -1] - srt[:, -2]
    if enrich:
        print("winner's log10 enrichment: median %.2f, min %.2f;  margin over the runner-up (decades): median %.2f, below 0.1 in %d games, below 0.3 in %d"
              % (np.median(srt[:, -1]), srt[:, -1].min(), np.median(margin), (margin < 0.1).sum(), (margin < 0.3).sum()))
    else:
        print("winner's share: median %.3f, min %.3f;  margin over the runner-up: median %.3f, below 0.05 in %d games, below 0.02 in %d"
              % (np.median(srt[:, -1]), srt[:, -1].min(), np.median(margin), (margin < 0.05).sum(), (margin < 0.02).sum()))

    # -------------------------------------------------------------- the bottom line: the colour bar of a-g, level with the legend of h
    CBX0, CBX1 = ML + 0.72 * PW, ML + 2 * PW + COLGAP - 0.22 * PW              # inches: under the inner halves of columns 1-2
    def yin(y):                                  # data y of a row-2 disk -> inches from the sheet's bottom edge
        return MB + FOOT + (y + LIM + EXTB) * PW / (2 * LIM)
    cax = fig.add_axes([CBX0 / FIG_W, (yin(YL) - 0.5 * CBH) / FIG_H, (CBX1 - CBX0) / FIG_W, CBH / FIG_H])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=ECMAP if enrich else CMAP), cax=cax, orientation="horizontal",
                      extend="neither" if linear else "min", extendfrac=0.04)
    if enrich:
        # powers of ten below 1, powers of two above: 16/n is what an atom of n members
        # reaches when it holds the whole population, so 16, 8, 4, 2 are the ceilings of
        # an atom of 1, 2, 4, 8 strategies
        cb.set_ticks([-4.0, -3.0, -2.0, -1.0, 0.0, np.log10(2), np.log10(4), np.log10(8), np.log10(16)])
        cb.ax.set_xticklabels(["$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "0.1", "1", "2", "4", "8", "16"])
        cb.ax.axvline(0.0, color="0.15", lw=1.0)
    elif linear:
        cb.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]); cb.ax.set_xticklabels(["0", "20%", "40%", "60%", "80%", "100%"])
    else:
        cb.set_ticks([1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99]); cb.ax.set_xticklabels(["0.1%", "1%", "10%", "50%", "90%", "99%"])
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=FS * 0.82, length=5, width=1.0)
    cb.outline.set_linewidth(1.0)
    label = "enrichment (a–g)" if enrich else "share of the population (a–g)" + ("" if linear else "\nlogit scale; darkest: below 0.1%")
    cax.text(-0.025 - (0.0 if linear else 0.045), 0.5, label, transform=cax.transAxes, ha="right", va="center",   # clear of the "min" arrow
             fontsize=FS * 0.95, color="0.12", linespacing=1.25)

    # -------------------------------------------------------------- checks: nothing off the sheet, nothing overlapping
    fig.canvas.draw()
    R = fig.canvas.get_renderer()
    def on_sheet(t):
        bb = t.get_window_extent(R)
        assert 0 < bb.x0 and bb.x1 < FIG_W * fig.dpi and 0 < bb.y0 and bb.y1 < FIG_H * fig.dpi, "text runs off the sheet: %r" % t.get_text()
    for t_ in fig.texts:                        # figure-level texts are not in any axes (none at present)
        on_sheet(t_)
    for ax_ in fig.axes:
        for t_ in list(ax_.texts) + (ax_.get_xticklabels() + ax_.get_yticklabels() if ax_.axison else []):
            if t_.get_text() and t_.get_visible():
                on_sheet(t_)
        heads = [t_ for t_ in ax_.texts if t_.get_text() and t_.get_va() == "bottom" and t_.get_ha() in ("center", "left")]
        for i_ in range(len(heads)):
            for j_ in range(i_ + 1, len(heads)):
                assert not heads[i_].get_window_extent(R).overlaps(heads[j_].get_window_extent(R)), \
                    "texts overlap: %r / %r" % (heads[i_].get_text(), heads[j_].get_text())
    for i_ in range(len(legend) - 1):
        nxt = ax.transData.transform((X0 + (i_ + 1) * PITCH, YL))[0]
        assert legend[i_].get_window_extent(R).x1 < nxt - 2, "legend entry %d runs into the next swatch" % i_

    out = os.path.join(HERE, stem + ".pdf") if linear else os.path.join(HERE, "other", os.path.basename(stem) + "_logit.pdf")
    # stem may carry a subfolder ("SI/SIFig2"); the logit variants always go to other/
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=150)
    plt.close(fig)
    print("%s: %.2f x %.2f in, pdf %.1f MB" % (os.path.basename(out), FIG_W, FIG_H, os.path.getsize(out) / 1e6))
    return arg
