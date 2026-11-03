#!/usr/bin/env python3
"""
atomshares.py -- the machinery of Figures 3 and 4 and SI Figures 2 and 3: the seven atoms of the three properties
(efficient, nash, competitive) over the disk of games in a Wright-Fisher run of DiskM2WF.

    Figure 3      fig3.py      run F1:  N = 1000,  beta = 100,  mu = 1e-2,  eps = 1e-4   share of the population
    Figure 4      fig4.py      run F1                                                    enrichment
    SI Figure 2   sifig2.py    run F2:  N = 100,   beta = 3,    mu = 1e-4,  eps = 1e-4   share of the population
    SI Figure 3   sifig3.py    run F2                                                    enrichment
    (Martin, 2026-09-17: main text = Figures 1, 2, 3, 4; the N = 100 pages are SI Figures 2 and 3)
    From 2026-09-19 (Martin's plan) those four pages are retired to other/ (fig3_old.py, fig4_old.py, sifig2_old.py,
    sifig3_old.py, still drawn by draw()).  The SI pages are now draw_runs(): two runs on one sheet of four rows,
    the seven atom shares and the most abundant atom for N = 100 (rows 1-2) and N = 1000 (rows 3-4) --
    SI Figure 2  sifig2.py    memory one (DiskM1WF's runs, read through wfdata.py)
    SI Figure 3  sifig3.py    memory two (DiskM2WF's runs)

The entry points call draw(which, stem, scale, mode); layout, colours and checks live here, so the
figures cannot drift apart.  mode="share" is the share of the population held by the atom (Figure 3,
SI Figure 2; viridis, linear); mode="enrichment" is that share divided by the atom's share of the 65536
strategies (r_abc = p_abc / q_abc of Literature/3Traits/3a.tex), log10 on a diverging map centred at
1 (RdBu: red depleted, blue enriched, the convention of the Worth9 and Beta6 drafts), fixed range
10^-4 to 10^4; panel h is then the most enriched atom.

    row 1:   a  000     b  100     c  010     d  001
    row 2:   e  110     f  011     g  111     h  the most abundant atom at each game
    bottom:  the colour bar of a-g under columns 1-2, level with the swatch legend of h; no notes on the sheet
             (Martin, 2026-09-17: the run parameters and the meaning of grey belong in the figure legend)

Panels a-g: the time-averaged share of the population held by the strategies of one atom at each of
the 512 games (the Voronoi cells of the Fermat sunflower), viridis on a linear scale from 0 to 1
(Martin's choice, 2026-09-17; scale="logit" gives 0.1% to 99% with the darkest colour for shares
below 0.1%, written to other/ -- the subfolder for variants that are not final figures).  Cells where the atom has no strategy at that game are grey, as in Figure 2.  Panel h
colours each game by the atom with the largest share, with a swatch legend under the disk: blues for
the two efficient atoms (100 light, 110 dark), reds for the two competitive ones (001 light, 011
dark), purple for 111, gold for 010 and dark grey for 000 -- the hues of the Euler diagram of
Figure 1, saturated.  The atom of each strategy at each game is decided in the eps -> 0 limit
(DiskM2WF/Opt/exact.masks, through figF4.collect, which also checks that the eight shares sum to one
and that 101 is empty at all 512 games).  Frame and furniture are those of Figure 2; the projection
is the same disk (K = 4).  Reads DiskM2WF/Opt read-only.
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
from matplotlib.colors import FuncNorm, Normalize
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, os.pardir, "Figures")
OPT = "/Users/martin/Documents/CL1/DiskM2WF/Opt"
sys.path.insert(0, FIGDIR)
sys.path.insert(0, OPT)
import disklib                                                     # noqa: E402
from disklib import NAMES, arc_text                                # noqa: E402
import figF4, vor, finalpanels                                     # noqa: E402

matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 15, "axes.linewidth": 1.2})
K = 4.0
assert finalpanels.K == K


def phi(p):
    return disklib.phi(p, K)


ORDER = ["000", "100", "010", "001", "110", "011", "111"]
NAME = {"000": "none of the three", "100": "efficient only", "010": "stable only", "001": "competitive only",
        "110": "efficient and stable only", "011": "stable and competitive only", "111": "efficient, stable and competitive"}
GREY = "0.84"                                    # no strategy in the atom at that game (Figure 2's grey)
WINCOL = {"000": "0.36", "100": "#8fb8e8", "010": "#dba622", "001": "#f0a48c",     # panel h: the most abundant atom
          "110": "#2f5fb3", "011": "#c43d3b", "111": "#7a4fa8"}
HEAD_DX = {"111": 0.03}   # data units: the head of the 111 panels, two letters longer since "Nash" became "stable",
                           # would come within 0.08 in of its panel letter if centred; 0.03 (0.055 in) to the right
                           # gives back the 0.13 in it had
CMAP = plt.get_cmap("viridis").copy()
CMAP.set_under("#270031")                        # logit scale: present, but holding less than FLOOR
FLOOR, CEIL = 1e-3, 0.99
LO, HI = -4.0, 4.0                               # enrichment: log10 range of the diverging map, centred at 1
ECMAP = plt.get_cmap("RdBu")
NSTRAT = 65536.0

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

    mode="share": the share of the population in each atom (Figure 3, SI Figure 2); mode="enrichment": the
    share divided by the atom's share of the 65536 strategies, log10, on the diverging map (Figure 4, SI Figure 3).
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
        norm = Normalize(vmin=LO, vmax=HI)
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
    if enrich:                                   # log10 of share / (|atom| / 65536), undefined where the atom is empty
        V = np.where(GAP, np.nan, np.log10(np.maximum(S[:, IDX] * NSTRAT / np.maximum(N[:, IDX], 1), 1e-300)))
        nclip = int((np.abs(V[~GAP]) > HI).sum())
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
        cb.set_ticks(list(range(-4, 5)))
        cb.ax.set_xticklabels(["$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "0.1", "1", "10", "100", "$10^{3}$", "$10^{4}$"])
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


def draw_runs(runs, stem, labels):
    """Two (or more) runs on one sheet, two rows of four per run: the shares of the seven atoms (viridis, linear, grey
    where the atom has no strategy) and the most abundant atom, as draw() draws one run; the colour bar and the swatch
    legend once, at the bottom.  runs = [(space, which), ...] read through wfdata.load ("m1" or "m2", "F1" or "F2");
    labels = the text printed rotated at the left of each block.  Writes HERE/<stem>.pdf and .png."""
    import wfdata                                                  # here, not at the top: wfdata imports this module
    NR = 2 * len(runs)
    LABW = 0.75
    W = ML + LABW + NCOL * PW + (NCOL - 1) * COLGAP + MR
    H = MB + NR * PH + (NR - 1) * ROWGAP + MT
    fig = plt.figure(figsize=(W, H))

    def xy(row, col):
        return (ML + LABW + col * (PW + COLGAP), MB + (NR - 1 - row) * (PH + ROWGAP))

    def cell(row, col):
        x, y = xy(row, col)
        return fig.add_axes([x / W, y / H, PW / W, PH / H])

    pol, keep = wfdata.cells()
    norm = Normalize(vmin=0.0, vmax=1.0)
    letters = iter("abcdefghijklmnopqrstuvwx")
    for b, (space, which) in enumerate(runs):
        S, N, E, UV, D = wfdata.load(space, which)                 # (512, 7) in ORDER, in cell order 0..511
        GAP = N == 0
        arg = np.where(GAP, -np.inf, S).argmax(1)
        print("%s %s (N = %d, beta = %d, mu = %s): share of the population per atom, 512 games" % (space, which, D["N"], D["beta"], D["muplain"]))
        print("%-4s %8s %8s %8s %8s  %6s  %6s" % ("atom", "mean", "median", "min", "max", "empty", "winner"))
        for i, code in enumerate(ORDER):
            row, col = divmod(i, NCOL)
            ax = cell(2 * b + row, col)
            fc = [GREY if GAP[k, i] else CMAP(norm(S[k, i])) for k in keep]
            ax.add_collection(PolyCollection(pol, facecolors=fc, edgecolors=fc, linewidths=0.3, antialiased=True, zorder=1))
            head = r"$\mathbf{%s}$   %s" % (code, NAME[code])
            furniture(ax, next(letters), head)
            if code in HEAD_DX:                                    # the longest head, clear of the panel letter
                (t_,) = [t_ for t_ in ax.texts if t_.get_text() == head]
                t_.set_x(HEAD_DX[code])
            ok = S[~GAP[:, i], i]
            print("%-4s %8.4f %8.4f %8.5f %8.4f  %6d  %6d" % (code, S[:, i].mean(), np.median(S[:, i]), ok.min() if len(ok) else np.nan,
                                                            S[:, i].max(), GAP[:, i].sum(), (arg == i).sum()))
        ax = cell(2 * b + 1, 3)
        fc = [WINCOL[ORDER[arg[k]]] for k in keep]
        ax.add_collection(PolyCollection(pol, facecolors=fc, edgecolors=fc, linewidths=0.3, antialiased=True, zorder=1))
        furniture(ax, next(letters), "most abundant atom")
        srt = np.sort(np.where(GAP, -np.inf, S), 1)
        print("most abundant atom: " + ", ".join("%s %d" % (c, (arg == i).sum()) for i, c in enumerate(ORDER) if (arg == i).any())
              + ";  its share: median %.3f; decided by < 0.05 at %d games" % (np.median(srt[:, -1]), ((srt[:, -1] - srt[:, -2]) < 0.05).sum()))
        x0, y0 = xy(2 * b, 0); x1, y1 = xy(2 * b + 1, 0)
        yc = 0.5 * (y0 + y1) + PH * (LIM + EXTB) / (2 * LIM + EXTT + EXTB)
        fig.text((x0 - 0.45) / W, yc / H, labels[b], rotation=90, ha="center", va="center", fontsize=FS * 1.25, color="0.12")
    # the bottom line, in the last row's own bottom margin: the swatch legend under the last panel h, the colour bar under columns 1-2
    legend = []
    for i, code in enumerate(ORDER):
        xi = X0 + i * PITCH
        ax.add_patch(Rectangle((xi, YL - 0.5 * SW), SW, SW, fc=WINCOL[code], ec="0.25", lw=0.6, zorder=12, clip_on=False))
        legend.append(ax.text(xi + SW + 0.02, YL, code, fontsize=LFS, fontweight="bold", color="0.12", ha="left", va="center", zorder=12))
    CBX0, CBX1 = ML + LABW + 0.72 * PW, ML + LABW + 2 * PW + COLGAP - 0.22 * PW
    ylast = xy(NR - 1, 0)[1]
    ycb = ylast + (YL + LIM + EXTB) * PW / (2 * LIM)
    cax = fig.add_axes([CBX0 / W, (ycb - 0.5 * CBH) / H, (CBX1 - CBX0) / W, CBH / H])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]); cb.ax.set_xticklabels(["0", "20%", "40%", "60%", "80%", "100%"])
    cb.ax.minorticks_off(); cb.ax.tick_params(labelsize=FS * 0.82, length=5, width=1.0); cb.outline.set_linewidth(1.0)
    cax.text(-0.025, 0.5, "share of the population (the atom panels)", transform=cax.transAxes, ha="right", va="center", fontsize=FS * 0.95, color="0.12")
    # checks: nothing off the sheet, no two heads overlapping, the legend entries apart
    fig.canvas.draw()
    R = fig.canvas.get_renderer()
    def on_sheet(t):
        bb = t.get_window_extent(R)
        assert 0 < bb.x0 and bb.x1 < W * fig.dpi and 0 < bb.y0 and bb.y1 < H * fig.dpi, "text runs off the sheet: %r" % t.get_text()
    for t_ in fig.texts:
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
    out = os.path.join(HERE, stem + ".pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=150)
    plt.close(fig)
    print("%s: %.2f x %.2f in, pdf %.1f MB" % (os.path.basename(out), W, H, os.path.getsize(out) / 1e6))
