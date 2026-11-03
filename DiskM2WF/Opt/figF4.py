"""
F4 -- the seven classes of the (efficient, Nash, rival) decomposition, and
where each of them holds its abundance on the disk.

    python3 figF4.py [outdir]          (default ../FiguresFinal)
    python3 figF6.py [outdir]          the same decomposition on the F2 run

Martin, 2026-09-14: "We have 3 sets of strategies: efficient (eps -> 0), NE
(eps -> 0), rival.  These 3 sets define 7 classes ... note that 101 is not
possible (because an efficient rival is a NE).  Make a figure with 7 panels,
show for each of the 7 classes the total frequency on the disk."

THE DECOMPOSITION.  `Opt/exact.py` returns three boolean masks over the 65536
memory-two strategies at each game (u, v): `eff`, `ne` and `riv`.  Label each
strategy with the three bits in that order, so 110 is "efficient and a Nash
equilibrium but not a rival" -- the PARTNER of the rest of this project -- and
111 is the FRIENDLY RIVAL.  Every strategy falls in exactly one of the eight
labels, so the eight shares of the Wright-Fisher abundance sum to 1 at every
cell.  That is asserted here cell by cell, not assumed.

WHY SEVEN AND NOT EIGHT.  101 -- efficient and rival but not an equilibrium --
is empty.  Martin gives the reason: an efficient rival is a Nash equilibrium.
This figure does not take that on trust.  `collect()` checks, at every one of
the 512 cells, that the 101 SET is empty, and raises rather than draw if it
ever is not.  On both runs it is empty at all 512.

THE EIGHTH PANEL held a Venn key until 2026-09-14, when Martin asked for "the
most abundant type (in terms of total abundance)" there instead.  It is now
the argmax over the seven -- the question F1's panel c asks of a coarser
five-type scheme, and this map is a strict REFINEMENT of that one: F1's
"rival" is 001 + 011 and its "neither" is 000 + 010, the other three match one
to one.  No figure has a Venn key now; F8 draws the lattice of named unions
instead, which is what the Venn was trying to say.

WHAT EACH PANEL SHOWS.  The 512 Voronoi cells of the Fermat sunflower, each
coloured by the fraction of that cell's F1 abundance that sits in the class --
the same cells, frame and two lines as F1's panel c.  All seven panels share
one logarithmic colour scale, because the shares run from below 10^-6 to 0.99
and a linear scale would leave four of the seven blank.  The linear headline is
printed in each panel's subtitle instead.

READING IT.  The two large classes are 110 (partner, disk-mean 58%) and 111
(friendly rival, 36%); between them they hold 94% of the abundance everywhere.
But none of the five small classes is uniformly small: 100 reaches 83% of a
cell in N, 011 reaches 53% in W, and 000 -- which is 84% of the strategy space
-- is present in every single cell and dominant in none.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpat
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir, "Figures"))
import disk, exact, finalpanels, vor, wfload      # noqa: E402

# the eight labels, bit 2 = efficient, bit 1 = NE, bit 0 = rival
CODE = {0: "000", 4: "100", 2: "010", 1: "001", 6: "110", 3: "011", 7: "111",
        5: "101"}
# Martin's titles, 2026-09-14, one line each.  "unbeatable" is his word for
# what exact.py calls `riv`: no co-player ever earns more.  The two names are
# joined in the caption so the figure and the rest of the project connect.
NAME = {0: "none of the three",
        4: "efficient only",
        2: "nash only",
        1: "unbeatable only",
        6: "efficient + nash",
        3: "nash + unbeatable",
        7: "all three",
        5: "efficient + unbeatable, not nash"}
# Martin's order: none, the three singletons, then the pairs and the triple
ORDER = [0, 4, 2, 1, 6, 3, 7]
IMPOSSIBLE = 5
# which figure this decomposition is called when drawn on which WF run.  The
# two runs are the same 512 games under different process parameters, so the
# colour scale is deliberately FIXED rather than fitted per dataset: F4 and F6
# are meant to be laid side by side.
STEM = {"F1": "F4", "F2": "F6"}

FLOOR = 1e-4               # the bottom of the log scale
CMAP = plt.get_cmap("viridis").copy()
# A VERY DARK PURPLE, which is viridis's own bottom (#440154) darkened to 35%
# of its linear luminance.  Martin, 2026-09-14: "is there a very very dark
# purple?  (the category distinction between 10^-4 or lower is not that
# relevant).  black looks a bit odd there."  So this should read as the ramp
# continuing one step below its bottom, not as a separate category: dE 9.4
# from viridis(0) in normal vision, and still dE 22 from black, so it is
# plainly purple and plainly not a hard break.  It was light grey until this
# morning, which read as MORE rather than less because viridis runs dark to
# light, and then black, which read as a category of its own.
CMAP.set_under("#270031")  # present, but holding less than FLOOR
# A class that does not EXIST at a game is a different statement from a class
# that exists and is rare, and the first version of this figure painted both
# the same grey.  Martin, 2026-09-14: "what is the meaning of light grey
# areas?"  They meant both.  Now they mean only the second, and this is the
# first -- shared with F5/F7, which need the distinction even more sharply
# because there the average over an empty class is not a number at all.
# Explained in the captions, not on the figures, since 2026-09-14.
GAP = dict(facecolor="white", edgecolor="#9a9a9a", lw=0.35)

# Panel 8's categorical palette.  This map is a strict REFINEMENT of F1's
# panel c -- that figure's five types are exactly these seven with 001+011
# merged into "rival" and 000+010 merged into "neither" -- so partner keeps
# F1's blue, the friendly rivals keep its yellow, 000 keeps its grey and 100
# stays in the green family.  The other three are chosen, not inherited:
# Opt/cvd.py puts the worst pair of this palette at dE 10.3 across normal,
# deuteranopia, protanopia and tritanopia, where the five-class palette F1 and
# F3 share manages only 1.4 (PROJECT.md 14.5).
BIGCOL = {0: "#9a9a9a", 4: "#7fbf7f", 2: "#a6761d", 1: "#543005",
          6: "#1f5fbf", 3: "#000000", 7: "#f2c200"}


def collect(which="F1"):
    """per cell: the eight class shares and the eight set sizes

    Two premises are checked here rather than asserted: that 101 is empty, and
    that the eight shares are a partition of the abundance.
    """
    d = exact.load()
    G = wfload.load_run(finalpanels.idxpath(which), ie=3)
    pts, S, N, UV = [], [], [], []
    for g in G.values():
        ipt = vor.index_of(g["tag"])
        if ipt is None:
            continue
        u, v = g["u"], g["v"]
        eff, riv, ne = exact.masks(d, u, v)
        pi = np.asarray(g["pi"], float)
        pi = pi / pi.sum()
        lab = eff.astype(np.int8) * 4 + ne.astype(np.int8) * 2 + riv
        n = np.array([int((lab == c).sum()) for c in range(8)])
        if n[IMPOSSIBLE]:
            raise ValueError(
                "%s (u=%.4f, v=%.4f): the 101 class is NOT empty -- %d "
                "strategies are efficient and rival but not Nash.  The figure's"
                " premise fails; do not draw." % (g["tag"], u, v, n[IMPOSSIBLE]))
        s = np.array([pi[lab == c].sum() for c in range(8)])
        if abs(s.sum() - 1.0) > 1e-9:
            raise ValueError("%s: the eight shares sum to %.12f, not 1"
                             % (g["tag"], s.sum()))
        pts.append(int(ipt)); S.append(s); N.append(n); UV.append((u, v))
    return (np.array(pts), np.array(S), np.array(N), np.array(UV))


def panel(ax, cells, pts, val, code, name, n, uv, norm):
    pol, fc, gap = [], [], []
    for ipt, x, ni in zip(pts, val, n):
        q = cells[int(ipt)]
        if len(q) < 3:
            continue
        if ni == 0:                      # the class is empty at that game
            gap.append(q)
            continue
        pol.append(q)
        fc.append(CMAP(norm(x)) if x >= FLOOR else CMAP.get_under())
    if gap:
        ax.add_collection(PolyCollection(gap, zorder=1, **GAP))
    ax.add_collection(PolyCollection(pol, facecolors=fc, edgecolors="none",
                                     zorder=1))
    frame(ax)
    title(ax, code, name)
    ax.text(0.0, -1.07, "%d" % n.min() if n.min() == n.max()
            else "%d–%d" % (n.min(), n.max()),
            ha="center", va="top", fontsize=9.0, color="0.38")


def frame(ax):
    """every one of the eight disks is drawn on exactly this box"""
    disk.draw_frame(ax, labels=False, quad=False)
    ax.set_xlim(-1.16, 1.16); ax.set_ylim(-1.23, 1.25)
    ax.set_aspect("equal"); ax.axis("off")


def title(ax, code, name):
    ax.text(0.0, 1.07, r"$\mathbf{%s}$   %s" % (code, name), ha="center",
            va="bottom", fontsize=11.5, color="0.10")


def pct(x):
    if x >= 0.0995:
        return "%.0f%%" % (100 * x)
    if x >= 0.00995:
        return "%.1f%%" % (100 * x)
    return "%.2f%%" % (100 * x)


def winners(ax, cells, pts, S):
    """panel 8: which of the seven holds the most abundance at each cell

    The seven panels each answer "where is this class?"; this one answers
    "and which of them wins?", which is the question F1's panel c asks of a
    coarser five-type scheme.  Replaces the Venn key that stood here until
    2026-09-14, and F5/F7's eighth panel was changed to match (the argmax of
    their own quantity, not of this one).  Drawn on `frame`, the same box as
    the other seven, so all eight disks are the same size.
    """
    big = S[:, ORDER].argmax(axis=1)
    pol, fc = [], []
    for ipt, j in zip(pts, big):
        q = cells[int(ipt)]
        if len(q) >= 3:
            pol.append(q)
            fc.append(BIGCOL[ORDER[int(j)]])
    ax.add_collection(PolyCollection(pol, facecolors=fc, edgecolors="none",
                                     zorder=1))
    frame(ax)
    ax.text(0.0, 1.07, "most abundant class", ha="center", va="bottom",
            fontsize=11.5, fontweight="bold", color="0.10")
    # nothing under this disk: Martin, 2026-09-14, "remove '101 empty
    # everywhere'.  this will be explained in the figure legend."  It is, in
    # the caption's first paragraph and again in the notes block.
    return np.bincount(big, minlength=len(ORDER))


def main(outdir=None, which="F1"):
    outdir = outdir or os.path.join(HERE, os.pardir, "FiguresFinal")
    stem = STEM[which]
    D = finalpanels.DATASETS[which]
    pts, S, N, UV = collect(which)
    norm = LogNorm(vmin=FLOOR, vmax=1.0)
    cells = vor.cells()

    # No suptitle: Martin, 2026-09-14, "remove top line (this should go into
    # the Figure legend)".  The run's parameters live in the caption now.
    fig = plt.figure(figsize=(15.2, 9.4))
    gs = fig.add_gridspec(2, 4, left=0.010, right=0.990, top=0.985,
                          bottom=0.125, wspace=0.03, hspace=0.16)
    for i, c in enumerate(ORDER):
        ax = fig.add_subplot(gs[i // 4, i % 4])
        panel(ax, cells, pts, S[:, c], CODE[c], NAME[c], N[:, c], UV, norm)
    cnt = winners(fig.add_subplot(gs[1, 3]), cells, pts, S)

    # the ramp on the LEFT, panel 8's key to the right of it
    cax = fig.add_axes([0.045, 0.062, 0.27, 0.019])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                      orientation="horizontal", extend="min")
    cb.set_ticks([1e-4, 1e-3, 1e-2, 1e-1, 1.0])
    cb.set_ticklabels(["$10^{-4}$", "0.1%", "1%", "10%", "100%"])
    cb.ax.tick_params(labelsize=8.0, length=3.0, pad=3)
    cb.outline.set_linewidth(0.5)
    cb.set_label("frequency", fontsize=10.0, color="0.20", labelpad=5)
    # bold class, regular count -- mathtext, the same device as the titles
    h = [mpat.Patch(facecolor=BIGCOL[c], edgecolor="0.35", lw=0.4,
                    label=r"$\mathbf{%s}$  %d" % (CODE[c], cnt[i]))
         for i, c in enumerate(ORDER)]
    # matplotlib fills a legend column by column; reorder so that reading
    # ACROSS the rows gives ORDER, which is the order of the seven panels
    nc, n_ = 7, len(h)
    nr = -(-n_ // nc)
    h = [h[r * nc + c] for c in range(nc) for r in range(nr)
         if r * nc + c < n_]
    fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.405, 0.028),
               ncol=nc, fontsize=9.0, frameon=False, handlelength=1.2,
               columnspacing=1.3, handletextpad=0.45,
               title="most abundant class: cells won, of %d" % len(pts),
               title_fontproperties={"size": 9.5})

    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, stem + ".pdf")
    fig.savefig(p); fig.savefig(p[:-4] + ".png", dpi=200)
    plt.close(fig)
    print("wrote", p)
    print("  %s on the %s run: N=%d beta=%d u=%s" % (stem, which, D["N"],
                                                     D["beta"], D["muplain"]))
    print("  %-4s %-28s %9s %9s %9s   %s"
          % ("code", "class", "mean", "median", "max", "set size"))
    for c in ORDER + [IMPOSSIBLE]:
        n = N[:, c]
        print("  %-4s %-28s %9.5f %9.5f %9.5f   %d–%d"
              % (CODE[c], NAME[c].replace("\n", ", "), S[:, c].mean(),
                 np.median(S[:, c]), S[:, c].max(), n.min(), n.max()))
    print("  shares sum to 1 at every cell: %s   (max |sum-1| = %.2e)"
          % (np.allclose(S.sum(axis=1), 1.0), np.abs(S.sum(axis=1) - 1).max()))
    return pts, S, N, UV


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
