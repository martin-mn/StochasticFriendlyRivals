#!/usr/bin/env python3
r"""Figure 3 of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity":
the memory-two cube in the Prisoner's Dilemma wedge, and how rare friendly rivals are.

    PYTHONDONTWRITEBYTECODE=1 python3 fig3.py      ->  Fig3.pdf, Fig3.png (dpi 200)

Panel a  The friendly rivals of wedge W = {u+v<1, T>S} inside Memory-2, stochastic ([0,1]^16) form
         two separated families (proved: PRFull/check/wsep; whether each family is connected is open).
         Two 4x4 probability tables (row = most recent outcome, column = the
         outcome before, focal player's action first, * = free in (0,1)): left the unique maximal
         stratum of the binary-anchored family (the solid 5-cube whose 3-face {p1=1, p5=0}
         carries the eight binary friendly rivals), right the maximal 7-dimensional stratum
         1*1*001*1*1**0*0 of the binary-invisible family (the one holding the exact certificate
         of PartnersRivalsM2f/Text/tables.tex, Table 3a).  Outlined: p6, p8 (0 on every stratum of
         one family, 1 on every stratum of the other) and p9 (free on every stratum of the
         binary-invisible family).  Both tables and the counts 126 / 238 are recomputed below
         from the census and asserted; the counts are exact: every counted stratum holds an
         exactly certified friendly rival (PRFull/cannon/stratcert, one exact rational point per
         stratum), and every other efficient stratum of W is excluded by the exclusion
         certificates (PRFull/check/ubound), so no other stratum holds a friendly rival of W.
Panel b  Dimension histogram of the 364 wedge-W strata, stacked binary-anchored / binary-invisible.
Panel c  The one-probability witness p(t) = (0,0,0,t,0,0,0,0,0,1,0,1,1,0,1,0) -- a friendly rival
         in wedge N (NOT W: its self-play limit is perfect alternation) for every t in (0,1).
         The self-play limit frequencies are computed in exact rationals by
         PartnersRivalsM2f/check/exact.py and asserted: (0,1/2,1/2,0) on (0,1),
         (1/6,1/3,1/3,1/6) at t=0, (2/9,1/3,1/3,1/9) at t=1.  The competitive margin over all 65536
         binary co-players is checked separately by fig3_witness_check.py (-> .out).
Panel d  Rate per draw of a friendly rival in W, N, E, S under the three-outcome kernel at
         alpha = 1/3 exactly (log scale), with the independence level P(efficient) x P(competitive).
         Values of SI Table 3 of Text/ms.tex, cross-checked here against a recomputation from the
         512 task files of the sweep at alpha = 1/3 (PRFull cannon/sweep13, array 48928780;
         it replaced the sweep at alpha = 0.333 of PartnersRivalsM2f).  Beside it, a locator disk of
         the plane of games (reference furniture), each wedge shaded by the same rate.

Sources (all READ-ONLY):
  /Users/martin/Documents/CL1/PartnersRivalsM2f/check/st_eff.txt        support-stratum census
  /Users/martin/Documents/CL1/PartnersRivals/Figures/m2_masks.npz       the binary memory-two masks
  /Users/martin/Documents/CL1/PartnersRivalsM2f/check/exact.py          exact eps->0 engine
  /Users/martin/Documents/CL1/PRFull/cannon/sweep13/out_13/a_1o3_*.tsv  alpha = 1/3 sweep, 512 tasks
  /Users/martin/Documents/CL1/PartnersRivals1/Figures/disklib.py        disk projection, rim text
"""
import os
import sys
import collections
from fractions import Fraction as F

sys.dont_write_bytecode = True                   # nothing may appear in the read-only projects
M2F = "/Users/martin/Documents/CL1/PartnersRivalsM2f"
sys.path.insert(0, os.path.join(M2F, "check"))
sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivals1/Figures")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch

import exact                                     # noqa: E402  (PartnersRivalsM2f/check, read-only)
import disklib                                   # noqa: E402  (PartnersRivals1/Figures, read-only)

HERE = os.path.dirname(os.path.abspath(__file__))
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                            "font.size": 7.5, "axes.linewidth": 0.7, "xtick.labelsize": 7,
                            "ytick.labelsize": 7, "xtick.major.width": 0.7,
                            "ytick.major.width": 0.7, "xtick.major.size": 2.6,
                            "ytick.major.size": 2.6, "hatch.linewidth": 0.5})

# the atom palette of the reference paper (PartnersRivals1/FinalFigures/eulerlib.py, atomshares.py):
# 111 = efficient, stable and competitive = the friendly rivals
LILAC, PURPLE = "#cbbbe6", "#7a4fa8"             # TINT["111"], WINCOL["111"]
DEEP = "#4a2c6e"
VIS, INV = LILAC, PURPLE                         # binary-anchored / binary-invisible family
INK = "0.12"
ONE = "0.885"                                    # a coordinate pinned at 1 (0 is white)
HL_SEP, HL_P9 = "0.05", "#e08214"                # outlines: p6, p8 / p9
C_CC, C_DD, C_ALT = "#2166ac", "0.15", "#d95f02"  # panel c series

OUT = ("CC", "CD", "DC", "DD")
SYM = {0: "0", 1: "*", 2: "1"}

# =============================================================================================
# 1.  the wedge-W census (as PartnersRivalsM2f/Figures/fig2.py)
# =============================================================================================


def dig(g):
    d = []
    for _ in range(16):
        d.append(g % 3)
        g //= 3
    return d


def undig(d):
    g = 0
    for j in range(15, -1, -1):
        g = g * 3 + d[j]
    return g


def pat(g):
    return "".join(SYM[x] for x in dig(g))


def pat2g(s):
    return undig([{"0": 0, "*": 1, "1": 2}[ch] for ch in s])


def nfree(g):
    return sum(1 for x in dig(g) if x == 1)


CENSUS = set()
with open(os.path.join(M2F, "check", "st_eff.txt")) as fh:
    for line in fh:
        a = line.split()
        if len(a) >= 4 and a[1] == "C" and int(a[2]) > 0:     # self-play all-CC, >= 1 rival draw
            CENSUS.add(int(a[0]))
# the two strata the 64-draw census missed, recovered by the exhaustive 512-draw re-census
# (PartnersRivalsM2f/Text/tables.tex Table 2b: wedge W, census 362, deficit 2, tightened 364)
CORR = [pat2g("1*1*00111*1**0*0"), pat2g("1*1100111*1**0*0")]
assert not (set(CORR) & CENSUS)
FULL = CENSUS | set(CORR)

_M = np.load("/Users/martin/Documents/CL1/PartnersRivals/Figures/m2_masks.npz")
BINW = sorted(np.where(_M["effCC"] & _M["rivP"])[0].tolist())   # binary friendly rivals of W
assert BINW == [19079, 19087, 20103, 20111, 23175, 23183, 24199, 24207]


def visible(g):
    """the closed cell of stratum g contains a binary friendly rival"""
    d = dig(g)
    return any(all(d[j] == 1 or (d[j] == 2) == bool((b >> j) & 1) for j in range(16))
               for b in BINW)


VISS = sorted(g for g in FULL if visible(g))
INVS = sorted(g for g in FULL if not visible(g))
HV = collections.Counter(nfree(g) for g in VISS)
HI = collections.Counter(nfree(g) for g in INVS)
assert (len(CENSUS), len(FULL), len(VISS), len(INVS)) == (362, 364, 126, 238)
assert sum(HV.values()) + sum(HI.values()) == 364
assert max(HV) == 5 and HV[5] == 1 and max(HI) == 7 and HI[0] == 0

VIS_MAX = [g for g in VISS if nfree(g) == 5]
INV_MAX = sorted(pat(g) for g in INVS if nfree(g) == 7)
assert [pat(g) for g in VIS_MAX] == ["1*1*0*0101*1*010"]
assert INV_MAX == sorted(["11**001*1*1**0*0", "1*1*001*1*1**0*0"])
CELLS_V, CELLS_I = pat(VIS_MAX[0]), "1*1*001*1*1**0*0"     # right: the Table-3a stratum
_dv = dig(VIS_MAX[0])                               # the 8 binaries are vertices of its closed cell
assert all(all(_dv[j] == 1 or (_dv[j] == 2) == bool((b >> j) & 1) for j in range(16)) for b in BINW)
assert [j for j in range(16) if len({(b >> j) & 1 for b in BINW}) == 2] == [3, 10, 12]
assert all(((b >> 1) & 1, (b >> 5) & 1) == (1, 0) for b in BINW)   # on the 3-face {p1=1, p5=0}
for j, (vv, ii) in {6: ("0", "1"), 8: ("0", "1"), 9: ("1", "*")}.items():
    assert all(pat(g)[j] == vv for g in VISS) and all(pat(g)[j] == ii for g in INVS), j
ZERO_DIM = sorted(sum(1 << j for j, x in enumerate(dig(g)) if x == 2) for g in FULL if nfree(g) == 0)
assert ZERO_DIM == BINW                          # the 0-dimensional strata are the binaries

print("wedge W: census %d + recovered %d = %d strata;  binary-anchored %d, binary-invisible %d"
      % (len(CENSUS), len(CORR), len(FULL), len(VISS), len(INVS)))
print("  anchored  by dimension", dict(sorted(HV.items())))
print("  invisible by dimension", dict(sorted(HI.items())))
print("  maximal: anchored %s   invisible %s" % (CELLS_V, " ".join(INV_MAX)))

# =============================================================================================
# 2.  panel c: the witness, exact
# =============================================================================================


def Pseg(t):
    return [F(0), F(0), F(0), F(t), F(0), F(0), F(0), F(0),
            F(0), F(1), F(0), F(1), F(1), F(0), F(1), F(0)]


TGRID = [F(0)] + [F(k, 24) for k in range(1, 24)] + [F(1)]
SELF = {t: exact.freqs(Pseg(t), Pseg(t)) for t in TGRID}
W0, W1, WI = SELF[F(0)], SELF[F(1)], SELF[F(1, 2)]
assert all(SELF[t] == WI for t in TGRID if 0 < t < 1)
assert WI == [0, F(1, 2), F(1, 2), 0]
assert W0 == [F(1, 6), F(1, 3), F(1, 3), F(1, 6)]
assert W1 == [F(2, 9), F(1, 3), F(1, 3), F(1, 9)]
print("witness self-play: t=0 %s  0<t<1 %s  t=1 %s" % tuple(
    "(" + ",".join(str(x) for x in w) + ")" for w in (W0, WI, W1)))

# =============================================================================================
# 3.  panel d: rarity at alpha = 1/3
# =============================================================================================
ORDER = ["W", "N", "E", "S"]
# SI Table 3 of Text/ms.tex (alpha = 1/3): rate per draw, and the product of the marginals
TAB4 = {"W": (3.9846e-6, 2.638e-3), "S": (1.344e-2, 2.637e-3),
        "N": (1.007e-4, 1.186e-3), "E": (9.995e-5, 1.185e-3)}
keys = ("effCC", "effALT", "effZ", "rivP", "rivM", "FR_W", "FR_S", "FR_N", "FR_E", "rayP", "rayM")
NDRAW, CNT = 0.0, {k: 0 for k in keys}
SW13 = "/Users/martin/Documents/CL1/PRFull/cannon/sweep13/out_13"
NTASK = 0
for task in range(512):
    for line in open(os.path.join(SW13, "a_1o3_%d.tsv" % task)):
        f = line.rstrip("\n").split("\t")
        if f and f[0] == "ALPHA":
            assert f[1] == "1/3"
            NDRAW += float(f[2])
            NTASK += 1
            for k, v in zip(keys, f[3:]):
                CNT[k] += int(v)
assert NTASK == 512 and NDRAW == 1633280000, (NTASK, NDRAW)
MARG = {"W": ("effCC", "rivP"), "N": ("effALT", "rivP"), "E": ("effALT", "rivM"),
        "S": ("effCC", "rivM")}
OBS = {k: CNT["FR_" + k] / NDRAW for k in ORDER}
IND = {k: CNT[MARG[k][0]] / NDRAW * CNT[MARG[k][1]] / NDRAW for k in ORDER}
for k in ORDER:                                  # the sweep reproduces SI Table 3 to its 4 digits
    assert abs(OBS[k] / TAB4[k][0] - 1) < 6e-4 and abs(IND[k] / TAB4[k][1] - 1) < 6e-4, k
RATE = {k: TAB4[k][0] for k in ORDER}
INDL = {k: TAB4[k][1] for k in ORDER}
print("rarity (SI Table 3):", {k: (RATE[k], INDL[k], round(RATE[k] / INDL[k], 5)) for k in ORDER})

# =============================================================================================
# 4.  the figure  (layout in inches)
# =============================================================================================
FW, FH = 7.2, 8.0
fig = plt.figure(figsize=(FW, FH))
REN = fig.canvas.get_renderer()


def axes_in(x0, y0, w, h, **kw):
    return fig.add_axes([x0 / FW, y0 / FH, w / FW, h / FH], **kw)


def letter(x, y, s):
    fig.text(x / FW, y / FH, s, fontsize=11, fontweight="bold", ha="left", va="baseline",
             color="0.05")


def title(x, y, s, ha="left"):
    fig.text(x / FW, y / FH, s, fontsize=8.2, color=INK, ha=ha, va="baseline")


# ------------------------------------------------------------------------------ panel a
A0, AH = FH - 2.95, 2.95                          # panel a occupies the top 2.95 in
axA = axes_in(0, A0, FW, AH)
axA.set_xlim(0, FW)
axA.set_ylim(0, AH)
axA.axis("off")
letter(0.08, FH - 0.24, "a")
title(0.36, FH - 0.24, "Memory-2, stochastic: the friendly rivals of wedge $W$ form two separated families")

CW, CH = 0.43, 0.325                              # table cell, inches
YT = AH - 1.10                                    # top edge of the tables (axes coordinates)


def table(xc, cells, colour, name, nstr, line2):
    x0 = xc - 2 * CW + 0.10
    # header: swatch, name, count; then the one-line property
    axA.add_patch(Rectangle((x0 - 0.02, AH - 0.535), 0.15, 0.15, fc=colour, ec=DEEP, lw=0.6))
    tn = axA.text(x0 + 0.20, AH - 0.46, name, fontsize=8.0, fontweight="bold", color=INK,
                  ha="left", va="center")
    wn = tn.get_window_extent(renderer=REN).width / fig.dpi          # inches
    axA.text(x0 + 0.20 + wn + 0.10, AH - 0.46, "%d strata" % nstr, fontsize=7.6,
             color="0.30", ha="left", va="center")
    axA.text(x0 - 0.02, AH - 0.68, line2, fontsize=7.4, color=INK, ha="left", va="center")
    # axis names and outcome labels
    axA.text(x0 + 2 * CW, YT + 0.215, "outcome before", fontsize=6.4, color="0.40",
             ha="center", va="bottom")
    axA.text(x0 - 0.40, YT - 2 * CH, "most recent", fontsize=6.4, color="0.40", rotation=90,
             ha="center", va="center")
    for c in range(4):
        axA.text(x0 + (c + 0.5) * CW, YT + 0.045, OUT[c], fontsize=6.8, color="0.25",
                 ha="center", va="bottom")
    for r in range(4):
        yc = YT - (r + 0.5) * CH
        axA.text(x0 - 0.05, yc, OUT[r], fontsize=6.8, color="0.25", ha="right", va="center")
        for c in range(4):
            j = 4 * r + c
            v = cells[j]
            fc = colour if v == "*" else (ONE if v == "1" else "white")
            axA.add_patch(Rectangle((x0 + c * CW, yc - CH / 2), CW, CH, fc=fc, ec="0.62",
                                    lw=0.5, zorder=2))
            dark = (v == "*" and colour == INV)
            axA.text(x0 + (c + 0.5) * CW, yc + 0.030, r"$\ast$" if v == "*" else v,
                     fontsize=(10.5 if v == "*" else 8.2), ha="center", va="center",
                     color=("white" if dark else INK), zorder=4)
            axA.text(x0 + c * CW + CW - 0.035, yc - CH / 2 + 0.035, "%d" % j, fontsize=4.8,
                     ha="right", va="bottom", color=("0.92" if dark else "0.45"), zorder=4)
    for j, col in ((6, HL_SEP), (8, HL_SEP), (9, HL_P9)):
        r, c = divmod(j, 4)
        axA.add_patch(Rectangle((x0 + c * CW + 0.012, YT - (r + 1) * CH + 0.012),
                                CW - 0.024, CH - 0.024, fc="none", ec=col, lw=1.7, zorder=5))
    axA.add_patch(Rectangle((x0, YT - 4 * CH), 4 * CW, 4 * CH, fc="none", ec="0.35", lw=0.7,
                            zorder=3))
    return x0


XL, XR = 1.95, 5.45
table(XL, CELLS_V, VIS, "binary-anchored", len(VISS),
      "contains all 8 binary friendly rivals, dimension %d" % max(HV))
table(XR, CELLS_I, INV, "binary-invisible", len(INVS),
      "contains no binary strategy, dimension %d" % max(HI))

# (no key under the tables: the legend names the framed cells p6, p8 and p9)

# ------------------------------------------------------------------------------ panel b
BY0, BH = 2.95, 1.75
axb = axes_in(0.62, BY0, 2.78, BH)
letter(0.08, BY0 + BH + 0.30, "b")
title(0.36, BY0 + BH + 0.30, "support strata of the friendly rivals in $W$")
dims = np.arange(8)
hv = np.array([HV[d] for d in dims], float)
hi = np.array([HI[d] for d in dims], float)
axb.bar(dims, hv, width=0.72, color=VIS, ec=DEEP, lw=0.55, zorder=3,
        label="binary-anchored (%d)" % len(VISS))
axb.bar(dims, hi, bottom=hv, width=0.72, color=INV, ec=DEEP, lw=0.55, zorder=3,
        label="binary-invisible (%d)" % len(INVS))
for d in dims:
    axb.text(d, hv[d] + hi[d] + 1.6, "%d" % (hv[d] + hi[d]), fontsize=6.6, ha="center",
             va="bottom", color=INK)
axb.set_xticks(dims)
axb.set_xlim(-0.6, 7.6)
axb.set_ylim(0, 118)
axb.spines["left"].set_bounds(0, 100)
axb.set_yticks([0, 20, 40, 60, 80, 100])
axb.set_xlabel("dimension of the stratum", fontsize=7.5, labelpad=2)
axb.set_ylabel("number of strata", fontsize=7.5, labelpad=2)
axb.grid(axis="y", color="0.90", lw=0.5, zorder=0)
for sp in ("top", "right"):
    axb.spines[sp].set_visible(False)
axb.legend(loc="upper right", bbox_to_anchor=(1.03, 1.06), fontsize=6.8, frameon=False,
           handlelength=1.2, handleheight=0.9, borderaxespad=0.2, labelspacing=0.35)

# ------------------------------------------------------------------------------ panel c
axc = axes_in(4.28, BY0, 2.72, BH)
letter(3.78, BY0 + BH + 0.30, "c")
title(4.06, BY0 + BH + 0.30, "a one-probability witness in wedge $N$")
fig.text(4.06 / FW, (BY0 + BH + 0.12) / FH,
         r"$p(t)=(0,0,0,t,\;0,0,0,0,\;0,1,0,1,\;1,0,1,0)$", fontsize=6.9, color="0.30",
         ha="left", va="baseline")
MS = 5.2
YB = 0.60                                          # the friendly-rival bar
axc.plot([0, 1], [YB, YB], color=PURPLE, lw=3.2, solid_capstyle="butt", zorder=3)
axc.plot([0, 1], [YB, YB], "o", ms=MS, mfc="white", mec=PURPLE, mew=1.3, zorder=4,
         clip_on=False)
axc.text(0.5, YB + 0.028, "friendly rival", fontsize=6.8, color=DEEP, ha="center",
         va="bottom")
# interior: w_CD = w_DC = 1/2, w_CC = w_DD = 0 (two colours, alternating dashes)
axc.plot([0, 1], [0.5, 0.5], color=C_ALT, lw=1.8, zorder=3, solid_capstyle="butt")
axc.plot([0, 1], [0, 0], color=C_CC, lw=1.8, zorder=3, solid_capstyle="butt", clip_on=False)
axc.plot([0, 1], [0, 0], color=C_DD, lw=1.8, zorder=3, ls=(0, (3.0, 3.0)), clip_on=False)
for t in (0, 1):
    axc.plot([t], [0.5], "o", ms=MS, mfc="white", mec=C_ALT, mew=1.3, zorder=5, clip_on=False)
    axc.plot([t], [0], "o", ms=MS, mfc="white", mec="0.35", mew=1.3, zorder=5, clip_on=False)
# the endpoint values, attained
for t, Wt in ((0, W0), (1, W1)):
    axc.plot([t], [float(Wt[1])], "o", ms=MS, mfc=C_ALT, mec=C_ALT, zorder=6, clip_on=False)
    if Wt[0] == Wt[3]:
        axc.plot([t], [float(Wt[0])], "o", ms=MS + 1.0, fillstyle="left", mfc=C_CC,
                 markerfacecoloralt=C_DD, mec="none", mew=0, zorder=6, clip_on=False)
    else:
        axc.plot([t], [float(Wt[0])], "o", ms=MS, mfc=C_CC, mec=C_CC, zorder=6, clip_on=False)
        axc.plot([t], [float(Wt[3])], "o", ms=MS, mfc=C_DD, mec=C_DD, zorder=6, clip_on=False)
YT_C = sorted({F(0), F(1, 2)} | set(W0) | set(W1))
axc.set_yticks([float(y) for y in YT_C])
axc.set_yticklabels(["0" if y == 0 else "%d/%d" % (y.numerator, y.denominator) for y in YT_C])
axc.set_xticks([0, 0.25, 0.5, 0.75, 1])
axc.set_xticklabels(["0", "1/4", "1/2", "3/4", "1"])
axc.set_xlim(-0.06, 1.06)
axc.set_ylim(-0.045, 0.70)
axc.set_xlabel("$t = p_3$", fontsize=7.5, labelpad=2)
axc.set_ylabel("self-play frequency", fontsize=7.5, labelpad=2)
for sp in ("top", "right"):
    axc.spines[sp].set_visible(False)
axc.spines["left"].set_bounds(0, 0.5)
axc.grid(axis="y", color="0.92", lw=0.5, zorder=0)
axc.legend(handles=[Line2D([], [], color=C_ALT, lw=1.8, label="$w_{CD}=w_{DC}$"),
                    Line2D([], [], color=C_CC, lw=1.8, label="$w_{CC}$"),
                    Line2D([], [], color=C_DD, lw=1.8, label="$w_{DD}$")],
           loc="center", bbox_to_anchor=(0.5, 0.40), fontsize=6.8, frameon=False, ncol=1,
           handlelength=1.6, labelspacing=0.3)

# ------------------------------------------------------------------------------ panel d
DY0, DH = 0.45, 1.58
axd = axes_in(0.62, DY0, 4.10, DH)
letter(0.08, DY0 + DH + 0.30, "d")
title(0.36, DY0 + DH + 0.30, r"how rare friendly rivals are  ($\alpha=1/3$)")
LOGLO, LOGHI = -7.0, -1.2
RCMAP = LinearSegmentedColormap.from_list("fr", ["#f2edf9", LILAC, PURPLE, DEEP])
RNORM = Normalize(LOGLO, LOGHI)


def rcol(x):
    return RCMAP(RNORM(np.log10(x)))


XMIN, XMAX = 1e-6, 1.0
ys = np.arange(4)
for i, k in enumerate(ORDER):
    axd.barh(i, RATE[k] - XMIN, left=XMIN, height=0.50, color=rcol(RATE[k]), ec=DEEP, lw=0.6,
             zorder=3)
    axd.plot([INDL[k]] * 2, [i - 0.36, i + 0.36], color="0.05", lw=1.8, zorder=6,
             solid_capstyle="butt")
    m, e = RATE[k] / 10 ** np.floor(np.log10(RATE[k])), int(np.floor(np.log10(RATE[k])))
    axd.text(RATE[k] * 1.25, i, r"$%.2f\times10^{%d}$" % (m, e), fontsize=6.8, ha="left",
             va="center", color=INK, zorder=7,
             path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
    ya = i + 0.43
    axd.add_patch(FancyArrowPatch((INDL[k], ya), (RATE[k], ya), arrowstyle="-|>",
                                  mutation_scale=6, lw=0.8, color="0.25", zorder=6,
                                  shrinkA=0, shrinkB=0))
axd.set_xscale("log")
axd.set_xlim(XMIN, XMAX)
axd.set_ylim(3.62, -0.55)
axd.set_yticks(ys)
axd.set_yticklabels(ORDER, fontsize=8.5, fontweight="bold")
axd.tick_params(axis="y", length=0, pad=4)
axd.set_xlabel("friendly rivals per random strategy", fontsize=7.5, labelpad=2)
axd.grid(axis="x", color="0.90", lw=0.5, zorder=0)
for sp in ("top", "right", "left"):
    axd.spines[sp].set_visible(False)
for i, k in enumerate(ORDER):                      # fold-change labels, above the arrows
    fold = RATE[k] / INDL[k]
    lab = (r"$%.3g\times$ below" % (1 / fold)) if fold < 1 else (r"$%.3g\times$ above" % fold)
    if fold < 1:
        xm, yl, ha, va = np.sqrt(RATE[k] * INDL[k]), i + 0.40, "center", "bottom"
    else:
        xm, yl, ha, va = RATE[k] * 1.3, i + 0.43, "left", "center"
    axd.text(xm, yl, lab, fontsize=6.0, ha=ha, va=va, color="0.25", zorder=7,
             path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
axd.legend(handles=[Rectangle((0, 0), 1, 1, fc=rcol(1e-4), ec=DEEP, lw=0.6,
                              label="observed rate"),
                    Line2D([], [], color="0.05", lw=1.8,
                           label="P(efficient) $\\times$ P(competitive)")],
           loc="upper right", bbox_to_anchor=(1.0, 1.08), fontsize=6.6, frameon=False,
           handlelength=1.3, labelspacing=0.35)

# ---------------------------------------------------------------- the locator disk (panel d)
K = 4.0
LIM = 1.175
DD_ = 2.0                                           # disk axes side, inches
axk = axes_in(5.05, DY0 - 0.12, DD_, DD_)
DPU = (2 * LIM / DD_) / 72.0


def phi(p):
    return disklib.phi(p, K)


def ray(dx, dy, n=600):
    t = np.concatenate([[0.0], np.geomspace(1e-3, 1e5, n)])
    return phi(np.c_[1 + t * dx, t * dy])


def rim(a0, a1, n=200):
    a = np.deg2rad(np.linspace(a0, a1, n))
    return np.c_[np.cos(a), np.sin(a)]


UL, DL, URt, DRt = ray(-1, 1), ray(-1, -1), ray(1, 1), ray(1, -1)
WEDGE = {"W": np.vstack([UL, rim(135, 225), DL[::-1]]),
         "N": np.vstack([URt, rim(45, 135), UL[::-1]]),
         "E": np.vstack([DRt, rim(-45, 45), URt[::-1]]),
         "S": np.vstack([DL, rim(225, 315), DRt[::-1]])}
for k, poly in WEDGE.items():
    axk.add_patch(Polygon(poly, closed=True, fc=rcol(RATE[k]), ec="none", zorder=1))
GC, GLW, GA = "white", 0.5, 0.35
axk.plot([-1, 1], [0, 0], color=GC, lw=GLW, alpha=GA, zorder=5)
axk.plot([0, 0], [-1, 1], color=GC, lw=GLW, alpha=GA, zorder=5)
th = np.linspace(0, 2 * np.pi, 2001)
PEL = [pe.withStroke(linewidth=1.4, foreground="white")]
for r in (1, 2, 4, 10):
    rr = float(np.hypot(*phi([[r, 0.0]])[0]))
    axk.plot(rr * np.cos(th), rr * np.sin(th), color=GC, lw=GLW, alpha=GA, zorder=5)
    for sg in (-1, 1):
        axk.plot([sg * rr] * 2, [-0.014, 0.014], color="0.2", lw=0.7, zorder=8)
        axk.plot([-0.014, 0.014], [sg * rr] * 2, color="0.2", lw=0.7, zorder=8)
        lab = ("%g" % r) if sg > 0 else "$-$%g" % r
        axk.text(sg * rr, -0.035, lab, fontsize=4.2, color="0.2", ha="center", va="top",
                 zorder=9, path_effects=PEL)
        axk.text(-0.035, sg * rr, lab, fontsize=4.2, color="0.2", ha="right", va="center",
                 zorder=9, path_effects=PEL)
HALO = [pe.withStroke(linewidth=1.6, foreground="0.35")]
for seg in disklib.conic(1, 1, 1, phi, 2000.0, (-1, 1)):             # switch line u+v = 1
    axk.plot(seg[:, 0], seg[:, 1], color="crimson", lw=0.9, ls=(0, (5, 3)), zorder=7)
for seg in disklib.conic(1, -1, 1, phi, 2000.0, (-1, 1)):            # T = S, i.e. u-v = 1
    axk.plot(seg[:, 0], seg[:, 1], color="white", lw=0.9, ls=(0, (5, 3)), zorder=7,
             path_effects=HALO)
for seg in disklib.conic(1, 1, 0, phi, 2000.0, (1,)):                # donation game
    axk.plot(seg[:, 0], seg[:, 1], color="#12d7e8", lw=0.9, ls=(0, (6, 2.2, 1.2, 2.2)),
             zorder=7)
axk.add_patch(plt.Circle((0, 0), 1.0, fc="none", ec="0.25", lw=0.8, zorder=8))
for k, (u, v) in {"W": (-5.5, -1.6), "N": (0.9, 5.5), "E": (6.0, -1.2), "S": (0.9, -5.0)}.items():
    x, y = phi([[u, v]])[0]
    light = RNORM(np.log10(RATE[k])) > 0.55
    axk.text(x, y, k, fontsize=9, fontweight="bold", ha="center", va="center",
             color=("white" if light else INK), zorder=10)
axk.text(1.03, 0.0, "$u$", fontsize=7.5, ha="left", va="center")
axk.text(0.0, 1.03, "$v$", fontsize=7.5, ha="center", va="bottom")
NMFS = 6.3
for a0, nm in disklib.NAMES:
    _, _, ro = disklib.arc_text(axk, nm, a0, 1.0, NMFS, DPU, disklib.NM_GAP, disklib.NM_TRK,
                                color="0.10")
    assert ro < LIM, (nm, ro)
axk.set_xlim(-LIM, LIM)
axk.set_ylim(-LIM, LIM)
axk.set_aspect("equal")
axk.axis("off")

fig.savefig(os.path.join(HERE, "Fig3.pdf"))
fig.savefig(os.path.join(HERE, "Fig3.png"), dpi=200)
print("wrote Fig3.pdf, Fig3.png")
