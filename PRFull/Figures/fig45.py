#!/usr/bin/env python3
"""
fig45.py -- Figures 4 and 5 of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity", and SI Figures 1 and 2.

    python3 fig45.py        ->  Fig4.pdf/.png, Fig5.pdf/.png      (run F1: N = 1000, beta = 100, mu = 1e-2)
                                SIFig1.pdf/.png, SIFig2.pdf/.png  (run F2: N = 100,  beta = 3,   mu = 1e-4)
                                fig45_stats.txt                   (everything printed, and the census tables below)

The layout is that of Figures 4 and 5 of the binary paper (PartnersRivals1/Revised2/ms.tex): two rows, memory one above
memory two, the same 512 Voronoi cells of the Fermat sunflower, eps = 1e-4, 10 replicates per game, second half sampled.
Here both rows are the stochastic strategy spaces, the cubes whose vertices are the binary strategies:

    row 1  "Memory-1, stochastic"  the cube [0,1]^4, PartnersRivalsM1f's run (packs cannon/cw/packed/cw{1,2}_e4.tsv and
           the kernel masses kappa_e4.tsv), read by that project's own Figures/fig4c.continuous() and cwdata.load()
    row 2  "Memory-2, stochastic"  the cube [0,1]^16, PartnersRivalsM2f's run (check/wfc_a333 for F1; check/wfc_f2 for F2,
           with the kernel masses of wfc_a333, the kernel being the same), read by that project's Figures/figwfc.py
Both runs mutate from the three-outcome law with alpha = 1/3 per coordinate (0 w.p. 1/3, 1 w.p. 1/3, uniform otherwise).

Figure 4 / SI Figure 1: efficiency (mean per-round payoff / Emax), the most abundant atom, the most enriched atom (share over
the kernel's mass on the atom).  In the memory-two row's most-enriched panel a dot marks each game at which the ten
replicates leave the most enriched atom undetermined at 95 % (figwfc.undetermined_rep).  In both rows' most-abundant panels
(b, e; MARK_WA) the same dot marks each game at which not all ten replicates name the pooled most abundant atom (the atom of
largest share in each replicate's sampled half, np.argmax over ORDER, against the pooled winner wa): the pooled leader is a
finite-horizon summary there.  The same white dot is used because each panel carries one kind of mark only, and it means in
every panel that the ten replicates do not settle the atom that panel draws; the criterion differs (unanimity in b and e,
the 95 % enrichment test in f) and the legend states each.  The counts per panel are printed after each "drawing" line.
Figure 5 / SI Figure 2: the share of the population that is efficient (1**), stable (*1*), competitive (**1).  Grey would
mark a game at which the kernel puts no mass on the property (M1f's fig5c rule for row 1, figwfc's for row 2); no cell is
grey in either run, so fig5.draw's key "grey: no strategy has it" is hidden (asserted per sheet, printed).
DOT_M1 = True would also dot row 1's most-enriched panel where the same ten-replicate test (kernel masses Poisson in the
1e7 kappa draws) is undetermined; the default False follows the brief (dots on the memory-two row only); the counts for
row 1 are printed either way.

The drawing is the binary paper's fig4.draw and fig5.draw, imported READ-ONLY from PartnersRivals1/FinalFigures with their
output directory redirected here; its property names (efficient, stable, competitive) are kept -- figwfc.py's in-memory
renaming to the older "Nash" / "rival" is undone on import.  PNGs at 200 dpi.  Nothing outside PRFull is written
(no bytecode either).
"""
import glob
import importlib.util
import io
import os
import sys

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FF = "/Users/martin/Documents/CL1/PartnersRivals1/FinalFigures"
M1F = "/Users/martin/Documents/CL1/PartnersRivalsM1f/Figures"
M2F = "/Users/martin/Documents/CL1/PartnersRivalsM2f"
CDIR = {"F1": os.path.join(M2F, "check", "wfc_a333"), "F2": os.path.join(M2F, "check", "wfc_f2")}
NDIR = os.path.join(M2F, "check", "wfc_a333")                    # the alpha-measures: F1's 20 000 (+ 1e6 refinement) draws
REPM1 = {"F1": "cw1_e4.rep.tsv", "F2": "cw2_e4.rep.tsv"}
STEMS = {"F1": ("Fig4", "Fig5"), "F2": ("SIFig1", "SIFig2")}
NLAB = {"F1": "$N = 1000$", "F2": "$N = 100$"}
PARS = {"F1": "N = 1000, beta = 100, mu = 1e-2, eps = 1e-4, 10 x 1e7 generations",
        "F2": "N = 100, beta = 3, mu = 1e-4, eps = 1e-4, 10 x 1e8 generations"}
PNG_DPI = 200
DOT_M1 = False          # True: also dot the memory-one most-enriched panel where the same ten-replicate test is undetermined
MARK_WA = True          # dot the most-abundant panels (b, e) where not all ten replicates name the pooled most abundant atom


class Tee(io.TextIOBase):
    def __init__(self, *s):
        self.s = s

    def write(self, x):
        for s in self.s:
            s.write(x)
        return len(x)

    def flush(self):
        for s in self.s:
            s.flush()


LOG = io.StringIO()
sys.stdout = Tee(sys.__stdout__, LOG)

sys.path.insert(0, FF)
import wfdata                                                      # noqa: E402
import fig4                                                        # noqa: E402
import fig5                                                        # noqa: E402

PROPS0 = [tuple(p) for p in fig5.PROPS]
assert PROPS0 == [("efficient", "1**", ["100", "110", "111"]), ("stable", "*1*", ["010", "110", "011", "111"]),
                  ("competitive", "**1", ["001", "011", "111"])], "PartnersRivals1 changed the properties"
sys.path.insert(1, os.path.join(M2F, "Figures"))
import figwfc                                                      # noqa: E402  the memory-two loaders
fig5.PROPS = [tuple(p) for p in PROPS0]                            # figwfc renames them on import (Nash / rival): undo
# The kernel mass of the friendly rivals, 111, depends only on the wedge.  figwfc carries the rates of the sweep at
# alpha = 0.333 to 1/3 by the census ratio; the sweep at alpha = 1/3 exactly (PRFull cannon/sweep13, array 48928780,
# 1 633 280 000 draws; sweep13/sweep13.json) replaces them.  They differ from figwfc's by +1.0% (W), -0.03% (S),
# +0.57% (N) and +0.11% (E); at N = 1000 no most enriched atom and no determination changes, at N = 100 one undetermined
# game of N (ipt 264) changes from 111 to 110.
FR13 = {"W": 6508 / 1633280000, "S": 21949101 / 1633280000, "N": 164451 / 1633280000, "E": 163245 / 1633280000}
figwfc.FR333.clear()
figwfc.FR333.update(FR13)
spec = importlib.util.spec_from_file_location("m1f_fig4c", os.path.join(M1F, "fig4c.py"))
m1f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m1f)                                       # the memory-one loaders (imports M1f's cwdata)
cwdata = m1f.cwdata
assert fig5.PROPS == PROPS0 and m1f.fig4 is fig4 and figwfc.fig4 is fig4 and figwfc.fig5 is fig5 and m1f.wfdata is wfdata
ORDER = wfdata.ORDER                                               # 000 100 010 001 110 011 111
IX = {c: ORDER.index(c) for c in ORDER}

# --------------------------------------------------------------------------------------------- the games and their regions
GL = [l.split() for l in open(os.path.join(M2F, "cannon", "wfc", "games.dat")) if l.strip() and not l.startswith("#")]
TAG = {int(f[7].split("_")[1][1:]): f[7].split("_")[2] for f in GL}
UVG = {int(f[7].split("_")[1][1:]): (float(f[5]), float(f[6])) for f in GL}
assert sorted(TAG) == list(range(512))
QUAD = np.array([TAG[i][:2] for i in range(512)])
WEDGE = np.array([TAG[i][2] for i in range(512)])
UV = np.array([UVG[i] for i in range(512)])
u_, v_ = UV[:, 0], UV[:, 1]
assert all(figwfc.wedge(*UV[i]) == WEDGE[i] for i in range(512)), "wedge tags disagree with (u, v)"
QCALC = np.where(u_ < 0, np.where(v_ > 0, "PD", "SH"), np.where(v_ > 0, "SD", "HA"))
assert (QCALC == QUAD).all(), "quadrant tags disagree with (u, v)"
REGIONS = [("all", np.ones(512, bool)), ("PD (u<0<v)", QUAD == "PD"), ("SD (u,v>0)", QUAD == "SD"),
           ("SH (u,v<0)", QUAD == "SH"), ("HA (v<0<u)", QUAD == "HA"), ("below switch", np.isin(WEDGE, ["W", "S"])),
           ("above switch", np.isin(WEDGE, ["N", "E"])), ("wedge W", WEDGE == "W"), ("wedge S", WEDGE == "S"),
           ("wedge N", WEDGE == "N"), ("wedge E", WEDGE == "E")]
WEDGES = [(n, m) for n, m in REGIONS if n.startswith("wedge")]


# --------------------------------------------------------------------------------------------- the data
def load_m1(run):
    """row 1: the memory-one cube.  E, most abundant, most enriched exactly as M1f's fig4c.continuous(); property
    shares and their grey rule exactly as M1f's fig5c.continuous_shares() (not imported: its import renames fig5.PROPS)"""
    print("=" * 110)
    print("Memory-1, stochastic, run %s (%s): PartnersRivalsM1f cannon/cw/packed, via its fig4c.continuous()" % (run, PARS[run]))
    E, wa, we = m1f.continuous(False, run=run)
    d = cwdata.load(False, uv_ref=wfdata.load("m1", run)[3], run=run)
    assert np.allclose(d["UV"], UV, rtol=0, atol=1e-12), "the memory-one run is not on the games of games.dat"
    S, K, M = d["S"], d["K"], d["M"]
    PS = np.stack([S[:, [IX[a] for a in atoms]].sum(1) for _, _, atoms in PROPS0], 1)
    PK = np.stack([K[:, [IX[a] for a in atoms]].sum(1) for _, _, atoms in PROPS0], 1)
    PN = np.where((PK == 0) & (PS > 0), cwdata.kfloor(M), PK)
    enr, nfl = cwdata.enrichment(S, K, M)
    # the ten replicates (cw?_e4.rep.tsv): per-replicate shares, in ORDER
    R = {}
    for l in open(os.path.join(cwdata.PACKDIR, REPM1[run])):
        if l.startswith("#") or not l.strip():
            continue
        f = l.split("\t")
        s = np.array([float(x) for x in f[6:14]])
        assert s[7] == 0.0
        R.setdefault(int(f[0]), []).append(s[:7])
    R = np.array([R[i] for i in range(512)])
    assert R.shape == (512, 10, 7)
    mean_ok = np.allclose(R.mean(1), S, rtol=1e-6, atol=1e-9)
    print("   replicate file %s: 512 x 10 lines; the mean of the replicate shares %s the pooled share"
          % (REPM1[run], "equals" if mean_ok else "DIFFERS FROM"))
    return dict(E=E, wa=wa, we=we, S=S, K=K, M=M, PS=PS, PN=PN, enr=enr, R=R, nfl=nfl, meas=K)


def load_m2(run):
    """row 2: the memory-two cube, exactly as figwfc.main() reads it"""
    print("=" * 110)
    print("Memory-2, stochastic, run %s (%s): PartnersRivalsM2f %s, alpha-measures from %s, via its figwfc.py"
          % (run, PARS[run], os.path.relpath(CDIR[run], M2F), os.path.relpath(NDIR, M2F)))
    S, M, E, nr, mdraw, UVc, K, H = figwfc.load_continuum(CDIR[run], nd=NDIR)
    assert (nr == 10).all(), "replicates per game: %s" % np.bincount(nr)
    assert (mdraw > 0).all()
    assert np.allclose(UVc, UV, rtol=0, atol=1e-9), "the memory-two run is not on the games of games.dat"
    wa, we = figwfc.winners_c(S, M, E, "Memory-2, stochastic %s (alpha = 1/3)" % run, K)
    und_p = figwfc.undetermined(S, M, K, H, mdraw)
    R = figwfc.load_reps(CDIR[run])
    und = figwfc.undetermined_rep(R, M, H, mdraw, S=S)
    print("   most enriched NOT determined at 95%% (ten-replicate t interval) at %d games (Poisson in the pooled samples "
          "would say %d)" % (len(und), len(und_p)))
    print("   undetermined games (ipt: winner vs contenders): " + "; ".join(
        "%d %s vs %s" % (i, ORDER[w], "/".join(ORDER[j] for j in r)) for i, w, r in und))
    print("   neutral draws per game: %s" % ", ".join("%d x %d" % (c, v) for v, c in zip(*np.unique(mdraw, return_counts=True))))
    PS = np.stack([S[:, [IX[a] for a in atoms]].sum(1) for _, _, atoms in PROPS0], 1)
    PN = np.stack([M[:, [IX[a] for a in atoms]].sum(1) for _, _, atoms in PROPS0], 1)
    ok = M > 0
    enr = np.where(ok, S / np.where(ok, M, 1.0), -np.inf)
    assert (enr.argmax(1) == we).all()
    dotted = {i for i, _, _ in und}
    hz = [(i, ORDER[j]) for i, j in zip(*np.nonzero((S > 0) & ~ok))]
    if hz:
        print("   held by the population but of zero measured alpha-measure (excluded from the winner, as figwfc does; the "
              "game is then undetermined): " + ", ".join("ipt %d atom %s share %.2g%s" % (i, a, S[i, IX[a]], "" if i in dotted
                                                                                        else " NOT DOTTED") for i, a in hz))
    return dict(E=E, wa=wa, we=we, S=S, M=M, PS=PS, PN=PN, enr=enr, R=R, und=und, und_p=und_p, meas=M, H=H, mdraw=mdraw,
                dotted=True)


def load_binary(space, run):
    """the binary run of the reference paper, for comparison only (not drawn): wfdata.load, fig4.winners"""
    S, N, E, UVb, D = wfdata.load(space, run)
    assert np.allclose(UVb, UV, rtol=0, atol=1e-9)
    wa, we = fig4.winners(S, N, E, D, "binary %s %s (reference paper, uniform mutants)" % (space, run))
    PS = np.stack([S[:, [IX[a] for a in atoms]].sum(1) for _, _, atoms in PROPS0], 1)
    PN = np.stack([N[:, [IX[a] for a in atoms]].sum(1) for _, _, atoms in PROPS0], 1)
    return dict(E=E, wa=wa, we=we, S=S, PS=PS, PN=PN, meas=N)


# --------------------------------------------------------------------------------------------- the census tables
def census(name, d, extra=True):
    E, wa, we, S, PS, PN = d["E"], d["wa"], d["we"], d["S"], d["PS"], d["PN"]
    print("-" * 110)
    print("%s" % name)
    print("  efficiency (mean per-round payoff / Emax)")
    print("    %-14s %5s %8s %8s %8s %6s %6s" % ("region", "games", "min", "median", "max", "<0.9", "<0"))
    for rn, m in REGIONS:
        e = E[m]
        print("    %-14s %5d %8.4f %8.4f %8.4f %6d %6d" % (rn, m.sum(), e.min(), np.median(e), e.max(), (e < 0.9).sum(), (e < 0).sum()))
    for lab, w in (("most abundant atom", wa), ("most enriched atom", we)):
        print("  %s: number of games" % lab)
        print("    %-14s %5s " % ("region", "games") + " ".join("%6s" % c for c in ORDER))
        for rn, m in REGIONS:
            print("    %-14s %5d " % (rn, m.sum()) + " ".join("%6d" % (w[m] == i).sum() for i in range(7)))
    print("  share of the population with each property: mean / median / games above 1/2 / games where no strategy has it")
    print("    %-14s %5s " % ("region", "games") + "".join("%34s" % ("%s %s" % (c, n)) for n, c, _ in PROPS0))
    for rn, m in REGIONS:
        cells = []
        for j in range(3):
            v = PS[m, j]
            cells.append("%8.4f %8.4f %6d %6d  " % (v.mean(), np.median(v), (v > 0.5).sum(), (PN[m, j] == 0).sum()))
        print("    %-14s %5d " % (rn, m.sum()) + "".join("%34s" % c for c in cells))
    print("  share of each atom: median per region (mean in brackets)")
    print("    %-14s %5s " % ("region", "games") + " ".join("%19s" % c for c in ORDER))
    for rn, m in REGIONS:
        print("    %-14s %5d " % (rn, m.sum()) + " ".join("%9.3g (%7.3g)" % (np.median(S[m, i]), S[m, i].mean()) for i in range(7)))
    s111 = S[:, IX["111"]]
    print("  friendly rivals (111): median share per wedge: " + ", ".join(
        "%s %.4g" % (rn.split()[1], np.median(s111[m])) for rn, m in WEDGES) + "; max per wedge: " + ", ".join(
        "%s %.4g" % (rn.split()[1], s111[m].max()) for rn, m in WEDGES))
    part = S[:, IX["110"]] + S[:, IX["111"]]
    print("  partners (110 + 111): median share per region: " + ", ".join("%s %.4f" % (rn, np.median(part[m])) for rn, m in REGIONS))
    print("  partners lead (most abundant atom 110 or 111) at: " + ", ".join(
        "%s %d" % (rn, np.isin(wa[m], [IX["110"], IX["111"]]).sum()) for rn, m in REGIONS))
    print("  an efficient atom (100, 110, 111) leads at: " + ", ".join(
        "%s %d" % (rn, np.isin(wa[m], [IX["100"], IX["110"], IX["111"]]).sum()) for rn, m in REGIONS))
    print("  most enriched != most abundant at %d games" % (we != wa).sum())
    if "R" in d:
        R = d["R"]
        lead_r = R.argmax(2)                                       # (512, 10)
        agree = (lead_r == wa[:, None]).all(1)
        print("  replicate agreement on the most abundant atom: all ten replicates name the pooled leader at %d of 512 games;"
              " the others by region: %s" % (agree.sum(), ", ".join("%s %d" % (rn, (~agree[m]).sum()) for rn, m in REGIONS[1:])))
        nmaj = np.array([(lead_r[i] == wa[i]).sum() for i in range(512)])
        print("    replicates naming the pooled leader, games by count: %s"
              % ", ".join("%d/10: %d" % (k, (nmaj == k).sum()) for k in range(11) if (nmaj == k).any()))
        meas = d["meas"]
        ok = meas > 0
        enr_r = np.where(ok[:, None, :], R / np.where(ok, meas, 1.0)[:, None, :], -np.inf)
        agree_e = (enr_r.argmax(2) == we[:, None]).all(1)
        print("  replicate agreement on the most enriched atom (replicate share / the same kernel mass): all ten at %d of 512"
              % agree_e.sum())
    if "und" in d:
        ud = np.zeros(512, bool)
        ud[[i for i, _, _ in d["und"]]] = True
        print("  most enriched undetermined at 95%% (ten-replicate test%s): %d games; by region: %s"
              % ("; dots in the figure" if d.get("dotted") else "; NOT dotted in the figure", ud.sum(),
                 ", ".join("%s %d" % (rn, ud[m].sum()) for rn, m in REGIONS[1:])))
        print("    determined games, most enriched atom by region:")
        for rn, m in REGIONS:
            mm = m & ~ud
            print("      %-14s %5d " % (rn, mm.sum()) + " ".join("%s %d" % (ORDER[i], (we[mm] == i).sum()) for i in range(7)
                                                             if (we[mm] == i).any()))


def m1_determinacy(d, run):
    """not drawn: figwfc's ten-replicate test applied to the memory-one row (kernel masses from M = 1e7 draws, Poisson)"""
    K, M, R, S = d["K"], d["M"], d["R"], d["S"]
    H = np.round(K * M)
    assert np.allclose(H, K * M, atol=1e-6), "kappa is not a count over M"
    und = figwfc.undetermined_rep(R, K, H, np.full(512, float(M)), S=None)
    ok = K > 0
    w_int = np.where(ok, R.mean(1) / np.where(ok, K, 1.0), -np.inf).argmax(1)
    print("   figwfc's ten-replicate test applied to the memory-one row (kernel masses Poisson in %d draws): most "
          "enriched undetermined at 95%% at %d games%s" % (M, len(und), (": " + "; ".join(
              "%d %s vs %s" % (i, ORDER[w], "/".join(ORDER[j] for j in r)) for i, w, r in und)) if und else ""))
    print("    (its winner, from the replicate mean over atoms of positive kappa, differs from fig4c's at %d games)"
          % (w_int != d["we"]).sum())
    return und


# --------------------------------------------------------------------------------------------- drawing
import matplotlib.figure                                          # noqa: E402  (fig4 has set the Agg backend)
from matplotlib.path import Path                                   # noqa: E402

fig4.HERE = fig5.HERE = HERE                                       # in memory only: fig4.draw and fig5.draw write here
_save0 = matplotlib.figure.Figure.savefig
DOTS = {"xy": None, "xy1": None, "xa1": None, "xa2": None, "hide_grey": False}


def _dot(ax, xy):
    """the one mark of the Figure-4 sheets: a small white dot with a dark rim, legible on every atom colour"""
    ax.scatter(xy[:, 0], xy[:, 1], s=9.0, c="white", edgecolors="0.1", linewidths=0.45, zorder=5)


def _save(fig, fname, *a, **k):
    """PNG at 200 dpi; on a Figure-4 sheet, the undetermined dots on the memory-two most-enriched panel (axes 5) and the
    replicate-disagreement dots on the two most-abundant panels (axes 1 and 4)"""
    if DOTS["xy"] is not None and not getattr(fig, "_dotted", False):
        P = [q.get_position() for q in fig.axes[:6]]
        x_lo, x_hi = min(q.x0 for q in P), max(q.x0 for q in P)
        y_lo, y_hi = min(q.y0 for q in P), max(q.y0 for q in P)
        assert P[5].x0 == x_hi and P[5].y0 == y_lo, "axes 5 is not the bottom-right panel"
        _dot(fig.axes[5], DOTS["xy"])
        if DOTS["xy1"] is not None and len(DOTS["xy1"]):
            assert P[2].x0 == x_hi and P[2].y0 == y_hi, "axes 2 is not the top-right panel"
            _dot(fig.axes[2], DOTS["xy1"])
        for j, key, y_ in ((1, "xa1", y_hi), (4, "xa2", y_lo)):    # the most-abundant panels, middle column
            if DOTS[key] is not None and len(DOTS[key]):
                assert x_lo < P[j].x0 < x_hi and P[j].y0 == y_, "axes %d is not a middle-column panel of its row" % j
                _dot(fig.axes[j], DOTS[key])
        fig._dotted = True
    if DOTS["hide_grey"]:                                          # the grey key of fig5.draw, when no cell is grey
        hit = [t for a_ in fig.axes for t in a_.texts if t.get_text() == "grey: no strategy has it"]
        assert len(hit) == 1, "fig5.draw's grey key not found"
        hit[0].set_visible(False)
    if str(fname).endswith(".png"):
        k["dpi"] = PNG_DPI
    return _save0(fig, fname, *a, **k)


def draw(run, d1, d2):
    L1, L2 = "Memory-1, stochastic, " + NLAB[run], "Memory-2, stochastic, " + NLAB[run]
    s4, s5 = STEMS[run]
    # the dots: games whose most enriched atom the ten memory-two replicates leave undetermined
    pol, keep = wfdata.cells()
    for j, ipt in enumerate(keep):
        site = UV[ipt] / np.sqrt(16.0 + UV[ipt] @ UV[ipt])
        assert Path(np.asarray(pol[j])).contains_point(site), "site %d is not in its Voronoi cell" % ipt
    und = [i for i, _, _ in d2["und"]]
    assert all(i in keep for i in und), "an undetermined game has no drawn cell"
    DOTS["xy"] = np.array([UV[i] / np.sqrt(16.0 + UV[i] @ UV[i]) for i in und]).reshape(-1, 2)
    if DOT_M1:
        u1 = [i for i, _, _ in d1["und"]]
        assert all(i in keep for i in u1)
        DOTS["xy1"] = np.array([UV[i] / np.sqrt(16.0 + UV[i] @ UV[i]) for i in u1]).reshape(-1, 2)
    # the replicate-disagreement dots: games at which not all ten replicates name the pooled most abundant atom
    ma = []
    for key, d in (("xa1", d1), ("xa2", d2)):
        lead_r = d["R"].argmax(2)                                  # (512, 10), in ORDER, as census() and SI Data 1
        dis = np.nonzero(~(lead_r == d["wa"][:, None]).all(1))[0]
        assert all(i in keep for i in dis), "a game with replicate disagreement has no drawn cell"
        ma.append(dis)
        if MARK_WA:
            DOTS[key] = np.array([UV[i] / np.sqrt(16.0 + UV[i] @ UV[i]) for i in dis]).reshape(-1, 2)
    matplotlib.figure.Figure.savefig = _save
    try:
        print("drawing %s (%d dots)" % (s4, len(und)))
        print("   most-abundant panels: not all ten replicates name the pooled most abundant atom at %d games (memory one, b) "
              "and %d (memory two, e) -- %s" % (len(ma[0]), len(ma[1]), "dotted" if MARK_WA else "NOT dotted (MARK_WA False)"))
        fig4.draw([(L1, d1["E"], d1["wa"], d1["we"]), (L2, d2["E"], d2["wa"], d2["we"])], s4)
        DOTS["xy"] = DOTS["xy1"] = DOTS["xa1"] = DOTS["xa2"] = None
        ngrey = [int((d["PN"][keep] == 0).sum()) for d in (d1, d2)]
        DOTS["hide_grey"] = ngrey == [0, 0]
        print("grey cells in %s: memory one %d, memory two %d%s" % (s5, ngrey[0], ngrey[1],
                                                                    " -- the grey key is hidden" if DOTS["hide_grey"] else ""))
        K1, K2 = ("m1s", run), ("m2s", run)
        lab = dict(wfdata.LABEL)
        lab.update({K1: L1, K2: L2})
        wfdata.LABEL = lab                                         # in memory only
        base = fig5.shares
        fig5.shares = lambda space, which: ((d1["PS"], d1["PN"]) if (space, which) == K1 else
                                            (d2["PS"], d2["PN"]) if (space, which) == K2 else base(space, which))
        print("drawing %s" % s5)
        try:
            fig5.draw([K1, K2], s5)
        finally:
            fig5.shares = base
    finally:
        matplotlib.figure.Figure.savefig = _save0
        DOTS.update(xy=None, xy1=None, xa1=None, xa2=None, hide_grey=False)
    print("(the files above are PRFull/Figures/%s.pdf/.png and %s.pdf/.png, PNG at %d dpi -- not FinalFigures)" % (s4, s5, PNG_DPI))


def main():
    print("fig45.py -- Figures 4, 5 (run F1) and SI Figures 1, 2 (run F2) of PRFull; every number below is computed from the data")
    print()
    print("THE 512 GAMES (tags of PartnersRivalsM2f/cannon/wfc/games.dat, checked against the signs of u, v, u+v-1, u-v-1):")
    print("   quadrants: PD (u<0<v) %d, SD (u,v>0) %d, SH (u,v<0) %d, HA (v<0<u) %d"
          % tuple((QUAD == q).sum() for q in ("PD", "SD", "SH", "HA")))
    print("   wedges: W %d, S %d, N %d, E %d; above the switch line u+v>1: %d; below: %d"
          % tuple([(WEDGE == w).sum() for w in "WSNE"] + [np.isin(WEDGE, ["N", "E"]).sum(), np.isin(WEDGE, ["W", "S"]).sum()]))
    print("   PD games by wedge: " + ", ".join("%s %d" % (w, ((QUAD == "PD") & (WEDGE == w)).sum()) for w in "WSNE"))
    print("   SD games by wedge: " + ", ".join("%s %d" % (w, ((QUAD == "SD") & (WEDGE == w)).sum()) for w in "WSNE"))
    print("   SH games by wedge: " + ", ".join("%s %d" % (w, ((QUAD == "SH") & (WEDGE == w)).sum()) for w in "WSNE"))
    print("   HA games by wedge: " + ", ".join("%s %d" % (w, ((QUAD == "HA") & (WEDGE == w)).sum()) for w in "WSNE"))
    print("   donation-game line u+v=0 (v>0): %d games within |u+v| < 0.05; games exactly on u+v=1: %d, on u-v=1: %d, on an axis: %d"
          % (((abs(u_ + v_) < 0.05) & (v_ > 0)).sum(), (u_ + v_ == 1).sum(), (u_ - v_ == 1).sum(), ((u_ == 0) | (v_ == 0)).sum()))
    pol, keep = wfdata.cells()
    print("   drawn cells (Voronoi polygons with >= 3 vertices): %d of 512" % len(keep))
    print()
    out = {}
    for run in ("F1", "F2"):
        d1, d2 = load_m1(run), load_m2(run)
        d1["und"] = m1_determinacy(d1, run)
        d1["dotted"] = DOT_M1
        out[run] = (d1, d2)
        print("=" * 110)
        draw(run, d1, d2)
    for run in ("F1", "F2"):
        d1, d2 = out[run]
        print()
        print("#" * 110)
        print("CENSUS, run %s (%s)" % (run, PARS[run]))
        census("Memory-1, stochastic [0,1]^4, run %s" % run, d1)
        census("Memory-2, stochastic [0,1]^16, run %s" % run, d2)
    print()
    print("#" * 110)
    print("FOR COMPARISON ONLY (not drawn): the binary runs of the reference paper, the same tables")
    for run in ("F1", "F2"):
        for space, nm in (("m1", "Memory-1, binary (16)"), ("m2", "Memory-2, binary (65536)")):
            print("=" * 110)
            census("%s, run %s" % (nm, run), load_binary(space, run))
    sys.stdout.flush()
    with open(os.path.join(HERE, "fig45_stats.txt"), "w") as fh:
        fh.write(LOG.getvalue())
    sys.__stdout__.write("wrote %s\n" % os.path.join(HERE, "fig45_stats.txt"))


if __name__ == "__main__":
    main()
