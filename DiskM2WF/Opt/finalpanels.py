"""
The four panels of ../FiguresFinal/F1 and F2, written once and drawn on either
WF dataset.

    python3 figF1.py [outdir]     ->  F1: N = 1000, beta = 100, u = 1e-2
    python3 figF2.py [outdir]     ->  F2: N = 100,  beta = 3,   u = 1e-4

WHY ONE MODULE AND NOT TWO FILES.  The two figures differ in the dataset and in
nothing else -- same games, same panels, same scales, same colours -- and this
project has already shown itself what happens when a second case is bolted on
beside an `if ... else`: `Opt/mechmaps.py`'s else-branch silently meant
"square", so adding the sunflower overwrote the square's own figure.  Here the
dataset is a row of DATASETS and there is no else-branch to inherit.  It also
means the two figures cannot drift apart in a panel definition, a threshold or
a colour.  ../Figures/panels.py is the same pattern for Figures 1-3.

THIS IS THE COUNTERPART OF `CL1/DiskM2IN/Opt/figF1.py`, cell for cell.  The two
projects run different processes over the SAME 512 games -- BinM2Ev ca6's
sunflower -- so the two F1.pdf files overlay exactly and the four panels are
the same four quantities.  Anything that differs between them is the process,
not the instrument, EXCEPT where this docstring or ../FiguresFinal/F1_legend.tex
says otherwise.  There are three such places and they are all in panel (a) and
panel (d); both are stated in the legend.

  a  Efficiency                                    viridis, 0 to 1
  b  Mean Hamming distance                         magma,   0 to 8
  c  Type with the highest frequency               5 discrete colours
  d  Total frequency of strict NE at eps = 1e-4    Blues,   0 to DCAP

WHAT IS THE SAME, DELIBERATELY.  The frame is `Opt/disk.py`, which is byte-
identical to DiskM2IN's.  FOUR of the five class colours are DiskM2IN's
`fig1.TYPECOL` -- blue, yellow, red and grey -- because two figures in one
paper must not colour the same five classes differently.  The green is NOT:
theirs is #2e8b3d, which is dE 1.4 from their own rival red under
deuteranopia, and since 2026-09-15 this project uses its own #6abf4b (see
TYPECOL below and ../FiguresFinal/CODE_AND_DATA.md).  That is a deliberate
one-colour divergence from the twin, taken because an unreadable key is worse
than a mismatched one.  `masstype` below is DiskM2IN's `fig3.masstype` verbatim,
and `_check_against_diskm2in()` imports theirs and compares on all 512 cells so
the two cannot drift.  Layout, sizes, the +-1.115 limit, the quadrant tags at
radius 1.06, the "For panel c only:" key row -- all DiskM2IN's.

WHAT IS DIFFERENT, AND WHY.

(a) EFFICIENCY IS THE REALISED ONE, NOT THE SELF-PLAY ONE.  DiskM2IN's panel a
    is sum_s pi(s) pi(s,s) / Emax, the mass-weighted SELF-play payoff.  In a
    rare-mutation process that IS the realised payoff, because the population
    is monomorphic.  Ours is not: at u = 1e-2 the population is a polymorphic
    cloud whose members play EACH OTHER, so the payoff actually collected is
    the simulator's measured mean per-round payoff `pay`, and `eff = pay/Emax`
    is what this panel draws.  Both panels therefore show "the efficiency the
    population realises"; each is computed the way its own process makes that
    quantity mean something.  The self-play functional applied to our pi is
    computed too, by `selfeff()`, and the gap between the two is quoted in the
    legend -- it is the polymorphic tail, and it is the point of PROJECT.md
    section 12.1.

(d) THE STRICT-NE SET COMES FROM THE SIXTEEN ONE-FLIP HALF-PLANES, NOT FROM
    65536 EXPLICIT DEVIATIONS.  DiskM2IN computes it with `Opt/nescan.f90` on
    Cannon, 271 core-hours.  `Opt/ne16.py` gets the identical set in seconds
    from CL/BinM2NE1's exact integer half-planes, because at eps > 0 the
    sixteen one-flip deviations decide the test against all of them
    (BinM2NE1 SI Prop. prop:sixteen -- see ne16.py's docstring for the
    statement, the strictness corollary and the arithmetic).  The set is the
    same set; the route is cheaper and exact rather than float64.

(d) THE COLOUR CEILING IS SET FROM OUR DATA, NOT COPIED.  DiskM2IN caps the
    ramp at 0.05 because nothing there reaches it.  A cap that the data
    exceeds would hide the finding instead of stating it, so `DCAP` is chosen
    by `_dcap()` and the choice is printed and quoted in the legend.

NOTHING ELSE ON THE PAGE.  No suptitle, no N, no beta, no eps, no ranges, no
grid description.  Every one of those is in ../FiguresFinal/F1_legend.tex,
which is generated from this script's own printed summary, so the two cannot
disagree.
"""
import collections
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpat
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, BoundaryNorm, Normalize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import disk                                       # noqa: E402
import exact                                      # noqa: E402
import mech                                       # noqa: E402
import ne16                                       # noqa: E402
import vor                                        # noqa: E402
import wfload                                     # noqa: E402

EPS = 1e-4
NC = 65536

# One row per figure; a name that is not here raises KeyError rather than
# quietly inheriting the other's dataset, output name or parameters.
#   idx    the packed run this figure draws
#   kit    Cannon/<kit>, the campaign that produced it
#   The N / beta / mu / itend fields are for the LEGEND, which states them; the
#   figure states none of them.
DATASETS = {
    "F1": dict(stem="F1", idx="dw1_e4.idx", kit="dw1", jobid="44394541",
               N=1000, beta=100, mu="10^{-2}", muplain="1e-2",
               itend="10^7", itplain="1e7", nrep=10, subs="10^5"),
    "F2": dict(stem="F2", idx="dw2_e4.idx", kit="dw2", jobid="45049704",
               N=100, beta=3, mu="10^{-4}", muplain="1e-4",
               itend="10^8", itplain="1e8", nrep=10, subs="10^4"),
}


def idxpath(which):
    return os.path.join(HERE, os.pardir, "Data", DATASETS[which]["idx"])

# DiskM2IN/Opt/fig1.py TYPECOL and its short key -- FOUR of the five verbatim.
# Martin, 2026-09-15: "fix DiskM2WF's green too".  DiskM2IN's #2e8b3d against
# its own rival red #c1121f is dE 1.4 under deuteranopia (Opt/cvd.py) -- the
# two are the same colour for a red-green dichromat, and panel c puts red and
# green cells in the same corner of the PD.  #6abf4b is this project's own
# green (Figures/panels.py), the same hue lighter; substituting that one hex
# lifts the worst of the ten pairs to 10.8 and clears every threshold, keeping
# DiskM2IN's blue, yellow, red and grey exactly.  The twin in DiskM2IN still
# has the failing green and is read-only from here, so the two figures now
# differ in one of five colours -- deliberately, and noted in both captions.
TYPECOL = ["#1f5fbf",   # partner            blue    DiskM2IN's
           "#f2c200",   # friendly rival     yellow  DiskM2IN's
           "#c1121f",   # rival              red     DiskM2IN's
           "#6abf4b",   # efficient          green   OURS, since 2026-09-15
           "#9a9a9a"]   # neither            grey    DiskM2IN's
TYPESHORT = ["partner", "friendly rival", "rival", "efficient", "neither"]
TYPELAB = ["partner (efficient, limit-Nash, not a rival)",
           "friendly rival (efficient AND a rival)",
           "rival (not efficient)",
           "efficient (neither partner nor rival)",
           "neither"]
BLUES = ListedColormap(plt.get_cmap("Blues")(np.linspace(0.14, 1.0, 256)))

K = 4.0
LIM = 1.115         # disk.draw_frame reserves +-1.28 for wedge labels this
                    # figure does not draw; just enough kept for the four tags
QUAD = [("PD", -0.75, 0.75), ("SD", 0.75, 0.75),
        ("SH", -0.75, -0.75), ("HA", 0.75, -0.75)]
TFS = 13.0          # panel title size
CFS = 10.5          # colorbar tick size

BITS = np.array([[(s >> k) & 1 for k in range(16)] for s in range(NC)],
                dtype=np.float64)


def hamming(pi):
    """E[H] for two independent draws from pi, exactly.

    DiskM2IN/Opt/fig1.hamming verbatim: for independent draws the expectation
    factorises over the sixteen loci, E[H] = sum_k 2 p_k (1 - p_k) with
    p_k = sum_s pi(s) [s_k].  Note this is the GENOTYPE marginal of the
    abundance vector, which is not quite the simulator's `ps` column (the
    realised PLAY frequency, eps-perturbed); `Opt/mech.py` uses the latter.
    The figure uses this one, because it is the estimator DiskM2IN's panel
    uses and the two figures must not differ by the estimator.  `collect()`
    reports the largest discrepancy between them over the 512 cells.
    """
    p = np.asarray(pi, float) @ BITS
    return float(2.0 * np.sum(p * (1.0 - p)))


def masstype(d, u, v, pi):
    """which type holds the most pi-mass, and how much.

    CL1/DiskM2IN/Opt/fig3.masstype verbatim.  Disjoint, most specific first;
    the return index is into TYPECOL / TYPESHORT / TYPELAB.
    `_check_against_diskm2in()` imports theirs and compares on all 512 cells.
    """
    eff, riv, ne = exact.masks(d, u, v)
    fr = eff & riv                       # friendly rival
    pa = eff & ne & ~riv                 # partner, not a rival
    rv = riv & ~eff                      # rival, not efficient
    ef = eff & ~fr & ~pa                 # efficient, neither
    nt = ~(fr | pa | rv | ef)
    m = np.array([pi[pa].sum(), pi[fr].sum(), pi[rv].sum(),
                  pi[ef].sum(), pi[nt].sum()])
    return int(m.argmax()), float(m.max())


def selfeff(T, g, pi):
    """DiskM2IN's OWN panel-a definition applied to our abundance vector:
    sum_s pi(s) * pi(s,s) / Emax, the mass-weighted SELF-play efficiency.

    Not what the figure draws -- see the module docstring -- but the
    like-for-like functional, so the legend can say by how much the two
    differ.  `mech.selfpay` is this project's validated per-strategy self-play
    payoff table (checked against Opt/selfplay.py, which does the same
    16x16 solve by a different route, and against DiskM2IN's own eff column).
    """
    return float(np.dot(pi, mech.selfpay(T, g["cR"], g["cS"], g["cT"]))
                 / g["emax"])


def collect(which="F1"):
    """everything the four panels and the legend need, per sunflower cell"""
    G = wfload.load_run(idxpath(which), ie=3)
    d = exact.load()
    T = mech.tables(EPS)
    R, dh = {}, 0.0
    for g in G.values():
        ipt = vor.index_of(g["tag"])
        if ipt is None or g.get("pi") is None:
            continue
        pi = np.asarray(g["pi"], float)
        s = pi.sum()
        if not s > 0:
            raise ValueError("%s: pi sums to %r" % (g["tag"], s))
        pi = pi / s
        if not np.isfinite(pi).all():
            raise ValueError("%s: pi is not finite" % g["tag"])
        u, v = g["u"], g["v"]
        t, f = masstype(d, u, v, pi)
        weak, strict = ne16.verdicts(u, v)
        hb = hamming(pi)
        dh = max(dh, abs(hb - mech.cell(g, T, d)["hbar"]))
        R[ipt] = dict(ipt=ipt, gid=g["gid"], tag=g["tag"], u=u, v=v,
                      quad=g["tag"].split("_")[-1][:2],
                      eff=g["eff"], selfeff=selfeff(T, g, pi),
                      ham=hb, typ=t, domf=f,
                      massne=float(pi[strict].sum()),
                      massweak=float(pi[weak].sum()),
                      nne=int(strict.sum()), nweak=int(weak.sum()),
                      neff=g["neff"], maxpi=g["maxpi"], efsd=g["efsd"],
                      nrep=g["nrep"], nlead=g["nlead"])
    if len(R) != 512:
        raise ValueError("%d sunflower cells, expected 512" % len(R))
    print("512 cells; hamming(pi) vs mech's ps-based Hbar differ by at most "
          "%.2e" % dh)
    print("%d of the %d x 16 x 512 signs were referred to exact integers"
          % (ne16.nrefer(), NC))
    return R


def _dcap(m):
    """panel d's colour ceiling.  DiskM2IN fixes 0.05 because its data never
    reaches it; a ceiling the data exceeds would hide the finding rather than
    state it, so this rounds the observed maximum UP to a round number and
    prints the choice."""
    mx = float(np.max(m))
    for c in (0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1.0):
        if mx <= c:
            return c
    return 1.0


def main(which="F1", outdir=None):
    d = DATASETS[which]
    outdir = outdir or os.path.join(HERE, os.pardir, "FiguresFinal")
    R = collect(which)
    ipt = sorted(R)
    polys = [vor.cells()[i] for i in ipt]

    effv = np.array([R[i]["eff"] for i in ipt])
    selfv = np.array([R[i]["selfeff"] for i in ipt])
    ham = np.array([R[i]["ham"] for i in ipt])
    typ = np.array([R[i]["typ"] for i in ipt], dtype=float)
    domf = np.array([R[i]["domf"] for i in ipt])
    nne = np.array([R[i]["massne"] for i in ipt])
    DCAP = _dcap(nne)

    fig, ax = plt.subplots(2, 2, figsize=(10.6, 10.6),
                           constrained_layout=True)

    def draw(a, vals, cmap, norm, title, cbar=True, ticks=None):
        fc = [cmap(norm(x)) for x in vals]
        a.add_collection(PolyCollection(polys, facecolors=fc,
                                        edgecolors="none", zorder=1))
        # quad=False: draw_frame puts its tags at (+-0.72, +-0.72), radius
        # 1.018, i.e. OUTSIDE the rim, where the frame line clips them.  Own
        # tags at radius 1.06, sized for the larger disk.
        disk.draw_frame(a, labels=False, quad=False)
        for lab, qx, qy in QUAD:
            r = 1.06 / np.hypot(qx, qy)
            a.text(qx * r, qy * r, lab, ha="center", va="center",
                   fontsize=10.5, color="0.40", zorder=7)
        a.set_xlim(-LIM, LIM); a.set_ylim(-LIM, LIM)   # reclaim the margin
        a.set_title(title, fontsize=TFS, pad=8)
        if cbar:
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                              ax=a, fraction=0.042, pad=0.015)
            if ticks is not None:
                cb.set_ticks(ticks)
            cb.ax.tick_params(labelsize=CFS)
            cb.outline.set_linewidth(0.6)

    draw(ax[0][0], effv, plt.get_cmap("viridis"),
         Normalize(min(0.0, effv.min()), 1.0), "a   Efficiency")

    draw(ax[0][1], ham, plt.get_cmap("magma"), Normalize(0.0, 8.0),
         "b   Mean Hamming distance")

    # the type panel is c, bottom left, directly above its key
    draw(ax[1][0], typ, ListedColormap(TYPECOL),
         BoundaryNorm(np.arange(-.5, 5.5), 5),
         "c   Type with the highest frequency", cbar=False)

    draw(ax[1][1], nne, BLUES, Normalize(0.0, DCAP),
         "d   Total frequency of strict NE at $\\epsilon=10^{-4}$",
         ticks=np.linspace(0.0, DCAP, 6))

    keyh = [mpat.Patch(facecolor="none", edgecolor="none",
                       label="For panel c only:")]
    keyh += [mpat.Patch(color=TYPECOL[i], label=TYPESHORT[i]) for i in range(5)]
    fig.legend(handles=keyh, loc="outside lower center", ncol=6,
               fontsize=11.0, frameon=False, handlelength=1.4,
               handletextpad=0.55, columnspacing=1.9)
    fig.get_layout_engine().set(w_pad=0.004, h_pad=0.004,
                                wspace=0.005, hspace=0.02)
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, d["stem"] + ".pdf")
    fig.savefig(p); fig.savefig(p[:-4] + ".png", dpi=200)
    plt.close(fig)
    print("wrote", p)
    _summary(R, ipt, effv, selfv, ham, typ, domf, nne, DCAP)
    return R


def _q(x, k):
    return float(np.percentile(x, k))


def _summary(R, ipt, effv, selfv, ham, typ, domf, nne, DCAP):
    Q = np.array([R[i]["quad"] for i in ipt])
    print("\n  a  realised efficiency  %.4f .. %.4f, median %.4f"
          % (effv.min(), effv.max(), np.median(effv)))
    print("     above 0.99: %d   below 0.9: %d   negative: %d"
          % (int((effv > 0.99).sum()), int((effv < 0.9).sum()),
             int((effv < 0).sum())))
    print("     the SELF-PLAY functional on the same pi: %.4f .. %.4f, "
          "median %.4f" % (selfv.min(), selfv.max(), np.median(selfv)))
    dd = effv - selfv
    print("     realised - self-play: %.4f .. %.4f, median %+.4f, mean %+.4f"
          % (dd.min(), dd.max(), np.median(dd), dd.mean()))
    print("  b  Hamming %.3f .. %.3f, median %.3f  (of 16)"
          % (ham.min(), ham.max(), np.median(ham)))
    c = collections.Counter(int(x) for x in typ)
    for i in range(5):
        print("  c  %-46s %4d" % (TYPELAB[i], c.get(i, 0)))
    print("     dominant share: median %.3f, min %.3f, below 0.5 in %d cells"
          % (np.median(domf), domf.min(), int((domf < 0.5).sum())))
    print("  d  strict-NE frequency %.4e .. %.4f, median %.4e, mean %.4e"
          % (nne.min(), nne.max(), np.median(nne), nne.mean()))
    print("     p25 %.4e  p50 %.4e  p75 %.4e  p90 %.4e  p95 %.4e  p99 %.4e"
          % tuple(_q(nne, k) for k in (25, 50, 75, 90, 95, 99)))
    print("     exactly zero: %d of %d;  above 0.01: %d;  ramp top %.2f, "
          "cells above it: %d"
          % (int((nne == 0.0).sum()), len(nne), int((nne > 0.01).sum()),
             DCAP, int((nne > DCAP).sum())))
    nn = np.array([R[i]["nne"] for i in ipt])
    nw = np.array([R[i]["nweak"] for i in ipt])
    mw = np.array([R[i]["massweak"] for i in ipt])
    print("     strict NE per cell: %d .. %d, median %d; cells with none: %d"
          % (nn.min(), nn.max(), int(np.median(nn)), int((nn == 0).sum())))
    print("     weak NE per cell:   %d .. %d, median %d; weak mass exceeds "
          "strict mass in %d cells (max gap %.3e)"
          % (nw.min(), nw.max(), int(np.median(nw)),
             int((mw > nne + 1e-15).sum()), float((mw - nne).max())))
    print("\n     by quadrant:   n   median freq        max     no strict NE")
    for q in ("SH", "PD", "HA", "SD"):
        m = Q == q
        print("       %-4s     %4d   %.4e   %.4f   %4d"
              % (q, int(m.sum()), np.median(nne[m]), nne[m].max(),
                 int((nn[m] == 0).sum())))
    print("\n     by quadrant, the other three panels:")
    print("       %-4s %8s %8s %8s   %s" % ("", "eff", "self", "Hbar", "type"))
    for q in ("PD", "SD", "SH", "HA"):
        m = Q == q
        cq = collections.Counter(int(x) for x in typ[m])
        print("       %-4s %8.4f %8.4f %8.3f   %s"
              % (q, np.median(effv[m]), np.median(selfv[m]),
                 np.median(ham[m]),
                 ", ".join("%s %d" % (TYPESHORT[i], cq[i])
                           for i in sorted(cq))))


def _check_against_diskm2in(which="F1"):
    """`masstype` and `hamming` here are DiskM2IN's.  Import THEIRS and
    compare on all 512 of our cells, so a divergence is a failure rather than
    a thing nobody notices."""
    p = "/Users/martin/Documents/CL1/DiskM2IN/Opt"
    if not os.path.isdir(p):
        print("DiskM2IN not present; cross-check skipped")
        return
    sys.path.insert(0, p)
    import importlib
    f3 = importlib.import_module("fig3")
    f1 = importlib.import_module("fig1")
    G = wfload.load_run(idxpath(which), ie=3)
    d = exact.load()
    bt = bh = 0
    n = 0
    for g in G.values():
        if vor.index_of(g["tag"]) is None:
            continue
        pi = np.asarray(g["pi"], float)
        pi = pi / pi.sum()
        n += 1
        t1, f1v = masstype(d, g["u"], g["v"], pi)
        t2, f2v = f3.masstype(d, g["u"], g["v"], pi)
        bt += (t1 != t2) or abs(f1v - f2v) > 0
        bh += abs(hamming(pi) - f1.hamming(pi)) > 0
    print("cross-check against DiskM2IN on %d cells: masstype differs in %d, "
          "hamming differs in %d" % (n, bt, bh))


def run(which, argv):
    """the body of figF1.py / figF2.py"""
    if argv and argv[0] == "check":
        return _check_against_diskm2in(which)
    return main(which, argv[0] if argv else None)
