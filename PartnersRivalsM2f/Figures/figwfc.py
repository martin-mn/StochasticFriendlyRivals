#!/usr/bin/env python3
"""
figwfc.py -- Figures 4 and 5 of PartnersRivals1/Revised2 for the FULL memory-two cube [0,1]^16.

    python3 figwfc.py <dir with r_1o3_*.txt and n_1o3_*.txt>      ->  Fig4c.pdf/.png, Fig5c.pdf/.png (here)

Three rows per sheet, on the same 512 Voronoi cells of the Fermat sunflower, all N = 1000, beta = 100, mu = 1e-2,
eps = 1e-4, 10 x 1e7 generations -- the layout of PartnersRivalsM1f's Figures 5 and 6, which set the memory-one cube
between the 16 and the 65536 binary strategies; here the memory-two cube sits beneath the two spaces it contains:
    row 1   the memory-one cube [0,1]^4: PartnersRivalsM1f's run, computed by that project's own fig4c.continuous()
    row 2   the 65536 binary memory-two strategies: Revised2's run F1, read through PartnersRivals1/FinalFigures/wfdata.py
            exactly as Revised2's fig4.py and fig5.py read it
    row 3   the memory-two cube [0,1]^16, the same process (cannon/wfc/wfc.c), mutants from the three-outcome law at 1/3
Uniform mutation over the binary genomes is that law at alpha = 1/2, so row 1 is the alpha = 1/2 member of the family
and row 2 the alpha = 1/3 member; cannon/wfc's probe checks that wfc.c at alpha = 1/2 reproduces row 1.

Definitions are Revised2's, carried to the continuum:
    efficiency           mean per-round payoff over the second half of each replicate / Emax (exact, not sampled)
    atom share           fraction of individual-generations of the second half held by strategies of the atom;
                         binary: the full histogram; continuum: 4000 uniform samples per replicate, 40 000 per game (111 in W: exact)
    most abundant atom   argmax of the share over the atoms that have strategies at the game
    most enriched atom   argmax of share / (share of the atom in the strategy space).  Binary: |a|/65536, the
                         mutation supply of uniform mutation.  Continuum: the alpha-measure of the atom at the game,
                         the mutation supply of the alpha-law -- from 20 000 draws per game (wfc neutral), except
                         for 111, whose measure is the friendly-rival rate of Table 4 of Text/ms.tex (1.633e9 draws
                         at alpha = 0.333, carried to 1/3 by the census ratio; 111 depends only on the wedge, every
                         friendly rival being Nash).  In wedge W the 111 SHARE is exact too (wfc.c counts every
                         friendly rival's individual-generations), because there 111 is too rare to sample
    property shares      sums of atom shares, as fig5.py
The drawing is Revised2's: fig4.draw and fig5.draw, imported read-only, with their output directory redirected here.
Nothing under PartnersRivals1/ or PartnersRivalsM1f/ is written.
"""
import glob
import os
import sys

sys.dont_write_bytecode = True                                    # import Revised2's modules without writing under PartnersRivals1
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FF = "/Users/martin/Documents/CL1/PartnersRivals1/FinalFigures"
sys.path.insert(0, FF)
import wfdata                                                     # noqa: E402
import fig4                                                       # noqa: E402
import fig5                                                       # noqa: E402

ORDER = wfdata.ORDER                                              # 000 100 010 001 110 011 111
# the column NAMES are this paper's vocabulary (Text/ms.tex: efficient, Nash, rival), as PartnersRivalsM1f's fig5c.py does,
# not PartnersRivals1's current "stable" / "competitive" (renamed there 2026-09-25); the codes and atoms are PR1's, asserted
NAMES = {"1**": "efficient", "*1*": "Nash", "**1": "rival"}
assert [(c, a) for _, c, a in fig5.PROPS] == [("1**", ["100", "110", "111"]), ("*1*", ["010", "110", "011", "111"]),
                                             ("**1", ["001", "011", "111"])], "PartnersRivals1 changed the properties"
fig5.PROPS = [(NAMES[c], c, a) for _, c, a in fig5.PROPS]         # in memory only: fig5.draw reads the names from here
# Table 4 of Text/ms.tex measured the friendly-rival rate at alpha = 0.333 (Cannon 45831015, 1.633e9 draws); this run is at
# alpha = 1/3 exactly, so each rate is carried to 1/3 by the exhaustive stratum census's own ratio (check/st_eff.txt:
# predicted rate at 1/3 over predicted rate at 0.333), a shift of 0.9 % in W, 0.1 % in S, 0.7 % in N and E.
FR0333 = {"W": 3.9111e-06, "S": 1.3426e-02, "N": 9.9431e-05, "E": 9.9162e-05}
TO13 = {"W": 1.00880, "S": 1.00120, "N": 1.00686, "E": 1.00686}
FR333 = {w: FR0333[w] * TO13[w] for w in FR0333}
BITS = ["%d%d%d" % (k >> 2, (k >> 1) & 1, k & 1) for k in range(8)]         # wfc's atom index -> code


def wedge(u, v):
    below, tgs = (u + v < 1.0), (u - v < 1.0)
    return ("W" if tgs else "S") if below else ("N" if tgs else "E")


def load_continuum(d, at="1o3", nd=None):
    """per ipt 0..511: shares S (512,7) in ORDER, alpha-measure M (512,7), efficiency E (512,), replicate counts,
    and K (512,7), the number of samples behind each share (inf where the share is exact: 111 in wedge W)"""
    S = np.zeros((512, 7)); cnt = np.zeros(512); E = np.zeros(512); nr = np.zeros(512, int); UV = np.zeros((512, 2))
    X111 = np.zeros(512); T111 = np.zeros(512); inW = np.zeros(512, bool)
    for fn in glob.glob(os.path.join(d, "r_%s_*.txt" % at)):
        for l in open(fn):
            f = l.split()
            if not f or f[0] != "REP":
                continue
            ipt = int(f[2]); ns = float(f[10]); a = [float(x) for x in f[12:20]]
            assert a[5] == 0, "atom 101 occupied at ipt %d" % ipt
            UV[ipt] = float(f[3]), float(f[4])
            for k in range(8):
                if k != 5:
                    S[ipt, ORDER.index(BITS[k])] += a[k]
            cnt[ipt] += ns; E[ipt] += float(f[7]); nr[ipt] += 1
            if int(f[26]):                                        # wedge W: 111 counted exactly, not sampled
                inW[ipt] = True; X111[ipt] += float(f[27]); T111[ipt] += float(f[30])
    # the alpha-measure: the main run's 20 000 draws per game (n_), pooled with any refinement draws (nr_, cannon/wfc
    # run_neu: a further 1e6 at the games whose most-enriched atom those 20 000 left undetermined)
    C = np.zeros((512, 7)); mdraw = np.zeros(512); UVn = np.zeros((512, 2))
    nd = nd or d                                                  # the kernel's mass depends on the game and the kernel only,
    for fn in glob.glob(os.path.join(nd, "n_%s_*.txt" % at)) + glob.glob(os.path.join(nd, "nr_%s_*.txt" % at)):  # not on N
        f = open(fn).read().split()
        assert f[0] == "NEU"
        ipt = int(f[2]); m = float(f[7]); c = [float(x) for x in f[8:16]]
        assert c[5] == 0
        C[ipt] += [c[BITS.index(code)] for code in ORDER]
        mdraw[ipt] += m; UVn[ipt] = float(f[3]), float(f[4])
    M = C / np.where(mdraw > 0, mdraw, np.nan)[:, None]
    for ipt in range(512):
        if mdraw[ipt] > 0:
            M[ipt, ORDER.index("111")] = FR333[wedge(*UVn[ipt])]
    H = C.copy(); H[:, ORDER.index("111")] = np.inf                # hits behind each measure; 111's is Table 4's 1.6e9 draws
    K = S.copy()                                                  # raw sample counts, before normalising
    S = S / cnt[:, None]
    i111 = ORDER.index("111")
    S[inW, i111] = X111[inW] / T111[inW]                          # the exact share replaces the sampled one
    K[inW, i111] = np.inf
    print("wedge W: %d games, 111 counted exactly; its share median %.3g, max %.3g"
          % (inW.sum(), np.median(S[inW, i111]), S[inW, i111].max()))
    return S, M, E / nr, nr, mdraw, UV, K, H


def undetermined(S, M, K, H, mdraw, nsamp=40000.0):
    """cells whose most-enriched atom is not determined at 95 %: some other atom held by the population has an upper
    enrichment bound at or above the winner's lower bound.  Share bounds from the sample count (exact where K = inf),
    measure bounds from the hit count (exact where H = inf), both Poisson."""
    from scipy.stats import chi2
    lo = lambda k: np.where(k > 0, chi2.ppf(0.025, 2 * np.maximum(k, 1e-12)) / 2, 0.0)
    hi = lambda k: chi2.ppf(0.975, 2 * (k + 1)) / 2
    out = []
    for i in range(512):
        k, h, m = K[i], H[i], mdraw[i]
        fk, fh = np.isinf(k), np.isinf(h)
        s_lo = np.where(fk, S[i], lo(np.where(fk, 0, k)) / nsamp); s_hi = np.where(fk, S[i], hi(np.where(fk, 0, k)) / nsamp)
        m_lo = np.where(fh, M[i], lo(np.where(fh, 0, h)) / m); m_hi = np.where(fh, M[i], hi(np.where(fh, 0, h)) / m)
        held = S[i] > 0
        e_lo = np.where(held, s_lo / np.maximum(m_hi, 1e-300), 0.0)
        e_hi = np.where(held, np.where(m_lo > 0, s_hi / np.maximum(m_lo, 1e-300), np.inf), 0.0)
        ok = M[i] > 0; enr = np.where(ok, S[i] / np.where(ok, M[i], 1.0), -np.inf); w = int(np.argmax(enr))
        riv = [j for j in range(7) if j != w and held[j] and e_hi[j] >= e_lo[w]]
        if riv:
            out.append((i, w, riv))
    return out


def load_reps(d, at="1o3"):
    """per-replicate shares R (512, 10, 7) in ORDER, 111 exact in wedge W, as load_continuum pools them"""
    R = {}
    for fn in glob.glob(os.path.join(d, "r_%s_*.txt" % at)):
        for l in open(fn):
            f = l.split()
            if not f or f[0] != "REP":
                continue
            ipt, ns, a = int(f[2]), float(f[10]), [float(x) for x in f[12:20]]
            s = np.zeros(7)
            for k in range(8):
                if k != 5:
                    s[ORDER.index(BITS[k])] = a[k] / ns
            if int(f[26]):
                s[ORDER.index("111")] = float(f[27]) / float(f[30])
            R.setdefault(ipt, []).append(s)
    return np.array([R[i] for i in range(512)])


def undetermined_rep(R, M, H, mdraw, S=None):
    """as undetermined(), but with the share bounds from the ten replicates, the process's independent unit: a t interval
    on the mean of the per-replicate shares (9 degrees of freedom, floored at 0).  The samples within a replicate are
    not independent -- the replicate-to-replicate variance of an atom's count is a median 21 times the binomial one at
    N = 100 and 4.4 times at N = 1000 -- so Poisson bounds on the pooled sample count are too narrow.  Measure bounds
    stay Poisson in the hit count: the kernel draws are independent."""
    from scipy.stats import chi2, t as tdist
    lo = lambda k: np.where(k > 0, chi2.ppf(0.025, 2 * np.maximum(k, 1e-12)) / 2, 0.0)
    hi = lambda k: chi2.ppf(0.975, 2 * (k + 1)) / 2
    n = R.shape[1]; q = tdist.ppf(0.975, n - 1)
    mean = R.mean(1); half = q * R.std(1, ddof=1) / np.sqrt(n)
    if S is not None:
        assert np.allclose(mean, S, rtol=1e-9, atol=1e-15), "the replicate mean is not the pooled share"
    out = []
    for i in range(512):
        h, m = H[i], mdraw[i]; fh = np.isinf(h)
        s_lo, s_hi = np.maximum(mean[i] - half[i], 0.0), mean[i] + half[i]
        m_lo = np.where(fh, M[i], lo(np.where(fh, 0, h)) / m); m_hi = np.where(fh, M[i], hi(np.where(fh, 0, h)) / m)
        held = mean[i] > 0
        e_lo = np.where(held, s_lo / np.maximum(m_hi, 1e-300), 0.0)
        e_hi = np.where(held, np.where(m_lo > 0, s_hi / np.maximum(m_lo, 1e-300), np.inf), 0.0)
        ok = M[i] > 0; enr = np.where(ok, mean[i] / np.where(ok, M[i], 1.0), -np.inf); w = int(np.argmax(enr))
        riv = [j for j in range(7) if j != w and held[j] and e_hi[j] >= e_lo[w]]
        if riv:
            out.append((i, w, riv))
    return out


def winners_c(S, M, E, name, K=None):
    ok = M > 0
    share = np.where(ok, S, -np.inf)
    enr = np.where(ok, S / np.where(ok, M, 1.0), -np.inf)
    wa, we = share.argmax(1), enr.argmax(1)
    if K is not None:                                             # cells whose most-enriched winner rests on few samples
        srt = np.sort(enr, 1)
        fragile = [i for i in range(512) if K[i, we[i]] < 20 and srt[i, -1] < 2.0 * srt[i, -2]]
        print("   most enriched decided by < 20 samples and < 2x over the runner-up at %d games%s"
              % (len(fragile), (": ipt " + " ".join(map(str, fragile[:40]))) if fragile else ""))
    print("%s: efficiency %.3f .. %.3f, median %.3f, negative at %d games, below 0.9 at %d"
          % (name, E.min(), E.max(), np.median(E), (E < 0).sum(), (E < 0.9).sum()))
    print("   most abundant: " + ", ".join("%s %d" % (ORDER[i], (wa == i).sum()) for i in range(7) if (wa == i).any()))
    print("   most enriched: " + ", ".join("%s %d" % (ORDER[i], (we == i).sum()) for i in range(7) if (we == i).any()))
    unmeasured = ((S > 0) & ~ok).sum()
    if unmeasured:
        print("   WARNING: %d (game, atom) pairs held by the population but of zero measured alpha-measure" % unmeasured)
    return wa, we


RUNS = {"F1": dict(n="$N = 1000$", stems=("Fig4c", "Fig5c")), "F2": dict(n="$N = 100$", stems=("SIFig4c", "SIFig5c"))}


def main(d, run="F1", nd=None):
    """run F1: N = 1000, beta = 100, mu = 1e-2 (Figures 3, 4); run F2: N = 100, beta = 3, mu = 1e-4 (Figures S1, S2).
    nd: where the neutral (kernel-mass) files are, if not with the run -- F2 uses F1's, the kernel being the same."""
    Sc, Mc, Ec, nr, mdraw, UVc, Kc, Hc = load_continuum(d, nd=nd)
    assert (nr == 10).all(), "replicates per game: %s" % np.bincount(nr)
    assert (mdraw > 0).all(), "games without a neutral measure: %d" % (mdraw == 0).sum()
    Sb, Nb, Eb, UVb, Db = wfdata.load("m2", run)
    assert np.allclose(UVb, UVc, atol=1e-9), "the continuum run is not on the 512 games of the binary run"
    wab, web = fig4.winners(Sb, Nb, Eb, Db, "binary M2, %s (alpha = 1/2)" % run)
    wac, wec = winners_c(Sc, Mc, Ec, "continuum M2 %s (alpha = 1/3)" % run, Kc)
    und_p = undetermined(Sc, Mc, Kc, Hc, mdraw)                    # Poisson in the pooled samples: too narrow, for reference
    und = undetermined_rep(load_reps(d), Mc, Hc, mdraw, S=Sc)      # the ten replicates: the test the figure marks
    print("   most enriched NOT determined at 95%% at %d games (Poisson in the pooled samples would say %d)%s"
          % (len(und), len(und_p), (": " + "; ".join(
              "ipt %d %s vs %s" % (i, ORDER[w], "/".join(ORDER[j] for j in r)) for i, w, r in und)) if und else ""))
    print("   neutral draws per game: %s" % ", ".join("%d x %d" % (c, v) for v, c in zip(*np.unique(mdraw, return_counts=True))))

    # the memory-one cube, PartnersRivalsM1f's run, row for row as that project's own fig4c.py computes it (read-only)
    import importlib.util
    M1F = "/Users/martin/Documents/CL1/PartnersRivalsM1f/Figures"
    sys.path.insert(1, M1F)
    spec = importlib.util.spec_from_file_location("m1f_fig4c", os.path.join(M1F, "fig4c.py"))
    m1f = importlib.util.module_from_spec(spec); spec.loader.exec_module(m1f)
    E1, wa1, we1 = m1f.continuous(False, run=run)
    d1 = m1f.cwdata.load(False, uv_ref=wfdata.load("m1", run)[3], run=run)
    assert np.allclose(d1["UV"], UVc, atol=1e-9), "the memory-one run is not on the same 512 games"

    # three rows, as PartnersRivalsM1f's Figure 5 has three: the memory-two cube beneath the two spaces it contains
    nl = RUNS[run]["n"]
    L1, L2, L3 = ("Memory-1 continuous, " + nl, "Memory-2 binary, " + nl, "Memory-2 continuous, " + nl)
    fig4.HERE = fig5.HERE = HERE                                  # write here, never into FinalFigures or PartnersRivalsM1f
    # a dot on each game of the bottom row's most-enriched panel (the ninth axes fig4.draw makes, row by row) whose most
    # enriched atom the ten replicates leave undetermined at 95 %; fig4.draw saves and closes its figure itself, so the
    # dots go in at its savefig
    import matplotlib.figure
    from matplotlib.path import Path
    UD = np.array([UVc[i] / np.sqrt(16.0 + UVc[i] @ UVc[i]) for i, _, _ in und]).reshape(-1, 2)
    pol, keep = wfdata.cells()
    for j, ipt in enumerate(keep):                                 # the site is inside its own cell: the disk map is right
        site = UVc[ipt] / np.sqrt(16.0 + UVc[ipt] @ UVc[ipt])
        assert Path(np.asarray(pol[j])).contains_point(site), "site %d is not in its Voronoi cell" % ipt
    save0 = matplotlib.figure.Figure.savefig

    def save_dotted(fig, *a, **k):
        ax = fig.axes[8]
        P = [a_.get_position() for a_ in fig.axes[:9]]
        assert P[8].x0 == max(q.x0 for q in P) and P[8].y0 == min(q.y0 for q in P), "axes 8 is not the bottom-right panel"
        if not getattr(ax, "_dotted", False):
            ax.scatter(UD[:, 0], UD[:, 1], s=9.0, c="white", edgecolors="0.1", linewidths=0.45, zorder=5)
            ax._dotted = True
        return save0(fig, *a, **k)
    matplotlib.figure.Figure.savefig = save_dotted
    try:
        fig4.draw([(L1, E1, wa1, we1), (L2, Eb, wab, web), (L3, Ec, wac, wec)], RUNS[run]["stems"][0])
    finally:
        matplotlib.figure.Figure.savefig = save0

    K1, K3 = ("m1c", run), ("m2c", "A333" + run)
    wfdata.LABEL[K1], wfdata.LABEL[("m2", run)], wfdata.LABEL[K3] = L1, L2, L3
    base = fig5.shares

    def props(S, W):
        PS = np.stack([S[:, [ORDER.index(a) for a in atoms]].sum(1) for _, _, atoms in fig5.PROPS], 1)
        PN = np.stack([W[:, [ORDER.index(a) for a in atoms]].sum(1) for _, _, atoms in fig5.PROPS], 1)
        return PS, PN

    def shares(space, which):
        if (space, which) == K1:
            return props(d1["S"], d1["K"])
        if (space, which) == K3:
            PS, PN = props(Sc, Mc)
            print("continuum M2 %s (alpha = 1/3): share of the population with the property" % run)
            for j, (nm, code, _) in enumerate(fig5.PROPS):
                v = PS[:, j]
                print("   %-12s %s  mean %.3f  median %.3f  min %.4f  max %.3f  above 1/2 at %3d games  no strategy at %3d"
                      % (nm, code, v.mean(), np.median(v), v.min(), v.max(), (v > 0.5).sum(), (PN[:, j] == 0).sum()))
            return PS, PN
        return base(space, which)
    fig5.shares = shares
    fig5.draw([K1, ("m2", run), K3], RUNS[run]["stems"][1])


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", help="the run's r_1o3_*.txt (and, for F1, n_/nr_1o3_*.txt)")
    ap.add_argument("--run", default="F1", choices=["F1", "F2"])
    ap.add_argument("--neutral", default=None, help="directory of the n_/nr_ kernel-mass files, if not DIR")
    a = ap.parse_args()
    main(a.dir, a.run, a.neutral)
