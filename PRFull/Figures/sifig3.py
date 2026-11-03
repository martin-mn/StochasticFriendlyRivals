#!/usr/bin/env python3
"""
sifig3.py -- SI Figure 3 of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity": convergence of the Wright-Fisher runs.

    python3 sifig3.py       ->  SIFig3.pdf / SIFig3.png (200 dpi), sifig3_stats.txt (everything printed)

The instrumented reruns of run F1 (N = 1000, beta = 100, mu = 1e-2, eps = 1e-4, alpha = 1/3) at four times the production
length, 4e7 generations instead of 1e7, ten independent replicates per game, one replicate per Cannon task:

    memory one   PRFull/cannon/cw2 arm TS    results/TS_series.npz, written by its collect.py from out/TS/p<ipt>_r<rep>/t<isl>:
                 eff_<ipt> (10, 2000): efficiency per block of 2e4 generations; share_<ipt> (10, 2000, 7): the exact share of
                 each atom per block (individual-generations / (N x 2e4)), in the paper's ORDER
    memory two   PRFull/cannon/wfc2 arm TS   out/TS/r_1o3_<ipt+1>_<rep>.txt, rep 1..10, read as analyze.py's parse_ts():
                 TB lines: efficiency per block of 1e5 generations (it0..it1); TS lines: atom counts c0..c7 (wfc's index
                 k = 4E + 2N + R, i.e. the code BITS[k]) of NI = 10 individuals sampled every 2e4 generations

Four games, one per column, the same in every row: W ipt 29 (V_p029_PDW), S ipt 171 (V_p171_SHS), N ipt 500 (V_p500_PDN),
E ipt 415 (V_p415_HAE).  Rows: a, memory one, the share of the leading atom; b, memory one, efficiency; c, d the same for
memory two.  The leading atom is the most abundant atom of the PRODUCTION run at that game, pooled over its ten replicates:
the diagnostics' own choice, read from cw2/results/TS_games.tsv (lead_prod) and wfc2/analysis/TS_games.tsv (prod_lead), and
checked against collect_all.log and against collect.py's rule on production's replicate file.

Drawn: the ten replicates as thin lines, each the mean over consecutive windows of 1e6 generations (50 blocks of 2e4 at memory
one; 50 samples of 10 individuals, or 10 blocks of 1e5, at memory two); the mean over the replicates of the running mean
from generation 0, computed at the native resolution, thick black; the production window (5e6, 1e7] shaded; the start of
the second half, 2e7, dashed.  Printed: the check of the window means (5e6, 1e7] and (2e7, 4e7] against the diagnostics.

Style: the SI Figures 1 and 2 (fig45.py, through PartnersRivals1/FinalFigures/fig4.py and atomshares.py): 14.9 in wide,
DejaVu Sans, font.size 15, axes.linewidth 1.2, FS = 13.05, bold panel letters at FS * 1.7, text colour 0.12, row labels
at FS * 1.25; the constants are copied, not imported (atomshares.py pulls in the disk machinery).  Read-only on every
input; nothing outside PRFull/Figures is written (no bytecode either).  Local, seconds.
"""
import glob
import io
import os
import sys

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.lines import Line2D                                # noqa: E402
from matplotlib.patches import Patch                               # noqa: E402
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MaxNLocator    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PRF = os.path.dirname(HERE)
CW2 = os.path.join(PRF, "cannon", "cw2")
WFC2 = os.path.join(PRF, "cannon", "wfc2")
STEM = "SIFig3"
PNG_DPI = 200

ORDER = ["000", "100", "010", "001", "110", "011", "111"]          # the paper's order (fig45.py, collect.py, analyze.py)
RAW1 = [0, 4, 2, 1, 6, 3, 7]                                       # memory one (cw): raw index 4E+2N+R of each ORDER atom
NPOP = 1000
GAMES = [29, 171, 500, 415]                                        # W, S, N, E
QNAME = {"PD": "Prisoner's Dilemma", "SH": "Stag Hunt", "HA": "Harmony", "SD": "Snowdrift"}
ITEND = 40_000_000
P0, P1, L0 = ITEND // 8, ITEND // 4, ITEND // 2                    # production window (5e6, 1e7]; second half (2e7, 4e7]
SM = 1_000_000                                                     # the drawn replicate lines: means over 1e6 generations

# ---- style of SI Figures 1-2 (PartnersRivals1/FinalFigures/atomshares.py, fig4.py)
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 15, "axes.linewidth": 1.2})
FS = max(9.0, 13.0 * (4.3 / 13.46) * 1.9) * 1.45                   # atomshares.FS = 13.05
TXT = "0.12"
REPCOL, REPALPHA, REPLW = "#5f8fa8", 0.55, 0.9                     # the replicates: a muted blue-grey
CUMLW = 2.4
BAND = "0.88"                                                      # the production window
HALFCOL = "0.30"


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


# --------------------------------------------------------------------------------------------- the games and the leaders
def games():
    """wfc2/games.dat: line = task = ipt + 1; columns gid, ., R, S, T, u, v, tag"""
    G = {}
    for n, l in enumerate([l for l in open(os.path.join(WFC2, "games.dat")) if l.strip() and not l.startswith("#")]):
        f = l.split()
        ipt = int(f[7].split("_")[1][1:])
        assert ipt == n, "games.dat line %d is ipt %d" % (n + 1, ipt)
        G[ipt] = dict(tag=f[7], u=float(f[5]), v=float(f[6]), quad=f[7][-3:-1], wedge=f[7][-1])
    assert len(G) == 512
    return G


def tsv(fn, key):
    rows = [l.rstrip("\n").split("\t") for l in open(fn) if l.strip() and not l.startswith("#")]
    h = rows[0]
    return {int(r[h.index(key)]): dict(zip(h, r)) for r in rows[1:]}


def m1_leaders():
    """cw2/results/TS_games.tsv (lead_prod, E_P, E_L, E_prod); checked against collect.py's rule on the production file"""
    T = tsv(os.path.join(CW2, "results", "TS_games.tsv"), "ipt")
    R = {}
    for l in open(os.path.join(CW2, "ref", "packed", "cw1_e4.rep.tsv")):
        if l.startswith("#") or not l.strip():
            continue
        f = l.split("\t")
        R.setdefault(int(f[0]), []).append([float(x) for x in f[6:13]])
    lead, prod = {}, {}
    for i in GAMES:
        S = np.array(R[i]).mean(0)
        assert len(R[i]) == 10
        assert ORDER[int(S.argmax())] == T[i]["lead_prod"], "memory one ipt %d: lead_prod disagrees with production" % i
        lead[i] = T[i]["lead_prod"]
        prod[i] = dict(E=float(T[i]["E_prod"]), lead=float(S[ORDER.index(lead[i])]))
    return lead, T, prod


def m2_log():
    """the TS rows of wfc2/collect_all.log (the last time-series section): ipt -> the printed columns"""
    L = open(os.path.join(WFC2, "collect_all.log")).read().split("\n")
    k = max(j for j, l in enumerate(L) if "--- time-series diagnostics" in l)
    hdr = L[k + 1]
    assert hdr.split()[:11] == ["ipt", "tag", "split", "lead", "eff", "P", "eff", "L", "eff", "prd", "z"], hdr
    D = {}
    for l in L[k + 2:]:
        f = l.split()
        if not f or not f[0].isdigit():
            if l.strip().startswith("CRITERION"):
                break
            continue
        D[int(f[0])] = dict(tag=f[1], split=f[2], lead=f[3], effP=f[4], effL=f[5], effprd=f[6], leadP=f[8], leadL=f[9],
                            leadpr=f[10])
    return D


def m2_leaders():
    T = tsv(os.path.join(WFC2, "analysis", "TS_games.tsv"), "ipt")
    D = m2_log()
    lead = {}
    for i in GAMES:
        assert T[i]["prod_lead"] == D[i]["lead"], "memory two ipt %d: TS_games.tsv and collect_all.log disagree" % i
        lead[i] = T[i]["prod_lead"]
    return lead, D


# --------------------------------------------------------------------------------------------- the series
def load_m1(i, lead):
    z = np.load(os.path.join(CW2, "results", "TS_series.npz"))
    assert list(z["order"]) == ORDER
    eff, sh = z["eff_%d" % i], z["share_%d" % i]
    assert eff.shape == (10, 2000) and sh.shape == (10, 2000, 7), (eff.shape, sh.shape)
    L = ITEND // eff.shape[1]
    assert L == 20000
    # the npz against the raw task output, where it is still local (collect.py's parsing, independently repeated)
    dirs = sorted(glob.glob(os.path.join(CW2, "out", "TS", "p%03d_r*" % i)))
    if dirs:
        assert len(dirs) == 10
        E, S = [], []
        for d in dirs:
            e, C = [], []
            for l in open(os.path.join(d, "t%d" % (3 * i + 3))):
                if l.startswith("TS"):
                    f = l.split()
                    e.append(float(f[4]))
                    C.append([int(x) for x in f[7:15]])
            C = np.array(C)
            assert len(e) == 2000 and (C.sum(1) == NPOP * L).all() and (C[:, 5] == 0).all()
            E.append(e)
            S.append(C[:, RAW1] / float(NPOP * L))
        assert np.allclose(np.array(E), eff, rtol=0, atol=1e-13) and np.allclose(np.array(S), sh, rtol=0, atol=1e-13), \
            "memory one ipt %d: TS_series.npz differs from out/TS" % i
        raw = "npz = raw out/TS files"
    else:
        raw = "raw out/TS not local"
    t_end = L * np.arange(1, 2001)                                 # block k is generations (t_end - L, t_end]
    return dict(share=sh[:, :, ORDER.index(lead)], eff=eff, ts=t_end, te=t_end, ws=np.ones_like(sh[:, :, 0]),
                Ls=L, Le=L, raw=raw, inst=False)


def load_m2(i, lead):
    """analyze.py parse_ts(), one file per replicate: TB (it0, it1, eff), TS (it, NI, c[8])"""
    k = int(lead, 2)                                               # wfc's index of the atom: BITS[k] == lead
    share, cnt, eff = [], [], []
    for rep in range(1, 11):
        fn = os.path.join(WFC2, "out", "TS", "r_1o3_%d_%d.txt" % (i + 1, rep))
        TB, TS, TSR = [], [], None
        for l in open(fn):
            f = l.split()
            if not f:
                continue
            if f[0] == "TB":
                TB.append((int(f[3]), int(f[4]), float(f[6])))
            elif f[0] == "TS":
                c = [int(x) for x in f[5:13]]
                assert sum(c) == int(f[4])
                TS.append((int(f[3]), int(f[4]), c))
            elif f[0] == "TSR":
                TSR = f
            elif f[0] == "REP":
                assert int(f[2]) == i
        assert TSR is not None and int(TSR[2]) == 2000 and int(TSR[3]) == 10 and int(TSR[4]) == 400, TSR
        assert [x[0] for x in TS] == list(range(20000, ITEND + 1, 20000)), fn
        assert [(x[0], x[1]) for x in TB] == [(1 + 100000 * b, 100000 * (b + 1)) for b in range(400)], fn
        share.append([x[2][k] / x[1] for x in TS])
        cnt.append([x[1] for x in TS])
        eff.append([x[2] for x in TB])
    return dict(share=np.array(share), eff=np.array(eff), ts=np.arange(20000, ITEND + 1, 20000),
                te=np.arange(100000, ITEND + 1, 100000), ws=np.array(cnt, float), Ls=20000, Le=100000, raw="", inst=True)


def binned(x, w, L):
    """means over consecutive windows of SM generations (weights w: individuals per sample; 1 for exact blocks)"""
    m = SM // L
    r, n = x.shape
    assert n % m == 0
    return (x * w).reshape(r, n // m, m).sum(2) / w.reshape(r, n // m, m).sum(2)


def running(x, w):
    """the mean over replicates of the running mean from generation 0"""
    return (np.cumsum(x * w, 1) / np.cumsum(w, 1)).mean(0)


def window(x, w, t_end, L, a, b, inst=False):
    """per replicate, the mean over the blocks or samples inside (a, b], as the diagnostics: a block (t_end - L, t_end]
    inside (a, b] (collect.py's slice, analyze.py's it0 > a and it1 <= b), or a sample at an instant a < t_end <= b"""
    m = ((t_end > a) & (t_end <= b)) if inst else ((t_end - L >= a) & (t_end <= b))
    return (x[:, m] * w[:, m]).sum(1) / w[:, m].sum(1)


# --------------------------------------------------------------------------------------------- checks
def checks(G, D1, D2, lead1, lead2, T1, prod1, log2):
    print("=" * 118)
    print("CHECK: window means recomputed from the drawn series (mean over the 10 replicates of each replicate's window mean)")
    print("  P = the production window (5e6, 1e7]; L = the second half (2e7, 4e7]; 'diag' = what the diagnostics printed")
    print("  memory one: cw2/results/TS_games.tsv (E_P, E_L to 6 significant digits; its lead shares are not printed, so the")
    print("  production pooled share from cw2/ref/packed/cw1_e4.rep.tsv is given instead); memory two: wfc2/collect_all.log")
    print("  (eff to 4 decimals, lead to 3)")
    fmt = "  %-3s %-4s %-3s %-4s | %-8s %-8s %-3s | %-8s %-8s %-3s | %-7s %-7s %-3s | %-7s %-7s %-3s | %-8s %-7s"
    print(fmt % ("mem", "ipt", "wdg", "lead", "eff P", "diag", "ok", "eff L", "diag", "ok", "lead P", "diag", "ok",
                 "lead L", "diag", "ok", "eff prod", "lead pr"))
    allok = True
    rows = []
    for mem, D, lead in (("1", D1, lead1), ("2", D2, lead2)):
        for i in GAMES:
            d = D[i]
            ws = d["ws"]
            we = np.ones_like(d["eff"])
            eP = window(d["eff"], we, d["te"], d["Le"], P0, P1).mean()
            eL = window(d["eff"], we, d["te"], d["Le"], L0, ITEND).mean()
            sP = window(d["share"], ws, d["ts"], d["Ls"], P0, P1, d["inst"]).mean()
            sL = window(d["share"], ws, d["ts"], d["Ls"], L0, ITEND, d["inst"]).mean()
            # the same from the drawn 1e6 means (the window edges are multiples of 1e6, so these must be identical)
            be, bs = binned(d["eff"], we, d["Le"]), binned(d["share"], ws, d["Ls"])
            cut = lambda a, b: slice(a // SM, b // SM)
            assert np.isclose(be[:, cut(P0, P1)].mean(), eP, rtol=0, atol=1e-12)
            assert np.isclose(be[:, cut(L0, ITEND)].mean(), eL, rtol=0, atol=1e-12)
            if mem == "1":                                         # exact blocks: equal weights in every 1e6 window
                assert np.isclose(bs[:, cut(P0, P1)].mean(), sP, rtol=0, atol=1e-12)
                assert np.isclose(bs[:, cut(L0, ITEND)].mean(), sL, rtol=0, atol=1e-12)
            if mem == "1":
                t = T1[i]
                dEP, dEL = t["E_P"], t["E_L"]
                okP = "%.6g" % eP == dEP
                okL = "%.6g" % eL == dEL
                dsP = dsL = "-"
                okSP = okSL = True
                pe, pl = "%.4f" % prod1[i]["E"], "%.3f" % prod1[i]["lead"]
                s_eP, s_eL = "%.6g" % eP, "%.6g" % eL
            else:
                g = log2[i]
                dEP, dEL, dsP, dsL = g["effP"], g["effL"], g["leadP"], g["leadL"]
                okP, okL = "%.4f" % eP == dEP, "%.4f" % eL == dEL
                okSP, okSL = "%.3f" % sP == dsP, "%.3f" % sL == dsL
                pe, pl = g["effprd"], g["leadpr"]
                s_eP, s_eL = "%.6g" % eP, "%.6g" % eL
            ok = okP and okL and okSP and okSL
            allok &= ok
            yn = lambda b: "yes" if b else "NO"
            print(fmt % (mem, i, G[i]["wedge"], lead[i], s_eP, dEP, yn(okP), s_eL, dEL, yn(okL), "%.4f" % sP, dsP,
                         yn(okSP) if dsP != "-" else "-", "%.4f" % sL, dsL, yn(okSL) if dsL != "-" else "-", pe, pl))
            rows.append((mem, i, eP, eL, sP, sL))
    print("  every printed value reproduced to its printed precision: %s" % ("YES" if allok else "NO"))
    for mem, D in (("1", D1), ("2", D2)):
        print("  memory %s: %s" % (mem, "; ".join("ipt %d %s" % (i, D[i]["raw"]) for i in GAMES if D[i]["raw"]) or "-"))
    return allok


# --------------------------------------------------------------------------------------------- the figure
ML, MR, MT, MB = 1.78, 0.14, 0.86, 1.62
COLGAP = 0.80
ROWGAP, MEMGAP = 0.24, 0.62
PH = 2.72
FIG_W = 14.9                                                       # the width of SI Figures 1 and 2 (14.895 in)
PW = (FIG_W - ML - MR - 3 * COLGAP) / 4
FIG_H = MB + 4 * PH + 2 * ROWGAP + MEMGAP + MT


def eff_limits(ys):
    lo = min(float(np.min(y)) for y in ys)
    hi = max(float(np.max(y)) for y in ys)
    pad = 0.12 * (hi - lo)
    return lo - pad, hi + pad


def draw(G, D1, D2, lead1, lead2):
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    xt = np.arange(0, ITEND + 1, 10_000_000)
    rowtop = []
    for r in range(4):
        y = MB + (3 - r) * PH + (3 - r - (1 if r < 2 else 0)) * ROWGAP + (MEMGAP if r < 2 else 0)
        rowtop.append(y + PH)
        mem = 1 if r < 2 else 2
        D, lead = (D1, lead1) if mem == 1 else (D2, lead2)
        kind = "share" if r % 2 == 0 else "eff"
        for c, i in enumerate(GAMES):
            x = ML + c * (PW + COLGAP)
            ax = fig.add_axes([x / FIG_W, y / FIG_H, PW / FIG_W, PH / FIG_H])
            d = D[i]
            if kind == "share":
                v, w, L, t = d["share"], d["ws"], d["Ls"], d["ts"]
            else:
                v, w, L, t = d["eff"], np.ones_like(d["eff"]), d["Le"], d["te"]
            b = binned(v, w, L)
            tb = (np.arange(b.shape[1]) + 0.5) * SM
            cm = running(v, w)
            ax.axvspan(P0, P1, color=BAND, lw=0, zorder=0)
            ax.axvline(L0, color=HALFCOL, lw=1.1, ls=(0, (4, 3)), zorder=1)
            for rr in range(b.shape[0]):
                ax.plot(tb, b[rr], color=REPCOL, alpha=REPALPHA, lw=REPLW, zorder=2, solid_joinstyle="round")
            ax.plot(t, cm, color="black", lw=CUMLW, zorder=4, solid_capstyle="butt")
            ax.set_xlim(0, ITEND)
            ax.xaxis.set_major_locator(FixedLocator(xt))
            if r == 3:
                ax.set_xticklabels(["%d" % (q // 10_000_000) for q in xt])
                ax.set_xlabel(r"generations ($10^7$)", fontsize=FS, color=TXT, labelpad=4)
            else:
                ax.set_xticklabels([])
            if kind == "share":
                ax.set_ylim(0, 1)
                ax.yaxis.set_major_locator(FixedLocator([0, 0.25, 0.5, 0.75, 1.0]))
                if c == 0:
                    ax.set_yticklabels(["0", "", "0.5", "", "1"])
                else:
                    ax.set_yticklabels([])
                lo_hi = float(np.median(b[:, -20:]))               # the atom's name where the lines are not
                ytxt, va = (0.06, "bottom") if lo_hi > 0.5 else (0.94, "top")
                ax.text(0.965, ytxt, r"atom $\mathtt{%s}$" % lead[i], transform=ax.transAxes, ha="right", va=va,
                        fontsize=FS * 1.0, color=TXT, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
            else:
                lo, hi = eff_limits([b, cm[t > SM]])
                ax.set_ylim(lo, hi)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
                step = np.diff(ax.get_yticks())[0]
                nd = max(3, int(np.ceil(-np.log10(step) - 1e-9)))
                ax.yaxis.set_major_formatter(FormatStrFormatter("%%.%df" % nd))
            ax.tick_params(labelsize=FS * 0.82, length=4, width=1.0, colors=TXT, pad=3)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if c == 0:
                ax.set_ylabel("share of the leading atom" if kind == "share" else "efficiency", fontsize=FS, color=TXT,
                              labelpad=8)
            if r == 0:
                g = G[i]
                num = lambda q: ("%.2f" % q).replace("-", "\u2212")
                ax.set_title("$%s$, %s\n$(u, v)$ = (%s, %s)" % (g["wedge"], QNAME[g["quad"]], num(g["u"]), num(g["v"])),
                             fontsize=FS * 1.05, color=TXT, pad=10, linespacing=1.35)
        fig.text(0.50 / FIG_W, (y + PH) / FIG_H, "abcd"[r], fontsize=FS * 1.7, fontweight="bold", ha="left", va="top")
    # memory labels, rotated, one per pair of rows (as SI Figure 1's row labels)
    for k, lab in enumerate(("Memory-1, stochastic, $N = 1000$", "Memory-2, stochastic, $N = 1000$")):
        top, bot = rowtop[2 * k], rowtop[2 * k + 1] - PH
        fig.text(0.24 / FIG_W, 0.5 * (top + bot) / FIG_H, lab, rotation=90, ha="center", va="center",
                 fontsize=FS * 1.25, color=TXT)
    # the key
    h = [Line2D([], [], color=REPCOL, alpha=REPALPHA + 0.15, lw=REPLW * 1.6),
         Line2D([], [], color="black", lw=CUMLW),
         Patch(facecolor=BAND, edgecolor="none"),
         Line2D([], [], color=HALFCOL, lw=1.1, ls=(0, (4, 3)))]
    labs = [r"each replicate, means over $10^6$ generations",
            r"running mean from generation 0, averaged over the replicates",
            r"production window $(5\times10^6,\ 10^7]$",
            r"start of the second half, $2\times10^7$"]
    fig.legend(h, labs, loc="lower center", bbox_to_anchor=(0.5 * (ML + FIG_W - MR) / FIG_W, 0.10 / FIG_H), ncol=2,
               frameon=False, fontsize=FS * 0.95, handlelength=2.6, columnspacing=2.4, labelcolor=TXT)

    # nothing off the sheet
    fig.canvas.draw()
    R = fig.canvas.get_renderer()
    W, H = FIG_W * fig.dpi, FIG_H * fig.dpi
    for t_ in fig.texts + [t for a in fig.axes for t in a.texts + [a.title, a.xaxis.label, a.yaxis.label]
                           + a.get_xticklabels() + a.get_yticklabels()]:
        if t_.get_text() and t_.get_visible():
            bb = t_.get_window_extent(R)
            assert 0 < bb.x0 and bb.x1 < W and 0 < bb.y0 and bb.y1 < H, "text runs off the sheet: %r" % t_.get_text()
    out = os.path.join(HERE, STEM + ".pdf")
    fig.savefig(out)
    fig.savefig(out.replace(".pdf", ".png"), dpi=PNG_DPI)
    plt.close(fig)
    print("wrote PRFull/Figures/%s.pdf and .png (%.2f x %.2f in, PNG at %d dpi)" % (STEM, FIG_W, FIG_H, PNG_DPI))


def main():
    print("sifig3.py -- SI Figure 3 of PRFull: convergence of the F1 runs (N = 1000, beta = 100, mu = 1e-2, eps = 1e-4, "
          "alpha = 1/3), reruns of 4e7 generations, 10 replicates per game")
    G = games()
    lead1, T1, prod1 = m1_leaders()
    lead2, log2 = m2_leaders()
    for i in GAMES:
        g = G[i]
        print("  ipt %3d %-11s wedge %s  (u, v) = (%.4f, %.4f)  production leader: memory one %s, memory two %s"
              % (i, g["tag"], g["wedge"], g["u"], g["v"], lead1[i], lead2[i]))
        assert log2[i]["tag"] == g["tag"][-3:]
    D1 = {i: load_m1(i, lead1[i]) for i in GAMES}
    D2 = {i: load_m2(i, lead2[i]) for i in GAMES}
    ok = checks(G, D1, D2, lead1, lead2, T1, prod1, log2)
    print("=" * 118)
    draw(G, D1, D2, lead1, lead2)
    sys.stdout.flush()
    with open(os.path.join(HERE, "sifig3_stats.txt"), "w") as fh:
        fh.write(LOG.getvalue())
    sys.__stdout__.write("wrote %s\n" % os.path.join(HERE, "sifig3_stats.txt"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
