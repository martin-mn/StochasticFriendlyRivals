#!/usr/bin/env python3
"""
sitable4.py -- SI Table 4 of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity": the four strategy spaces in the two
populations, every number computed from the data.

    python3 sitable4.py   ->  ../Text/si_table4.tex   the table (\\input by ms.tex at the end of SI section "The evolutionary maps")
                          ->  ../Text/si_table4_body.tex   the same without its legend (deposited in release 1.0.0)
                              sitable4.out            every number with its provenance, the cross-checks against the fact
                                                      sheets (../notes/facts_m1.md 12.6, facts_m2.md 15, facts_ref.md B9)
                                                      and the loaders' own log

Columns: the 16 binary memory-one strategies, the memory-one cube [0,1]^4, the 65536 binary memory-two strategies, the
memory-two cube [0,1]^16.  Blocks: run F1 (N = 1000, beta = 100, mu = 1e-2) and run F2 (N = 100, beta = 3, mu = 1e-4),
both at eps = 1e-4, 10 replicates per game, the second half of each sampled, at the 512 sampled games.

The data are read exactly as the figures read them, through fig45.py (imported, not run: its main() is not called, so
nothing is drawn):
    binary memory one / two   fig45.load_binary(space, run) = PartnersRivals1/FinalFigures/wfdata.load + fig4.winners
    memory-one cube           fig45.load_m1(run)  = PartnersRivalsM1f Figures/fig4c.continuous() + cwdata.load
                              (packs cannon/cw/packed/cw{1,2}_e4.tsv, replicates cw{1,2}_e4.rep.tsv)
    memory-two cube           fig45.load_m2(run)  = PartnersRivalsM2f Figures/figwfc.load_continuum / winners_c /
                              load_reps (check/wfc_a333 for F1, check/wfc_f2 for F2)
The per-replicate atom shares of the binary memory-two runs are not in wfdata (pooled only); they are the seeded re-runs
f1_e4 / f2_e4 of PartnersRivals1/Data, read by FinalFigures/wfrep.load, whose pooled shares and efficiencies are asserted
here to reproduce wfdata's (as FinalFigures/sifig11.py asserts).  The binary memory-one runs (DiskM1WF mw1/mw2) keep per
replicate only pay, ef and the four marginal cooperation probabilities (DiskM1WF/Data/mw?_e4.rep.tsv), from which the
atom shares cannot be recovered; that cell is a dash.

The leading (most abundant) atom is each loader's own: fig4.winners (binary; atoms with no strategy excluded),
fig4c.continuous (memory-one cube; atoms with kappa = 0 and share 0 excluded), figwfc.winners_c (memory-two cube; atoms
of zero measured kernel mass excluded) -- the leaders of Figures 4 and 5 and SI Figures 1 and 2.  The regions are the
tags of PartnersRivalsM2f/cannon/wfc/games.dat, which fig45 checks against the signs of u, v, u+v-1, u-v-1.

Read-only outside PRFull; no bytecode is written; matplotlib's config directory is a temporary one, removed at exit.
"""
import atexit
import contextlib
import io
import os
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
_MPL = tempfile.mkdtemp(prefix="sitable4_mpl_")
os.environ["MPLCONFIGDIR"] = _MPL
atexit.register(shutil.rmtree, _MPL, True)

import numpy as np                                                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.normpath(os.path.join(HERE, os.pardir, "Text", "si_table4.tex"))
OUT = os.path.join(HERE, "sitable4.out")
NOTES = os.path.normpath(os.path.join(HERE, os.pardir, "notes"))

LOADLOG = io.StringIO()
_stdout = sys.stdout
sys.path.insert(0, HERE)
with contextlib.redirect_stdout(LOADLOG):
    import fig45                                                   # noqa: E402  (sets sys.stdout to its Tee: undone below)
sys.stdout = _stdout
import wfrep                                                       # noqa: E402  FinalFigures (on sys.path via fig45)

ORDER = fig45.ORDER                                                # 000 100 010 001 110 011 111
IX = fig45.IX
QUAD, WEDGE = fig45.QUAD, fig45.WEDGE
PD = QUAD == "PD"
ABOVE = np.isin(WEDGE, ["N", "E"])
INS = WEDGE == "S"
assert (PD.sum(), ABOVE.sum(), INS.sum(), (~INS).sum()) == (129, 207, 126, 386), "the regions of the 512 games moved"
PROPS = fig45.PROPS0                                               # efficient 1**, stable *1*, competitive **1
assert [p[0] for p in PROPS] == ["efficient", "stable", "competitive"]
assert PROPS[2][2] == ["001", "011", "111"] and PROPS[0][2] == ["100", "110", "111"]
PARTNER = [IX["110"], IX["111"]]
EFFIC = [IX[a] for a in PROPS[0][2]]
I111 = IX["111"]

SPACES = [("b1", "binary memory one"), ("c1", "memory-one cube"), ("b2", "binary memory two"), ("c2", "memory-two cube")]
RUNS = [("F1", "$N=1000$, $\\beta=100$, $\\mu=10^{-2}$", "N = 1000, beta = 100, mu = 1e-2"),
        ("F2", "$N=100$, $\\beta=3$, $\\mu=10^{-4}$", "N = 100, beta = 3, mu = 1e-4")]
REPRUN = {"F1": "f1_e4", "F2": "f2_e4"}                            # the seeded binary memory-two re-runs with .rep files

SRC = {
    "b1": "wfdata.load('m1', run) [DiskM1WF/Data/mw{1,2}_e4.tsv via DiskM1WF/Opt/m1load] + fig4.winners",
    "c1": "PartnersRivalsM1f fig4c.continuous() + cwdata.load [cannon/cw/packed/cw{1,2}_e4.tsv, kappa_e4.tsv]",
    "b2": "wfdata.load('m2', run) [DiskM2WF/Opt figF4.collect, wfload.load_run] + fig4.winners",
    "c2": "PartnersRivalsM2f figwfc.load_continuum/winners_c [check/wfc_a333 (F1), check/wfc_f2 (F2); masses wfc_a333]",
}
RSRC = {
    "b1": None,
    "c1": "cannon/cw/packed/cw{1,2}_e4.rep.tsv (fig45.load_m1: 512 x 10 replicate shares, mean = pooled share)",
    "b2": "PartnersRivals1/Data/{f1,f2}_e4.rep via FinalFigures/wfrep.load (seeded re-run; pooled S, E = wfdata's)",
    "c2": "check/wfc_a333 or wfc_f2 r_1o3_*.txt REP lines via figwfc.load_reps (111 in W exact)",
}


# --------------------------------------------------------------------------------------------- the data
def load(run):
    """per space: E (512,), wa (512,), S (512, 7), PS (512, 3), R (512, 10, 7) or None, and notes"""
    D = {}
    with contextlib.redirect_stdout(LOADLOG):
        print("=" * 110)
        print("LOADING run %s" % run)
        b1 = fig45.load_binary("m1", run)
        c1 = fig45.load_m1(run)
        b2 = fig45.load_binary("m2", run)
        c2 = fig45.load_m2(run)
        g = wfrep.load(REPRUN[run])
    # the memory-one cube: the replicate mean is the pooled share (fig45 prints it; assert it here)
    assert np.allclose(c1["R"].mean(1), c1["S"], rtol=1e-6, atol=1e-9), "memory-one cube: replicate mean != pooled"
    assert np.allclose(c2["R"].mean(1), c2["S"], rtol=1e-9, atol=1e-15), "memory-two cube: replicate mean != pooled"
    # the binary memory-two re-run reproduces Figure 4's run (sifig11.main_run's check)
    S0, N0, E0, UV0, _ = fig45.wfdata.load("m2", run)
    same = np.allclose(g["S"], S0, atol=1e-6) and np.allclose(g["E"], E0, atol=1e-6)
    assert same, "the seeded re-run %s is not the run of wfdata.load('m2', %s)" % (REPRUN[run], run)
    assert (g["nrep"] == 10).all() and not np.isnan(g["R"]).any()
    assert np.allclose(g["R"].mean(1), g["S"], atol=2e-6)
    assert (np.where(N0 > 0, g["S"], -np.inf).argmax(1) == b2["wa"]).all(), "re-run's pooled leader != wfdata's"
    note_b2 = ("seeded re-run %s: pooled shares within %.1e and efficiencies within %.1e of wfdata.load('m2', '%s')"
               % (REPRUN[run], np.abs(g["S"] - S0).max(), np.abs(g["E"] - E0).max(), run))
    for key, d, R in (("b1", b1, None), ("c1", c1, c1["R"]), ("b2", b2, g["R"]), ("c2", c2, c2["R"])):
        raw = d["S"].argmax(1)
        D[key] = dict(E=d["E"], wa=d["wa"], S=d["S"], PS=d["PS"], R=R, rawdiff=int((raw != d["wa"]).sum()))
    D["b2"]["note"] = note_b2
    return D


# --------------------------------------------------------------------------------------------- the rows
def rows_of(d):
    E, wa, S, PS, R = d["E"], d["wa"], d["S"], d["PS"], d["R"]
    v = {}
    v["med_eff"] = float(np.median(E))
    v["n_lt09"] = int((E < 0.9).sum())
    v["n_neg"] = int((E < 0).sum())
    v["med_eff_pd"] = float(np.median(E[PD]))
    v["pd_partner"] = int(np.isin(wa[PD], PARTNER).sum())
    v["above_eff"] = int(np.isin(wa[ABOVE], EFFIC).sum())
    v["s_111"] = int((wa[INS] == I111).sum())
    v["s_111_share"] = float(np.median(S[INS, I111]))
    v["out_111"] = int((wa[~INS] == I111).sum())
    v["med_comp"] = float(np.median(PS[:, 2]))
    v["mean_eff"] = float(PS[:, 0].mean())
    v["mean_stab"] = float(PS[:, 1].mean())
    v["mean_comp"] = float(PS[:, 2].mean())
    if R is None:
        v["agree10"] = None
    else:
        assert R.shape == (512, 10, 7)
        v["agree10"] = int((R.argmax(2) == wa[:, None]).all(1).sum())
    return v


ROWS = [  # key, label in the table, kind, provenance of the computation
    ("med_eff", "median efficiency", "eff", "np.median(E) over the 512 games; E = mean per-round payoff / Emax"),
    ("n_lt09", "efficiency below $0.9$", "n", "#games with E < 0.9"),
    ("n_neg", "negative mean payoff", "n", "#games with E < 0 (Emax > 0, so E < 0 iff the mean payoff is negative)"),
    ("med_eff_pd", "median efficiency in the PD", "eff", "np.median(E) over the 129 PD games (games.dat tag PD: u<0<v)"),
    ("pd_partner", "PD led by a partner atom", "n", "#PD games whose most abundant atom is 110 or 111"),
    ("above_eff", "above the switch line, led by an efficient atom", "n",
     "#games of wedges N, E (u+v>1; 207) whose most abundant atom is 100, 110 or 111"),
    ("s_111", "$S$, led by \\texttt{111}", "n", "#games of wedge S (126) whose most abundant atom is 111"),
    ("s_111_share", "$S$, median share of \\texttt{111}", "share", "np.median(S[:, 111]) over the 126 games of S"),
    ("out_111", "outside $S$, led by \\texttt{111}", "n", "#games outside S (386) whose most abundant atom is 111"),
    ("med_comp", "median share, competitive", "share",
     "np.median over 512 games of the share of 001 + 011 + 111 (fig45 PS[:, 2])"),
    ("mean_eff", "mean share, efficient", "share", "mean over 512 games of the share of 100 + 110 + 111"),
    ("mean_stab", "mean share, stable", "share",
     "mean over 512 games of the share of 010 + 110 + 011 + 111"),
    ("mean_comp", "mean share, competitive", "share", "mean over 512 games of the share of 001 + 011 + 111"),
    ("agree10", "all ten replicates name the leading atom", "n",
     "#games at which every replicate's own most abundant atom (argmax of its 7 atom shares) is the pooled leader"),
]


def tex_num(x, kind):
    if x is None:
        return "--"
    if kind == "n":
        return "$%d$" % x
    if kind == "eff":
        return "$%.3f$" % x
    if x == 0:
        return "$0$"
    if x >= 0.01:
        return "$%.3f$" % x
    m, e = ("%.1e" % x).split("e")
    return "$%s\\times10^{%d}$" % (m, int(e)) if x < 1e-3 else "$%s$" % ("%.2g" % x)


def txt_num(x, kind):
    if x is None:
        return "-"
    return "%d" % x if kind == "n" else "%.6g" % x


# --------------------------------------------------------------------------------------------- the table
# The legend of SI Table 4 is part of the paper: it lives in sitable4_legend.tex, which the code-and-data release
# 1.0.0 leaves out (it is added with the paper in release 1.1.0).  Without it the script writes the table without its
# legend; si_table4_body.tex, written in every case, is the table without the legend line.
LEGEND_FILE = os.path.join(HERE, "sitable4_legend.tex")
LEGEND = open(LEGEND_FILE).read().rstrip("\n") if os.path.exists(LEGEND_FILE) else None
TEX_BODY = os.path.normpath(os.path.join(HERE, os.pardir, "Text", "si_table4_body.tex"))


def table(V, legend=True):
    L = ["% SI Table 4 -- written by Figures/sitable4.py from the data; do not edit by hand (numbers and provenance:",
         "% Figures/sitable4.out)."]
    if legend and LEGEND is not None:
        L.append("\\sitable{SI Table 4: The four strategy spaces in the two populations.}{%s}" % LEGEND)
    elif legend:
        L.append("% (the legend, sitable4_legend.tex, is part of the paper and not in this release)")
    L += ["\\begin{center}\\small",
         "\\begin{tabular}{@{}lrrrr@{}}",
         "\\toprule",
         " & binary & memory-one & binary & memory-two\\\\",
         " & memory one & cube & memory two & cube\\\\"]
    for k, (run, lab, _) in enumerate(RUNS):
        L.append("\\midrule")
        L.append("\\multicolumn{5}{@{}l}{%s}\\\\" % lab)
        L.append("\\addlinespace[2pt]")
        for key, rl, kind, _ in ROWS:
            L.append("%s & %s\\\\" % (rl, " & ".join(tex_num(V[run][s][key], kind) for s, _ in SPACES)))
    L += ["\\bottomrule", "\\end{tabular}", "\\end{center}", ""]
    return "\n".join(L)


# --------------------------------------------------------------------------------------------- the cross-checks
# (run, space, key, stated value as written, source file:line).  A stated value is compared after rounding the computed
# one to the number of decimals it is written with; counts exactly.
F_M1 = "notes/facts_m1.md"
F_M2 = "notes/facts_m2.md"
F_RF = "notes/facts_ref.md"
CHECKS = [
    # facts_m1.md 12.6
    ("F1", "c1", "med_eff", "0.9742", F_M1 + ":1007"), ("F1", "b1", "med_eff", "0.8016", F_M1 + ":1007"),
    ("F1", "b2", "med_eff", "0.9834", F_M1 + ":1007"),
    ("F1", "c1", "n_lt09", "89", F_M1 + ":1008"), ("F1", "b1", "n_lt09", "324", F_M1 + ":1008"),
    ("F1", "b2", "n_lt09", "10", F_M1 + ":1008"),
    ("F1", "c1", "n_neg", "0", F_M1 + ":1009"), ("F1", "b1", "n_neg", "77", F_M1 + ":1009"),
    ("F1", "b1", "pd_partner", "20", F_M1 + ":1010"), ("F1", "c1", "pd_partner", "128", F_M1 + ":1010"),
    ("F1", "c1", "above_eff", "162", F_M1 + ":1011"),
    ("F1", "c1", "s_111", "126", F_M1 + ":1017"), ("F1", "c1", "s_111_share", "0.9612", F_M1 + ":1017"),
    ("F1", "b1", "s_111", "126", F_M1 + ":1017"), ("F1", "b1", "s_111_share", "0.9653", F_M1 + ":1017"),
    ("F1", "c1", "out_111", "0", F_M1 + ":1018"), ("F1", "b2", "out_111", "60", F_M1 + ":1018"),
    ("F1", "c1", "agree10", "487", F_M1 + ":1019"),
    ("F2", "c1", "n_lt09", "185", F_M1 + ":1025"),
    ("F2", "c1", "med_eff_pd", "0.6516", F_M1 + ":1027"), ("F2", "b1", "med_eff_pd", "0.2495", F_M1 + ":1027"),
    ("F2", "c1", "s_111", "126", F_M1 + ":1030"), ("F2", "c1", "s_111_share", "0.7759", F_M1 + ":1030"),
    ("F2", "c1", "agree10", "434", F_M1 + ":1031"),
    ("F1", "c1", "mean_eff", "0.865", F_M1 + ":1051"), ("F1", "c1", "mean_stab", "0.738", F_M1 + ":1051"),
    ("F1", "c1", "mean_comp", "0.267", F_M1 + ":1051"),
    ("F2", "c1", "mean_eff", "0.775", F_M1 + ":1051"), ("F2", "c1", "mean_stab", "0.842", F_M1 + ":1051"),
    ("F2", "c1", "mean_comp", "0.309", F_M1 + ":1051"),
    ("F2", "c1", "pd_partner", "93", F_M1 + ":1052-1053 (PD leaders {110: 93, 011: 33, 010: 3})"),
    ("F1", "b1", "med_eff", "0.802", F_M1 + ":1059"), ("F1", "c1", "med_eff", "0.974", F_M1 + ":1059"),
    ("F1", "b2", "med_eff", "0.983", F_M1 + ":1059"),
    ("F1", "b1", "n_lt09", "324", F_M1 + ":1059"), ("F1", "c1", "n_lt09", "89", F_M1 + ":1059"),
    ("F1", "b2", "n_lt09", "10", F_M1 + ":1059"),
    ("F1", "b1", "n_neg", "77", F_M1 + ":1059"), ("F1", "c1", "n_neg", "0", F_M1 + ":1059"),
    ("F1", "b2", "n_neg", "0", F_M1 + ":1059"),
    ("F1", "b1", "mean_eff", "0.35", F_M1 + ":1059-1060"), ("F1", "b1", "mean_stab", "0.54", F_M1 + ":1059-1060"),
    ("F1", "b1", "mean_comp", "0.59", F_M1 + ":1059-1060"),
    ("F1", "c1", "mean_eff", "0.87", F_M1 + ":1059-1060"), ("F1", "c1", "mean_stab", "0.74", F_M1 + ":1059-1060"),
    ("F1", "c1", "mean_comp", "0.27", F_M1 + ":1059-1060"),
    ("F1", "b2", "mean_eff", "0.95", F_M1 + ":1059-1060"), ("F1", "b2", "mean_stab", "0.94", F_M1 + ":1059-1060"),
    ("F1", "b2", "mean_comp", "0.37", F_M1 + ":1059-1060"),
    ("F2", "b1", "med_eff", "0.753", F_M1 + ":1061"), ("F2", "c1", "med_eff", "0.968", F_M1 + ":1061"),
    ("F2", "b2", "med_eff", "0.987", F_M1 + ":1061"),
    ("F2", "b1", "n_lt09", "331", F_M1 + ":1061"), ("F2", "c1", "n_lt09", "185", F_M1 + ":1061"),
    ("F2", "b2", "n_lt09", "91", F_M1 + ":1061"),
    ("F2", "b1", "n_neg", "57", F_M1 + ":1061"), ("F2", "c1", "n_neg", "0", F_M1 + ":1061"),
    ("F2", "b2", "n_neg", "1", F_M1 + ":1061"),
    ("F2", "b1", "mean_eff", "0.37", F_M1 + ":1062"), ("F2", "b1", "mean_stab", "0.62", F_M1 + ":1062"),
    ("F2", "b1", "mean_comp", "0.61", F_M1 + ":1062"),
    ("F2", "c1", "mean_eff", "0.78", F_M1 + ":1062"), ("F2", "c1", "mean_stab", "0.84", F_M1 + ":1062"),
    ("F2", "c1", "mean_comp", "0.31", F_M1 + ":1062"),
    ("F2", "b2", "mean_eff", "0.85", F_M1 + ":1062"), ("F2", "b2", "mean_stab", "0.86", F_M1 + ":1062"),
    ("F2", "b2", "mean_comp", "0.27", F_M1 + ":1062"),
    # facts_m2.md 15 (order there: m2 cube / m2 binary / m1 cube)
    ("F1", "c2", "med_eff", "0.987", F_M2 + ":861"), ("F1", "b2", "med_eff", "0.983", F_M2 + ":861"),
    ("F1", "c1", "med_eff", "0.974", F_M2 + ":861"),
    ("F1", "c2", "n_lt09", "7", F_M2 + ":861"), ("F1", "b2", "n_lt09", "10", F_M2 + ":861"),
    ("F1", "c1", "n_lt09", "89", F_M2 + ":861"),
    ("F1", "c2", "n_neg", "0", F_M2 + ":861-862"), ("F1", "b2", "n_neg", "0", F_M2 + ":861-862"),
    ("F1", "c1", "n_neg", "0", F_M2 + ":861-862"),
    ("F1", "c2", "pd_partner", "128", F_M2 + ":864"), ("F1", "b2", "pd_partner", "128", F_M2 + ":864"),
    ("F1", "c1", "pd_partner", "128", F_M2 + ":864"),
    ("F1", "c2", "med_eff_pd", "0.984", F_M2 + ":865"), ("F1", "b2", "med_eff_pd", "0.984", F_M2 + ":865"),
    ("F1", "c1", "med_eff_pd", "0.974", F_M2 + ":865"),
    ("F1", "c2", "above_eff", "207", F_M2 + ":866"), ("F1", "b2", "above_eff", "207", F_M2 + ":866"),
    ("F1", "c1", "above_eff", "162", F_M2 + ":866"),
    ("F1", "c2", "s_111", "126", F_M2 + ":868"), ("F1", "b2", "s_111", "126", F_M2 + ":868"),
    ("F1", "c1", "s_111", "126", F_M2 + ":868"),
    ("F1", "c2", "s_111_share", "0.867", F_M2 + ":868"), ("F1", "b2", "s_111_share", "0.919", F_M2 + ":868"),
    ("F1", "c1", "s_111_share", "0.961", F_M2 + ":868"),
    ("F1", "c2", "out_111", "4", F_M2 + ":869"), ("F1", "b2", "out_111", "60", F_M2 + ":869"),
    ("F1", "c1", "out_111", "0", F_M2 + ":869"),
    ("F1", "c2", "med_comp", "0.0016", F_M2 + ":873"), ("F1", "b2", "med_comp", "0.1101", F_M2 + ":873"),
    ("F1", "c1", "med_comp", "0.0266", F_M2 + ":873"),
    ("F1", "c2", "agree10", "477", F_M2 + ":878"), ("F1", "b2", "agree10", "487", F_M2 + ":878"),
    ("F2", "c2", "n_lt09", "124", F_M2 + ":885"), ("F2", "b2", "n_lt09", "91", F_M2 + ":885"),
    ("F2", "c1", "n_lt09", "185", F_M2 + ":885"),
    ("F2", "c2", "med_eff_pd", "0.816", F_M2 + ":885-886"), ("F2", "b2", "med_eff_pd", "0.832", F_M2 + ":885-886"),
    ("F2", "c1", "med_eff_pd", "0.652", F_M2 + ":885-886"),
    ("F2", "c2", "n_neg", "0", F_M2 + ":886-887"), ("F2", "b2", "n_neg", "1", F_M2 + ":886-887"),
    ("F2", "c2", "med_eff", "0.981", F_M2 + ":887"), ("F2", "b2", "med_eff", "0.987", F_M2 + ":887"),
    ("F2", "c1", "med_eff", "0.968", F_M2 + ":887"),
    ("F2", "c2", "pd_partner", "95", F_M2 + ":888"), ("F2", "b2", "pd_partner", "87", F_M2 + ":888"),
    ("F2", "c1", "pd_partner", "93", F_M2 + ":888"),
    ("F2", "c2", "above_eff", "181", F_M2 + ":890"), ("F2", "b2", "above_eff", "206", F_M2 + ":890"),
    ("F2", "c1", "above_eff", "144", F_M2 + ":890"),
    ("F2", "c2", "s_111", "23", F_M2 + ":892"), ("F2", "b2", "s_111", "42", F_M2 + ":892"),
    ("F2", "c1", "s_111", "126", F_M2 + ":892"),
    ("F2", "c2", "s_111_share", "0.211", F_M2 + ":893"),
    ("F2", "c2", "out_111", "0", F_M2 + ":896"), ("F2", "b2", "out_111", "47", F_M2 + ":896"),
    ("F2", "c2", "med_comp", "0.041", F_M2 + ":897 (4.1%)"), ("F2", "b2", "med_comp", "0.20", F_M2 + ":897 (20%)"),
    ("F2", "c1", "med_comp", "0.13", F_M2 + ":897 (13%)"),
    ("F2", "c2", "agree10", "424", F_M2 + ":898"), ("F2", "b2", "agree10", "412", F_M2 + ":898"),
    # facts_ref.md B9 (order there: M1 N=1000 / M1 N=100 / M2 N=1000 / M2 N=100)
    ("F1", "b1", "med_eff", "0.802", F_RF + ":685"), ("F2", "b1", "med_eff", "0.753", F_RF + ":685"),
    ("F1", "b2", "med_eff", "0.983", F_RF + ":685"), ("F2", "b2", "med_eff", "0.987", F_RF + ":685"),
    ("F1", "b1", "n_lt09", "324", F_RF + ":687"), ("F2", "b1", "n_lt09", "331", F_RF + ":687"),
    ("F1", "b2", "n_lt09", "10", F_RF + ":687"), ("F2", "b2", "n_lt09", "91", F_RF + ":687"),
    ("F1", "b1", "n_neg", "77", F_RF + ":688"), ("F2", "b1", "n_neg", "57", F_RF + ":688"),
    ("F1", "b2", "n_neg", "0", F_RF + ":688"), ("F2", "b2", "n_neg", "1", F_RF + ":688"),
    ("F1", "b1", "mean_eff", "0.35", F_RF + ":693"), ("F1", "b1", "mean_stab", "0.54", F_RF + ":693"),
    ("F1", "b1", "mean_comp", "0.59", F_RF + ":693"),
    ("F2", "b1", "mean_eff", "0.37", F_RF + ":693"), ("F2", "b1", "mean_stab", "0.62", F_RF + ":693"),
    ("F2", "b1", "mean_comp", "0.61", F_RF + ":693"),
    ("F1", "b2", "mean_eff", "0.95", F_RF + ":693"), ("F1", "b2", "mean_stab", "0.94", F_RF + ":693"),
    ("F1", "b2", "mean_comp", "0.37", F_RF + ":693"),
    ("F2", "b2", "mean_eff", "0.85", F_RF + ":693"), ("F2", "b2", "mean_stab", "0.86", F_RF + ":693"),
    ("F2", "b2", "mean_comp", "0.27", F_RF + ":693"),
    ("F1", "b1", "s_111", "126", F_RF + ":702"), ("F1", "b1", "s_111_share", "0.97", F_RF + ":702"),
    ("F2", "b1", "s_111_share", "1.00", F_RF + ":702-703"),
    ("F1", "b2", "s_111", "126", F_RF + ":713"),
    ("F1", "b2", "out_111", "60", F_RF + ":713-714 (27 Harmony + 14 PD above the line + 19 deep PD)"),
    ("F1", "b2", "med_eff_pd", "0.984", F_RF + ":721"), ("F2", "b2", "med_eff_pd", "0.832", F_RF + ":721"),
    ("F1", "b1", "med_eff_pd", "0.623", F_RF + ":721"), ("F2", "b1", "med_eff_pd", "0.249", F_RF + ":721"),
    ("F1", "b1", "pd_partner", "20", F_RF + ":722"), ("F2", "b1", "pd_partner", "15", F_RF + ":722"),
    ("F1", "b2", "pd_partner", "128", F_RF + ":722-723"), ("F2", "b2", "pd_partner", "87", F_RF + ":722-723"),
    ("F2", "b2", "s_111_share", "0.357", F_RF + ":723"),
    ("F2", "b2", "agree10", "412", F_RF + ":726"), ("F1", "b2", "agree10", "487", F_RF + ":726"),
]


def check(V):
    out, bad = [], []
    for run, sp, key, stated, src in CHECKS:
        x = V[run][sp][key]
        kind = [r[2] for r in ROWS if r[0] == key][0]
        if "." in stated:
            nd = len(stated.split(".")[1])
            ok = x is not None and round(x, nd) == round(float(stated), nd)
            shown = "%.*f" % (nd, x) if x is not None else "-"
            full = "%.6g" % x if x is not None else "-"
        else:
            ok = x is not None and x == int(stated)
            shown = full = txt_num(x, kind)
        line = "  %-4s %-3s %-12s stated %-8s computed %-8s (%-10s) %-5s %s" % (
            run, sp, key, stated, shown, full, "ok" if ok else "DIFF", src)
        out.append(line)
        if not ok:
            bad.append(line)
    return out, bad


STATS = os.path.join(HERE, "fig45_stats.txt")
STATS_HEAD = {("F1", "c1"): "Memory-1, stochastic [0,1]^4, run F1", ("F1", "c2"): "Memory-2, stochastic [0,1]^16, run F1",
              ("F2", "c1"): "Memory-1, stochastic [0,1]^4, run F2", ("F2", "c2"): "Memory-2, stochastic [0,1]^16, run F2",
              ("F1", "b1"): "Memory-1, binary (16), run F1", ("F1", "b2"): "Memory-2, binary (65536), run F1",
              ("F2", "b1"): "Memory-1, binary (16), run F2", ("F2", "b2"): "Memory-2, binary (65536), run F2"}


def stats_check(V):
    """the same rows read back from fig45_stats.txt (fig45.census, the text behind Figures 4, 5 and SI Figures 1, 2),
    compared at the precision printed there"""
    import re
    lines = open(STATS).read().split("\n")
    out, bad = [], []
    for (run, sp), head in STATS_HEAD.items():
        i0 = lines.index(head)
        blk = []
        for l in lines[i0 + 1:]:
            if l.startswith("-" * 20) or l.startswith("=" * 20) or l.startswith("#" * 20):
                break
            blk.append(l)
        sec = {}
        cur = None
        for l in blk:
            if l.startswith("  ") and not l.startswith("    "):
                cur = l.strip()
            elif cur is not None:
                sec.setdefault(cur.split(":")[0], []).append(l.split())
        eff = {" ".join(r[:-6]): r[-6:] for r in sec["efficiency (mean per-round payoff / Emax)"][1:]}
        ab = {" ".join(r[:-8]): r[-8:] for r in sec["most abundant atom"][1:]}
        ps = {" ".join(r[:-13]): r[-13:] for r in sec["share of the population with each property"][1:]}
        txt = "\n".join(blk)
        got = {"med_eff": eff["all"][2], "n_lt09": eff["all"][4], "n_neg": eff["all"][5],
               "med_eff_pd": eff["PD (u<0<v)"][2],
               "s_111": ab["wedge S"][7], "out_111": str(int(ab["all"][7]) - int(ab["wedge S"][7])),
               "med_comp": ps["all"][10], "mean_eff": ps["all"][1], "mean_stab": ps["all"][5], "mean_comp": ps["all"][9],
               "pd_partner": re.search(r"partners lead \(most abundant atom 110 or 111\) at: .*?PD \(u<0<v\) (\d+)",
                                       txt).group(1),
               "above_eff": re.search(r"an efficient atom \(100, 110, 111\) leads at: .*?above switch (\d+)", txt).group(1),
               "s_111_share": re.search(r"friendly rivals \(111\): median share per wedge: .*?S ([0-9.e+-]+),",
                                        txt).group(1)}
        m = re.search(r"all ten replicates name the pooled leader at (\d+) of 512", txt)
        if m:
            got["agree10"] = m.group(1)
        for key, stated in got.items():
            x = V[run][sp][key]
            if "." in stated or "e" in stated:
                sig = len(stated.replace(".", "").lstrip("0").split("e")[0]) if "e" in stated else None
                if sig is None:
                    nd = len(stated.split(".")[1])
                    ok = round(x, nd) == round(float(stated), nd)
                else:
                    ok = float("%.*g" % (sig, x)) == float(stated)
            else:
                ok = x is not None and x == int(stated)
            line = "  %-4s %-3s %-12s fig45_stats %-10s computed %-12s %s" % (run, sp, key, stated,
                                                                            "-" if x is None else "%.6g" % x,
                                                                            "ok" if ok else "DIFF")
            out.append(line)
            if not ok:
                bad.append(line)
    return out, bad


# --------------------------------------------------------------------------------------------- details for the .out
def details(run, D):
    L = []
    for sp, name in SPACES:
        d = D[sp]
        wa, S, R = d["wa"], d["S"], d["R"]
        L.append("  %s (%s)" % (name, sp))
        L.append("    data: %s" % SRC[sp])
        L.append("    replicates: %s" % (RSRC[sp] or "NONE -- wfdata.load('m1', run) is pooled only, and DiskM1WF/Data/"
                                                 "mw%s_e4.rep.tsv carries per replicate only pay, ef, ps1..ps4" % run[1]))
        if "note" in d:
            L.append("    " + d["note"])
        L.append("    most abundant atom (loader's rule), all 512: " + ", ".join(
            "%s %d" % (ORDER[i], (wa == i).sum()) for i in range(7) if (wa == i).any())
            + "; raw S.argmax differs at %d games" % d["rawdiff"])
        for rn, m in (("PD", PD), ("above switch", ABOVE), ("S", INS), ("outside S", ~INS)):
            L.append("      %-13s %3d games: " % (rn, m.sum()) + ", ".join(
                "%s %d" % (ORDER[i], (wa[m] == i).sum()) for i in range(7) if (wa[m] == i).any()))
        s111 = S[INS, I111]
        L.append("    111 share in S: median %.6g, min %.6g, max %.6g" % (np.median(s111), s111.min(), s111.max()))
        o = ~INS
        L.append("    111 share outside S: median %.3g, max %.3g" % (np.median(S[o, I111]), S[o, I111].max()))
        E = d["E"]
        L.append("    efficiency: min %.4f, mean %.4f, median %.4f, max %.4f; PD median %.4f" % (
            E.min(), E.mean(), np.median(E), E.max(), np.median(E[PD])))
        if R is not None:
            lead = R.argmax(2)
            nag = (lead == wa[:, None]).sum(1)
            L.append("    replicates naming the pooled leader, games by count: " + ", ".join(
                "%d/10: %d" % (k, (nag == k).sum()) for k in range(11) if (nag == k).any()))
            L.append("    games where not all ten agree, by quadrant: " + ", ".join(
                "%s %d" % (q, ((nag < 10) & (QUAD == q)).sum()) for q in ("PD", "SD", "SH", "HA")))
    return L


BISTABLE = "/Users/martin/Documents/CL1/DiskM1WF/Data/bistable.tsv"


def bistable_counts():
    n = {}
    for l in open(BISTABLE):
        if l.strip() and not l.startswith("#"):
            n[l.split("\t")[0]] = n.get(l.split("\t")[0], 0) + 1
    return "%d F1 and %d F2 lines" % (n.get("F1", 0), n.get("F2", 0))


def main():
    V, DD = {}, {}
    for run, _, _ in RUNS:
        DD[run] = load(run)
        V[run] = {sp: rows_of(DD[run][sp]) for sp, _ in SPACES}
    tex = table(V)
    with open(TEX, "w") as fh:
        fh.write(tex)
    body = table(V, legend=False)
    with open(TEX_BODY, "w") as fh:
        fh.write(body)
    chk, bad = check(V)
    L = ["sitable4.py -- SI Table 4 of PRFull (the four strategy spaces in the two populations), computed from the data.",
         "Written: %s, %s and %s." % (os.path.relpath(TEX, os.path.dirname(HERE)), os.path.relpath(TEX_BODY, os.path.dirname(HERE)),
                                    os.path.relpath(OUT, os.path.dirname(HERE))),
         "",
         "REGIONS (games.dat tags, checked by fig45 against the signs of u, v, u+v-1, u-v-1): PD (u<0<v) %d, above the "
         "switch line (wedges N, E) %d, wedge S %d, outside S %d, all 512." % (PD.sum(), ABOVE.sum(), INS.sum(),
                                                                               (~INS).sum()),
         ""]
    for run, _, par in RUNS:
        L.append("#" * 110)
        L.append("RUN %s (%s, eps = 1e-4, 10 replicates, second half sampled)" % (run, par))
        L.append("  %-50s %12s %12s %12s %12s" % ("row", "bin m1", "m1 cube", "bin m2", "m2 cube"))
        for key, rl, kind, _ in ROWS:
            L.append("  %-50s %12s %12s %12s %12s" % (key + "  " + rl.replace("\\texttt{", "").replace("}", "")
                                                      .replace("$", ""),
                                                      *[txt_num(V[run][sp][key], kind) for sp, _ in SPACES]))
        L.append("")
        L += details(run, DD[run])
        L.append("")
    L.append("#" * 110)
    L.append("HOW EACH ROW IS COMPUTED")
    for key, rl, kind, how in ROWS:
        L.append("  %-12s %s" % (key, how))
    L.append("  leading atom: fig4.winners (binary), fig4c.continuous (memory-one cube), figwfc.winners_c (memory-two cube)")
    L.append("  shares PS: fig45 (sums of the atom shares S over the atoms of each property, fig5.PROPS)")
    L.append("  table formatting: counts as integers; efficiencies to 3 decimals; shares to 3 decimals, below 0.01 to 2 "
             "significant digits")
    L.append("")
    L.append("#" * 110)
    L.append("CROSS-CHECKS against the fact sheets (a stated value is compared after rounding the computed one to its "
             "decimals)")
    L += chk
    L.append("")
    L.append("DISAGREEMENTS with the fact sheets: %d" % len(bad))
    L += bad
    L.append("")
    schk, sbad = stats_check(V)
    L.append("#" * 110)
    L.append("CROSS-CHECK against fig45_stats.txt (fig45.census over the same loaders; the text behind Figures 4, 5 and SI "
             "Figures 1, 2), at its printed precision")
    L += schk
    L.append("DISAGREEMENTS with fig45_stats.txt: %d" % len(sbad))
    L += sbad
    L.append("")
    L.append("#" * 110)
    L.append("NOTES")
    x = V["F1"]["c1"]["mean_eff"]
    L.append("  * memory-one cube, F1, mean efficient share: computed %.6f. facts_m1.md:1051 states 0.865 (agrees) and "
             "facts_m1.md:1059-1060 (the M1f README table, RM:97-103) 0.87, which is %.6f rounded twice (to 0.865, then "
             "0.87); at two decimals it is %.2f." % (x, x, x))
    L.append("  * binary memory one, replicate row: no per-replicate atom shares are available locally. wfdata.load('m1', "
             "run) returns pooled shares only; DiskM1WF/Data/mw1_e4.rep.tsv and mw2_e4.rep.tsv carry per replicate only "
             "pay, ef and ps1..ps4, the population's mean cooperation probability after each outcome (DiskM1WF/Cannon/mw/"
             "sf1.f:403-426), four marginals of sixteen shares. The w-files with each replicate's ten most abundant "
             "strategies (DiskM1WF/PROJECT.md:442) are not in the local tree. The binary paper states only a "
             "different criterion: at 59 (N = 100) and 11 (N = 1000) games the ten replicates fixed in different basins "
             "(some atom's between-replicate s.d. > 0.05; the binary paper [PR2026]; DiskM1WF/Data/"
             "bistable.tsv, %s)." % bistable_counts())
    L.append("  * binary memory two, replicate row: 487 (F1) and 412 (F2) reproduce the binary paper [PR2026] and "
             "facts_ref.md:726; they come from the seeded re-runs f1_e4 / f2_e4, not from wfdata, and are exact "
             "per-replicate atom shares (the .rep files), not the rank-1 proxy of the .win files.")
    L.append("  * the leading atom of every space is its loader's; the raw argmax of the pooled shares gives the same "
             "leader at every game in all eight runs (rawdiff = %s)." % ", ".join(
                 "%s %s %d" % (run, sp, DD[run][sp]["rawdiff"]) for run, _, _ in RUNS for sp, _ in SPACES))
    L.append("")
    L.append("#" * 110)
    L.append("THE TABLE, without its legend (%s)" % os.path.relpath(TEX_BODY, os.path.dirname(HERE)))
    L.append(body)
    L.append("#" * 110)
    L.append("LOADERS' LOG (fig45.load_binary / load_m1 / load_m2, wfrep.load)")
    L.append(LOADLOG.getvalue())
    with open(OUT, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L[:L.index("#" * 110, 6)]))
    for run, _, par in RUNS:
        print("%s (%s)" % (run, par))
        for key, rl, kind, _ in ROWS:
            print("  %-14s %s" % (key, "  ".join("%10s" % txt_num(V[run][sp][key], kind) for sp, _ in SPACES)))
    print("DISAGREEMENTS with the fact sheets: %d" % len(bad))
    for b in bad:
        print(b)
    print("DISAGREEMENTS with fig45_stats.txt: %d" % len(sbad))
    for b in sbad:
        print(b)
    print("wrote %s and %s" % (TEX, OUT))


if __name__ == "__main__":
    main()
