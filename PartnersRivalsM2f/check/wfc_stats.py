#!/usr/bin/env python3
"""wfc_stats.py -- every number the section 'What evolution selects' of Text/ms.tex quotes, from the data.

Three runs on the same 512 sunflower games, all N = 1000, beta = 100, mu = 1e-2, eps = 1e-4, 10 x 1e7 generations:
    m1c  the memory-one cube        PartnersRivalsM1f's run, read through its own cwdata.load() (read-only)
    m2b  the 65536 binary M2        PartnersRivals1's run F1, through its FinalFigures/wfdata.py (read-only)
    m2c  the memory-two cube        this project's cannon/wfc run 48366947 (+ 48391602), through Figures/figwfc.py
Quadrant and wedge of each game from the tag of games.dat (V_pNNN_<quadrant><wedge>).
Partners are the atoms 110 and 111, friendly rivals the atom 111, efficient atoms 100, 110, 111.
"""
import os
import sys

sys.dont_write_bytecode = True
import numpy as np

ROOT = "/Users/martin/Documents/CL1/PartnersRivalsM2f"
sys.path.insert(0, os.path.join(ROOT, "Figures"))
import figwfc as F                                                 # noqa: E402  (puts FinalFigures on sys.path)
import wfdata                                                      # noqa: E402
sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivalsM1f/Figures")
import cwdata                                                      # noqa: E402

# usage: wfc_stats.py [F1|F2]   F1 (default): N = 1000, the section's numbers; F2: N = 100, the small population
RUN = sys.argv[1] if len(sys.argv) > 1 else "F1"
CDIR = {"F1": "check/wfc_a333", "F2": "check/wfc_f2"}[RUN]
ORDER = wfdata.ORDER
IX = {c: ORDER.index(c) for c in ORDER}
tags = [l.split()[7] for l in open(os.path.join(ROOT, "cannon/wfc/games.dat")) if l.strip() and not l.startswith("#")]
ipt_of = {int(t.split("_")[1][1:]): t.split("_")[2] for t in tags}
QUAD = np.array([ipt_of[i][:2] for i in range(512)])
WEDGE = np.array([ipt_of[i][2] for i in range(512)])

# ---- the three runs: shares S (512,7) in ORDER, efficiency E, divisor for enrichment
Sc, Mc, Ec, nr, mdraw, UVc, Kc, Hc = F.load_continuum(os.path.join(ROOT, CDIR), nd=os.path.join(ROOT, "check/wfc_a333"))
Sb, Nb, Eb, UVb, Db = wfdata.load("m2", RUN)
d1 = cwdata.load(False, uv_ref=wfdata.load("m1", RUN)[3], run=RUN)
print("run %s: m1c %s, m2b %s, m2c %s" % (RUN, "PartnersRivalsM1f " + RUN, "PartnersRivals1 " + RUN, CDIR))
S1, K1, E1 = d1["S"], d1["K"], d1["E"]
assert np.allclose(UVb, UVc, atol=1e-9) and np.allclose(UVb, d1["UV"], atol=1e-9)
# wedge from the tags agrees with the one computed from (u, v)
assert all(F.wedge(*UVc[i]) == WEDGE[i] for i in range(512))

RUNS = {"m1c": (S1, E1), "m2b": (Sb, Eb), "m2c": (Sc, Ec)}
lead = {k: v[0].argmax(1) for k, v in RUNS.items()}               # most abundant atom (every held atom has mass > 0)
isP = lambda a: np.isin(a, [IX["110"], IX["111"]])
isFR = lambda a: a == IX["111"]
isEff = lambda a: np.isin(a, [IX["100"], IX["110"], IX["111"]])
comp = {k: v[0][:, [IX["001"], IX["011"], IX["111"]]].sum(1) for k, v in RUNS.items()}
effsh = {k: v[0][:, [IX["100"], IX["110"], IX["111"]]].sum(1) for k, v in RUNS.items()}
above = np.isin(WEDGE, ["N", "E"]); PD = QUAD == "PD"; SD = QUAD == "SD"; S_ = WEDGE == "S"

print("games: PD %d, SD %d, SH %d, HA %d; wedges W %d S %d N %d E %d; above the switch line %d"
      % (PD.sum(), SD.sum(), (QUAD == "SH").sum(), (QUAD == "HA").sum(), (WEDGE == "W").sum(), S_.sum(),
         (WEDGE == "N").sum(), (WEDGE == "E").sum(), above.sum()))
print()
print("%-44s %10s %10s %10s" % ("", "m1c", "m2b", "m2c"))
row = lambda name, f: print("%-44s %10s %10s %10s" % ((name,) + tuple(f(k) for k in ("m1c", "m2b", "m2c"))))
row("median efficiency", lambda k: "%.3f" % np.median(RUNS[k][1]))
row("min efficiency", lambda k: "%.3f" % RUNS[k][1].min())
row("games with efficiency < 0.9", lambda k: "%d" % (RUNS[k][1] < 0.9).sum())
row("games with efficiency < 0", lambda k: "%d" % (RUNS[k][1] < 0).sum())
row("PD: games where partners lead", lambda k: "%d" % isP(lead[k][PD]).sum())
row("PD: games where friendly rivals lead", lambda k: "%d" % isFR(lead[k][PD]).sum())
row("PD: median efficiency", lambda k: "%.3f" % np.median(RUNS[k][1][PD]))
row("above switch line: efficient atom leads", lambda k: "%d" % isEff(lead[k][above]).sum())
row("above switch line: partners lead", lambda k: "%d" % isP(lead[k][above]).sum())
row("above switch line: median efficient share", lambda k: "%.3f" % np.median(effsh[k][above]))
row("SD: games with efficiency < 0.9", lambda k: "%d" % (RUNS[k][1][SD] < 0.9).sum())
row("SD: median efficiency", lambda k: "%.3f" % np.median(RUNS[k][1][SD]))
row("wedge S: games where friendly rivals lead", lambda k: "%d" % isFR(lead[k][S_]).sum())
row("wedge S: median share of friendly rivals", lambda k: "%.3f" % np.median(RUNS[k][0][S_, IX["111"]]))
row("outside S: games where friendly rivals lead", lambda k: "%d" % isFR(lead[k][~S_]).sum())
for w in "WNE":
    m = WEDGE == w
    row("  wedge %s: friendly rivals lead" % w, lambda k: "%d" % isFR(lead[k][m]).sum())
    row("  wedge %s: median friendly-rival share" % w, lambda k: "%.2e" % np.median(RUNS[k][0][m, IX["111"]]))
row("median competitive share", lambda k: "%.4f" % np.median(comp[k]))
row("competitive share > 1/2 at", lambda k: "%d" % (comp[k] > 0.5).sum())
for w in "WNE":
    m = WEDGE == w
    row("  wedge %s: median competitive share" % w, lambda k: "%.4f" % np.median(comp[k][m]))
row("wedge S: median competitive share", lambda k: "%.4f" % np.median(comp[k][S_]))
print()
for k in ("m2b", "m2c"):
    print("%s most abundant: %s" % (k, ", ".join("%s %d" % (ORDER[i], (lead[k] == i).sum()) for i in range(7) if (lead[k] == i).any())))

# ---- replicate agreement on the leading atom, m2c: per-replicate shares from the REP lines
BITS = F.BITS
rep = {}
import glob
for fn in glob.glob(os.path.join(ROOT, CDIR, "r_1o3_*.txt")):
    for l in open(fn):
        f = l.split()
        if not f or f[0] != "REP":
            continue
        ipt = int(f[2]); ns = float(f[10]); a = [float(x) for x in f[12:20]]
        s = np.zeros(7)
        for kk in range(8):
            if kk != 5:
                s[IX[BITS[kk]]] = a[kk] / ns
        if int(f[26]):
            s[IX["111"]] = float(f[27]) / float(f[30])
        rep.setdefault(ipt, []).append(s.argmax())
agree = sum(all(r == lead["m2c"][i] for r in rep[i]) for i in range(512))
print("m2c: all ten replicates name the pooled leading atom at %d of 512 games" % agree)
dis = [i for i in range(512) if not all(r == lead["m2c"][i] for r in rep[i])]
print("   the others by quadrant: " + ", ".join("%s %d" % (q, sum(QUAD[i] == q for i in dis)) for q in ("PD", "SD", "SH", "HA")))
repb = {}
for l in open("/Users/martin/Documents/CL1/PartnersRivals/Data/%s_e4.rep" % RUN.lower()):
    if l.startswith("#") or not l.strip():
        continue
    f = l.split(); repb.setdefault(int(f[0]) - 2000, []).append(np.array([float(x) for x in f[3:10]]).argmax())
agb = sum(all(r == lead["m2b"][i] for r in repb[i]) for i in range(512))
print("m2b: all ten replicates name the pooled leading atom at %d of 512 games" % agb)

# ---- the most enriched atom of m2c, and what is not determined
ok = Mc > 0
enr = np.where(ok, Sc / np.where(ok, Mc, 1.0), -np.inf); we = enr.argmax(1)
und = F.undetermined_rep(F.load_reps(os.path.join(ROOT, CDIR)), Mc, Hc, mdraw, S=Sc)   # the ten replicates, as the figure
und_p = F.undetermined(Sc, Mc, Kc, Hc, mdraw)                                            # Poisson in the pooled samples
print("m2c most enriched: %s; not determined at 95%% at %d games (Poisson in the pooled samples: %d)"
      % (", ".join("%s %d" % (ORDER[i], (we == i).sum()) for i in range(7) if (we == i).any()), len(und), len(und_p)))
ud = {i for i, _, _ in und}
det = [i for i in range(512) if i not in ud]
print("   determined: %s; 111 there by wedge: %s; undetermined by wedge: %s"
      % (", ".join("%s %d" % (ORDER[j], sum(we[i] == j for i in det)) for j in range(7) if any(we[i] == j for i in det)),
         ", ".join("%s %d" % (w, sum(WEDGE[i] == w and we[i] == IX["111"] for i in det)) for w in "WSNE"),
         ", ".join("%s %d" % (w, sum(WEDGE[i] == w for i in ud)) for w in "WSNE")))
print("m2c neutral draws: %s" % ", ".join("%d games x %d" % (c, v) for v, c in zip(*np.unique(mdraw, return_counts=True))))

# ---- cost, from the REP and NEU lines themselves (the connection to Cannon is not needed)
ts = tc = 0.0
for fn in glob.glob(os.path.join(ROOT, CDIR, "r_1o3_*.txt")):
    for l in open(fn):
        f = l.split()
        if f and f[0] == "REP":
            ts += float(f[24]); tc += float(f[25])
tn = sum(float(open(fn).read().split()[18]) for fn in glob.glob(os.path.join(ROOT, CDIR, "n*_1o3_*.txt")))
print("m2c cost: simulation %.0f, classification %.0f, neutral draws %.0f core-hours; total %.0f"
      % (ts / 3600, tc / 3600, tn / 3600, (ts + tc + tn) / 3600))
