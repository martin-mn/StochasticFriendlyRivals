#!/usr/bin/env python3
"""analyze.py -- the tables of the wfc2 sensitivity campaign (Referee 1, issues 6-8), in the format of the paper's
section 'Evolution over the memory-two cube', and the comparison with production at the same games.  Local, seconds.

    python3 analyze.py [--data DIR] [ARM ...]    DIR: the collected out/ of the kit (default: ./out)
                                                 ARM: A1 A2 K0 E3 E5 E3F2 E5F2 B0 B0K0 TS (default: all with data)
    python3 analyze.py --selftest                production F1 read as if it were an arm: every change must be 0
    python3 analyze.py --probe [--data DIR]      the probe: bit-identical checks against the archive, and the costs

Read-only inputs (never written):
    production F1   PartnersRivalsM2f/check/wfc_a333  (r_1o3_*: REP lines; n_1o3_*, nr_1o3_*: the kernel masses)
    production F2   PartnersRivalsM2f/check/wfc_f2
    binary F1/F2    PartnersRivals/Data/f1_e4.rep, f2_e4.rep     (alpha = 1/2: uniform mutation over the binary genomes)
    frmass.json     this kit (frmass.py): the friendly-rival kernel mass per wedge at alpha = 1/3, 1/4, 2/5
Definitions are production's (Figures/figwfc.py, check/wfc_stats.py): the share of an atom is the mean over the
replicates of its sampled share, except 111 in wedge W, counted exactly; most abundant = argmax share over the atoms
with positive kernel mass; the enrichment is share / kernel mass; the most enriched atom is determined at 95 % when
the lower end of its interval (t interval of the replicate shares, 9 d.f., over the upper Poisson bound of the mass,
except 111, whose mass is fixed) lies above every other held atom's upper end; competitive = 001 + 011 + 111; partners = 110 + 111.
Writes analysis/<ARM>_games.tsv (one row per game) and analysis/summary.txt (everything printed).
"""
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M2 = "/Users/martin/Documents/CL1/PartnersRivalsM2f"
PROD = {"F1": M2 + "/check/wfc_a333", "F2": M2 + "/check/wfc_f2"}
BIN = {"F1": "/Users/martin/Documents/CL1/PartnersRivals/Data/f1_e4.rep",
       "F2": "/Users/martin/Documents/CL1/PartnersRivals/Data/f2_e4.rep"}
ORDER = ["000", "100", "010", "001", "110", "011", "111"]
IX = {c: i for i, c in enumerate(ORDER)}
BITS = ["%d%d%d" % (k >> 2, (k >> 1) & 1, k & 1) for k in range(8)]          # wfc's atom index -> code
COMP = [IX["001"], IX["011"], IX["111"]]
TAGS = [l.split()[7] for l in open(os.path.join(HERE, "games.dat")) if l.strip() and not l.startswith("#")]
WEDGE = np.array([t[-1] for t in TAGS]); QUAD = np.array([t[-3:-1] for t in TAGS])
FRM = json.load(open(os.path.join(HERE, "frmass.json")))
# arm -> (alpha tag, production regime, kernel masses: 'prod' = production's alpha = 1/3 draws, 'own' = the arm's)
ARMS = {"A1": ("1o10", "F1", "own"), "A2": ("1o5", "F1", "own"), "A3": ("2o5", "F1", "own"), "K0": ("1o3", "F1", "prod"),
        "E3": ("1o3", "F1", "prod"), "E5": ("1o3", "F1", "prod"),
        "E3F2": ("1o3", "F2", "prod"), "E5F2": ("1o3", "F2", "prod"),
        "B0": ("1o3", None, "prod"), "B0K0": ("1o3", None, "prod"), "TS": ("1o3", "F1", "prod")}
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)


# ------------------------------------------------------------------------------------------------ loading
def parse_rep(fn):
    """REP (and MUT) lines of one r_ file -> list of per-replicate dicts"""
    reps, mut = [], {}
    for l in open(fn):
        f = l.split()
        if not f:
            continue
        if f[0] == "MUT":
            if len(f) == 8:                             # the first K0 build's MUT line had no frwork field
                f = f[:7] + ["-1"] + f[7:]
            mut[int(f[1])] = dict(redraw=int(f[2]), step=int(f[3]), noop=int(f[4]), glob=int(f[5]),
                                  nfrscan=int(f[6]), frwork=int(f[7]), ex111u=int(f[8]))
        if f[0] != "REP":
            continue
        a = np.array([float(x) for x in f[12:20]]); ns = float(f[10])
        assert a[5] == 0, "atom 101 occupied in %s" % fn
        s = np.zeros(7)
        for k in range(8):
            if k != 5:
                s[IX[BITS[k]]] = a[k] / ns
        r = dict(ipt=int(f[2]), u=float(f[3]), v=float(f[4]), pay=float(f[6]), eff=float(f[7]), alpha=float(f[8]),
                 irep=int(f[9]), nsamp=ns, ndist=int(f[11]), s=s, samp111=s[IX["111"]], tsim=float(f[24]),
                 tcl=float(f[25]), inW=int(f[26]), ex111=float(f[27]), total=float(f[30]), line=l.rstrip("\n"))
        reps.append(r)
    for r in reps:
        r["mut"] = mut.get(r["irep"])
        if r["inW"]:                                    # wedge W: the exact 111 share replaces the sampled one ...
            u = r["mut"]["ex111u"] if r["mut"] else 0
            r["ex111u"] = u
            if u == 0:
                r["s"][IX["111"]] = r["ex111"] / r["total"]
            # ... unless K0's scan budget left some of it unresolved: then the sampled share stays, bounds are kept
            r["b111"] = (r["ex111"] / r["total"], (r["ex111"] + u) / r["total"])
    return reps


def load_arm(d, at):
    R = {}
    for fn in sorted(glob.glob(os.path.join(d, "r_%s_*.txt" % at))):
        for r in parse_rep(fn):
            R.setdefault(r["ipt"], []).append(r)
    for i in R:
        R[i].sort(key=lambda r: r["irep"])
        assert len({r["irep"] for r in R[i]}) == len(R[i]), "duplicate replicate at ipt %d in %s" % (i, d)
    return R


def load_mass(dirs, at, ipts):
    """kernel masses M (ipt -> 7) and hit counts H (inf for 111, whose mass is the census/sweep value)"""
    C, md = {}, {}
    for d in dirs:
        for fn in glob.glob(os.path.join(d, "n_%s_*.txt" % at)) + glob.glob(os.path.join(d, "nr_%s_*.txt" % at)):
            f = open(fn).read().split()
            assert f[0] == "NEU"
            i = int(f[2]); c = [float(x) for x in f[8:16]]
            assert c[5] == 0
            C[i] = C.get(i, np.zeros(7)) + np.array([c[BITS.index(code)] for code in ORDER]); md[i] = md.get(i, 0) + float(f[7])
    M, H, MD = {}, {}, {}
    for i in ipts:
        if i not in C:
            continue
        M[i] = C[i] / md[i]; M[i][IX["111"]] = FRM[at][WEDGE[i]]
        H[i] = C[i].copy(); H[i][IX["111"]] = np.inf; MD[i] = md[i]
    return M, H, MD


def load_binary(run):
    B = {}
    if not os.path.exists(BIN[run]):
        return B
    for l in open(BIN[run]):
        if l.startswith("#") or not l.strip():
            continue
        f = l.split(); B.setdefault(int(f[0]) - 2000, []).append((float(f[2]), np.array([float(x) for x in f[3:10]])))
    return B


# ------------------------------------------------------------------------------------------------ statistics
def game_stats(reps, M, H, md):
    """per-game statistics in production's definitions"""
    from scipy.stats import chi2, t as tdist
    E = np.array([r["eff"] for r in reps]); Sr = np.array([r["s"] for r in reps]); n = len(reps)
    S = Sr.mean(0)
    ok = M > 0 if M is not None else np.ones(7, bool)
    share = np.where(ok, S, -np.inf); lead = int(np.argmax(share))
    rl = [int(np.argmax(np.where(ok, x, -np.inf))) for x in Sr]
    g = dict(n=n, eff=E.mean(), effsd=E.std(ddof=1) if n > 1 else 0.0, S=S, Sr=Sr, E=E, lead=lead,
             agree=sum(x == lead for x in rl), comp=S[COMP].sum(), fr=S[IX["111"]], part=S[IX["110"]] + S[IX["111"]])
    if any("b111" in r for r in reps):
        g["b111"] = (np.mean([r["b111"][0] for r in reps]), np.mean([r["b111"][1] for r in reps]))
        g["unres"] = sum(r.get("ex111u", 0) > 0 for r in reps)
    if M is None:
        return g
    enr = np.where(ok, S / np.where(ok, M, 1.0), -np.inf); we = int(np.argmax(enr))
    g.update(enr=enr, we=we, fr_enr=enr[IX["111"]])
    # determination at 95 % (figwfc.undetermined_rep)
    lo = lambda k: chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0
    hi = lambda k: chi2.ppf(0.975, 2 * (k + 1)) / 2
    q = tdist.ppf(0.975, n - 1) if n > 1 else np.inf
    half = q * Sr.std(0, ddof=1) / math.sqrt(n) if n > 1 else np.full(7, np.inf)
    s_lo, s_hi = np.maximum(S - half, 0.0), S + half
    m_lo = np.array([M[j] if np.isinf(H[j]) else lo(H[j]) / md for j in range(7)])
    m_hi = np.array([M[j] if np.isinf(H[j]) else hi(H[j]) / md for j in range(7)])
    held = S > 0
    e_lo = np.where(held, s_lo / np.maximum(m_hi, 1e-300), 0.0)
    e_hi = np.where(held, np.where(m_lo > 0, s_hi / np.maximum(m_lo, 1e-300), np.inf), 0.0)
    g["det"] = not any(j != we and held[j] and e_hi[j] >= e_lo[we] for j in range(7))
    g["enr_lo"], g["enr_hi"] = e_lo, e_hi
    g["s_lo"], g["s_hi"] = s_lo, s_hi
    return g


def boot_median(games, key, nb=2000, seed=1):
    """95 % interval of the median over games of a per-game share, resampling replicates within each game"""
    rng = np.random.default_rng(seed); vals = []
    idx = IX["111"] if key == "fr" else None
    for _ in range(nb):
        v = []
        for g in games:
            Sr = g["Sr"]; pick = Sr[rng.integers(0, len(Sr), len(Sr))].mean(0)
            v.append(pick[idx] if idx is not None else pick[COMP].sum())
        vals.append(np.median(v))
    return np.percentile(vals, [2.5, 97.5])


def welch(a, b):
    if len(a) < 2 or len(b) < 2:
        return np.nan
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    return (a.mean() - b.mean()) / math.sqrt(va + vb) if va + vb > 0 else (0.0 if a.mean() == b.mean() else np.inf)


# ------------------------------------------------------------------------------------------------ report of one arm
def wedge_table(G, ipts, name):
    say("  %-4s %4s %8s %6s %6s %-14s %10s %10s %9s %8s %-22s %6s" %
        ("wedge", "n", "eff med", "111 ld", "110 ld", "other leads", "111 med", "comp med", "111 enr", "det", "most enriched (det)", "agree"))
    rows = {}
    for w in ["W", "S", "N", "E", "all"]:
        sel = [i for i in ipts if (w == "all" or WEDGE[i] == w)]
        if not sel:
            continue
        gs = [G[i] for i in sel]
        led = [ORDER[g["lead"]] for g in gs]
        oth = {}
        for c in led:
            if c not in ("111", "110"):
                oth[c] = oth.get(c, 0) + 1
        det = [g for g in gs if g.get("det")]
        wins = {}
        for g in det:
            wins[ORDER[g["we"]]] = wins.get(ORDER[g["we"]], 0) + 1
        fe = [g["fr_enr"] for g in gs if "fr_enr" in g]
        rows[w] = dict(n=len(sel), eff=np.median([g["eff"] for g in gs]), l111=led.count("111"), l110=led.count("110"),
                       fr=np.median([g["fr"] for g in gs]), comp=np.median([g["comp"] for g in gs]),
                       enr=np.median(fe) if fe else np.nan, det=len(det), wins=wins,
                       agree=sum(g["agree"] == g["n"] for g in gs))
        r = rows[w]
        say("  %-4s %4d %8.4f %6d %6d %-14s %10.3g %10.3g %9.3g %8d %-22s %6d" %
            (w, r["n"], r["eff"], r["l111"], r["l110"], ",".join("%s:%d" % kv for kv in sorted(oth.items())) or "-",
             r["fr"], r["comp"], r["enr"], r["det"], ",".join("%s:%d" % kv for kv in sorted(wins.items())) or "-", r["agree"]))
    return rows


def write_tsv(arm, G, P, ipts, B):
    os.makedirs(os.path.join(HERE, "analysis"), exist_ok=True)
    fn = os.path.join(HERE, "analysis", "%s_games.tsv" % arm)
    with open(fn, "w") as fo:
        fo.write("\t".join(["task", "ipt", "tag", "wedge", "nrep", "eff", "eff_sd", "lead", "agree", "fr_share",
                            "comp_share", "most_enriched", "determined", "fr_enr"] + ["s" + c for c in ORDER] +
                           ["prod_eff", "prod_lead", "prod_fr", "prod_comp", "welch_t_eff", "bin_eff", "bin_lead"]) + "\n")
        for i in ipts:
            g = G[i]; p = P.get(i)
            b = B.get(i)
            row = [i + 1, i, TAGS[i], WEDGE[i], g["n"], "%.6f" % g["eff"], "%.2e" % g["effsd"], ORDER[g["lead"]], g["agree"],
                   "%.4g" % g["fr"], "%.4g" % g["comp"], ORDER[g["we"]] if "we" in g else "", int(g.get("det", 0)),
                   "%.4g" % g.get("fr_enr", np.nan)] + ["%.4g" % x for x in g["S"]]
            if p:
                row += ["%.6f" % p["eff"], ORDER[p["lead"]], "%.4g" % p["fr"], "%.4g" % p["comp"], "%.2f" % welch(g["E"], p["E"])]
            else:
                row += [""] * 5
            if b:
                row += ["%.6f" % np.mean([x[0] for x in b]), ORDER[int(np.argmax(np.mean([x[1] for x in b], 0)))]]
            else:
                row += ["", ""]
            fo.write("\t".join(map(str, row)) + "\n")
    return fn


def report_arm(arm, data, prod_override=None):
    at, run, mass = ARMS.get(arm, ("1o3", "F1", "prod"))
    d = os.path.join(data, arm)
    R = load_arm(d, at)
    if not R:
        say("%s: no data in %s" % (arm, d)); return None
    ipts = sorted(R)
    nrep = sorted({len(v) for v in R.values()})
    say("")
    say("=" * 118)
    say("%s: %d games, replicates per game %s, alpha %s" % (arm, len(ipts), nrep, at.replace("o", "/")))
    hdr = [l for l in open(glob.glob(os.path.join(d, "r_%s_*.txt" % at))[0]) if l.startswith("#")]
    for h in hdr:
        say("  " + h.strip())
    mdirs = [d] if mass == "own" else [PROD["F1"]]
    M, H, MD = load_mass(mdirs, at, ipts)
    miss = [i for i in ipts if i not in M]
    if miss:
        say("  WARNING: no kernel masses at ipt %s: enrichments there are skipped" % miss)
    G = {i: game_stats(R[i], M.get(i), H.get(i), MD.get(i)) for i in ipts}
    if arm in ("B0", "B0K0"):
        return report_neutral(arm, R, G, M, MD, ipts)
    if arm == "TS":
        report_ts(data, R, G, ipts)
    # production at the same games
    P = {}
    if run:
        pdir = prod_override or PROD[run]
        PR = load_arm(pdir, "1o3")
        PM, PH, PMD = load_mass([PROD["F1"]], "1o3", ipts)
        P = {i: game_stats(PR[i], PM.get(i), PH.get(i), PMD.get(i)) for i in ipts if i in PR}
    B = load_binary(run) if run else {}
    say("  --- %s" % arm)
    rows = wedge_table(G, ipts, arm)
    if P:
        say("  --- production (%s, alpha = 1/3, global kernel, eps = 1e-4) at the same games" % run)
        prow = wedge_table(P, ipts, "prod")
    if B:
        say("  --- binary (alpha = 1/2, PR2026 %s) at the same games: median efficiency by wedge %s" %
            (run, ", ".join("%s %.4f" % (w, np.median([np.mean([x[0] for x in B[i]]) for i in ipts if WEDGE[i] == w]))
                            for w in "WSNE" if any(WEDGE[i] == w for i in ipts))))
    if any(r["mut"] for i in ipts for r in R[i]) and any("b111" in G[i] for i in ipts):
        un = [(i, G[i]["unres"]) for i in ipts if G[i].get("unres")]
        say("  W-111 exact counter: unresolved (K0 scan budget) at %d games%s" %
            (len(un), (": " + " ".join("ipt%d(%d reps)" % x for x in un)) if un else "; exact everywhere"))
    if P:
        ch = [i for i in ipts if i in P and G[i]["lead"] != P[i]["lead"]]
        t = np.array([welch(G[i]["E"], P[i]["E"]) for i in ipts if i in P])
        say("  vs production: most abundant atom unchanged at %d of %d games%s" %
            (len(ipts) - len(ch), len(ipts), ("; changed: " + ", ".join("ipt%d %s %s->%s (agree %d/%d, prod %d/%d)" %
             (i, TAGS[i][-3:], ORDER[P[i]["lead"]], ORDER[G[i]["lead"]], G[i]["agree"], G[i]["n"], P[i]["agree"], P[i]["n"])
             for i in ch)) if ch else ""))
        say("  vs production: efficiency, median change %+.2e, |Welch t| > 2 at %d of %d games (max %.1f)" %
            (np.median([G[i]["eff"] - P[i]["eff"] for i in ipts if i in P]), (np.abs(t) > 2).sum(), len(t), np.nanmax(np.abs(t)) if np.isfinite(t).any() else np.nan))
        for w in "WNE":
            gs = [P[i] for i in ipts if WEDGE[i] == w]
            if gs:
                bf, bc = boot_median(gs, "fr"), boot_median(gs, "comp")
                a_fr = np.median([G[i]["fr"] for i in ipts if WEDGE[i] == w]); a_c = np.median([G[i]["comp"] for i in ipts if WEDGE[i] == w])
                say("  wedge %s medians: 111 %.3g (production %.3g, 95%% [%.3g, %.3g]) %s; competitive %.3g (production %.3g, [%.3g, %.3g]) %s"
                    % (w, a_fr, prow[w]["fr"], bf[0], bf[1], "inside" if bf[0] <= a_fr <= bf[1] else "OUTSIDE",
                       a_c, prow[w]["comp"], bc[0], bc[1], "inside" if bc[0] <= a_c <= bc[1] else "OUTSIDE"))
    criteria(arm, G, P, B, ipts, rows)
    fn = write_tsv(arm, G, P, ipts, B)
    say("  per-game table: %s" % fn)
    cost(R)
    return G


def cost(R):
    ts = [r["tsim"] for v in R.values() for r in v]; tc = [r["tcl"] for v in R.values() for r in v]
    say("  cost: t_sim median %.0f s (max %.0f), t_class median %.1f s (max %.0f); %.1f core-h in REP lines" %
        (np.median(ts), max(ts), np.median(tc), max(tc), (sum(ts) + sum(tc)) / 3600))


def criteria(arm, G, P, B, ipts, rows):
    """the success criteria of the triage (reports 5 and 6), evaluated on the data"""
    S_ = [i for i in ipts if WEDGE[i] == "S"]
    out = []
    if arm in ("A1", "A2", "A3", "K0"):
        n111 = sum(bool(G[i]["lead"] == IX["111"] and G[i].get("det") and G[i]["we"] == IX["111"]) for i in S_)
        need = math.ceil(15 * len(S_) / 16)
        out.append(("111 leads and is the determined most enriched atom at >= %d of %d S games" % (need, len(S_)), n111, n111 >= need))
        for w, b in (("N", 30.0), ("E", 170.0)):
            e = rows.get(w, {}).get("enr", np.nan)
            out.append(("median 111 enrichment in %s at least 10x below binary %.0f" % (w, b), "%.3g" % e, e <= b / 10))
        det = [i for i in ipts if WEDGE[i] != "S" and G[i].get("det")]
        n110 = sum(G[i]["we"] == IX["110"] for i in det)
        out.append(("110 is the most enriched atom outside S where determined", "%d of %d" % (n110, len(det)), n110 == len(det)))
        if B:
            be = np.median([np.mean([x[0] for x in B[i]]) for i in ipts if i in B])
            out.append(("median efficiency at least binary's at the same games (%.4f)" % be, "%.4f" % rows["all"]["eff"], rows["all"]["eff"] >= be))
    if arm in ("E3", "E5", "E3F2", "E5F2") and P:
        same = sum(G[i]["lead"] == P[i]["lead"] for i in ipts if i in P)
        out.append(("most abundant atom unchanged at >= 95% of games", "%d of %d" % (same, len(P)), same >= 0.95 * len(P)))
        t = [abs(welch(G[i]["E"], P[i]["E"])) for i in ipts if i in P]
        nt = sum(bool(x > 2.26) for x in t)    # t_{0.975} at 9 d.f.
        out.append(("efficiency within replicate spread (|Welch t| <= 2.26) at >= 95% of games", "%d exceed" % nt, nt <= 0.05 * len(t)))
        if arm in ("E3", "E5"):
            n = sum(G[i]["lead"] == IX["111"] for i in S_)
            out.append(("friendly rivals lead every S game (N = 1000)", "%d of %d" % (n, len(S_)), n == len(S_)))
        else:
            mv = [(i, ORDER[P[i]["lead"]], ORDER[G[i]["lead"]]) for i in S_ if G[i]["lead"] != P[i]["lead"]]
            out.append(("S games whose leading atom moved (reported, as Revised2 did at N = 100)", " ".join("ipt%d %s->%s" % m for m in mv) or "none", None))
    for name, val, ok in out:
        say("  CRITERION %-86s %-14s %s" % (name, val, "" if ok is None else ("PASS" if ok else "FAIL")))


# ------------------------------------------------------------------------------------------------ B0: the neutral check
def report_neutral(arm, R, G, M, MD, ipts):
    """at beta = 0 the expected share of every atom is its kernel mass (both kernels leave the alpha-law invariant).
    Per (game, atom): z = (mean replicate share - mass) / sqrt(se_share^2 + se_mass^2), se_mass binomial in the kernel
    draws (0 for 111, whose mass is the census/sweep value), se_share the replicate s.e. floored at the binomial s.e.
    of the game's pooled samples (so that an unsampled rare atom gives a finite z).  The triage's literal criterion (every pair inside) is
    printed together with the multiplicity-aware one: the number of |z| > t_{0.975, n-1} must not exceed the 95 %
    quantile of Binomial(pairs, 0.05), and no |z| may exceed 4."""
    from scipy.stats import t as tdist, binom
    say("  neutral check: z = (mean share - kernel mass) / combined s.e., flagged beyond the t quantile (n-1 d.f.)")
    say("  %-5s %-5s %s" % ("ipt", "tag", "  ".join("%-24s" % ("%s share/mass z" % c) for c in ORDER)))
    npair = nout = 0; zmax = 0.0; vr = []
    for i in ipts:
        g = G[i]; Sr = g["Sr"]; n = len(Sr); q = tdist.ppf(0.975, n - 1) if n > 1 else np.inf
        m = M.get(i); md = MD.get(i, 0)
        nstot = sum(r["nsamp"] for r in R[i])
        cells = []
        for j in range(7):
            mu = Sr[:, j].mean(); ses = Sr[:, j].std(ddof=1) / math.sqrt(n) if n > 1 else np.inf
            # floor: the binomial s.e. of all the game's samples.  Without it an atom of small mass (111 in N and E,
            # ~1e-4, i.e. ~4 expected hits in 10 x 4000 samples) that no replicate happened to sample has ses = 0 and,
            # with sem = 0 for 111, z = inf: a certain FAIL of the multiplicity-aware criterion for a null outcome.
            ses = max(ses, math.sqrt(m[j] * (1 - m[j]) / nstot))
            sem = 0.0 if j == IX["111"] else math.sqrt(m[j] * (1 - m[j]) / md)
            if m[j] == 0 and mu == 0:
                cells.append("%-24s" % "absent"); continue
            z = (mu - m[j]) / math.sqrt(ses ** 2 + sem ** 2) if ses ** 2 + sem ** 2 > 0 else np.inf
            npair += 1; nout += abs(z) > q; zmax = max(zmax, abs(z))
            cells.append("%-24s" % ("%.3g/%.3g %+.1f%s" % (mu, m[j], z, "*" if abs(z) > q else "")))
            if 0.01 < m[j] < 0.99 and j != IX["111"] and n > 1:
                vr.append(Sr[:, j].var(ddof=1) / (m[j] * (1 - m[j]) / R[i][0]["nsamp"]))
        say("  %-5d %-5s %s" % (i, TAGS[i][-3:], "  ".join(cells)))
    lim = int(binom.ppf(0.95, npair, 0.05)) if npair else 0
    say("  CRITERION (triage, literal) every atom's interval contains its kernel mass at every game: %d of %d pairs outside  %s"
        % (nout, npair, "PASS" if nout == 0 else "FAIL"))
    say("  CRITERION (multiplicity-aware) outside <= %d (95%% of Binomial(%d, 0.05)) and max |z| <= 4: %d, max |z| %.1f  %s"
        % (lim, npair, nout, zmax, "PASS" if nout <= lim and zmax <= 4 else "FAIL"))
    if vr:
        say("  neutral replicate variance / binomial (atoms with mass in (1%%, 99%%)): median %.1f, p90 %.1f (production F1 at beta = 100: 4.4)"
            % (np.median(vr), np.percentile(vr, 90)))
    cost(R)
    return G


# ------------------------------------------------------------------------------------------------ TS: the diagnostics
def geyer_iat(x):
    """integrated autocorrelation time (in samples) by Geyer's initial positive sequence"""
    x = np.asarray(x, float); n = len(x)
    if n < 4 or x.var() == 0:
        return 1.0
    y = x - x.mean(); f = np.fft.rfft(y, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / (n * y.var())
    tau = -1.0
    for k in range(0, n - 1, 2):
        s = ac[k] + ac[k + 1]
        if s <= 0:
            break
        tau += 2 * s
    return max(tau, 1.0)


def parse_ts(fn):
    TB, TS, TSR = {}, {}, {}
    for l in open(fn):
        f = l.split()
        if not f:
            continue
        if f[0] == "TB":
            TB.setdefault(int(f[1]), []).append((int(f[3]), int(f[4]), float(f[6]), int(f[7]), float(f[8])))
        elif f[0] == "TS":
            TS.setdefault(int(f[1]), []).append((int(f[3]), int(f[4]), np.array([int(x) for x in f[5:13]])))
        elif f[0] == "TSR":
            TSR[int(f[1])] = f
    return TB, TS, TSR


def report_ts(data, R, G, ipts):
    """time-series diagnostics (report 6): windows P = [itend/8, itend/4] (production's sampled half, at 4x length) and
    L = [itend/2, itend] (the REP line's half); quarters Q3/Q4 = the two halves of P"""
    prod = load_arm(PROD["F1"], "1o3")
    split = {i for i in prod if len({int(np.argmax(r["s"])) for r in prod[i]}) > 1}
    say("  --- time-series diagnostics: P = production window at 4x length, L = the second half; Q3/Q4 = halves of P")
    say("  %-5s %-4s %-5s %-4s %8s %8s %8s %6s %7s %7s %7s %6s %6s %6s %s" %
        ("ipt", "tag", "split", "lead", "eff P", "eff L", "eff prd", "z Q34", "lead P", "lead L", "lead pr", "z Q34", "ESSmin", "ESSmed", "111 takeovers/escapes per 1e7 (exposure)"))
    fails = []; fails_pl = []
    for i in ipts:
        fns = glob.glob(os.path.join(data, "TS", "r_1o3_%d_*.txt" % (i + 1)))
        lead = int(np.argmax(np.mean([r["s"] for r in prod[i]], 0))) if i in prod else G[i]["lead"]
        rows = []
        for fn in fns:
            TB, TS, _ = parse_ts(fn)
            for irep in TB:
                tb = TB[irep]; ts = TS[irep]
                itend = tb[-1][1]
                P0, P1, L0 = itend // 8, itend // 4, itend // 2
                eP = [x[2] for x in tb if x[0] > P0 and x[1] <= P1]; eL = [x[2] for x in tb if x[0] > L0]
                q = len(eP) // 2; e3, e4 = np.array(eP[:q]), np.array(eP[q:])
                sP = [x for x in ts if P0 < x[0] <= P1]; sL = [x for x in ts if x[0] > L0]
                lp = np.array([x[2][[k for k in range(8) if BITS[k] == ORDER[lead]][0]] / x[1] for x in sP])
                ll = np.array([x[2][[k for k in range(8) if BITS[k] == ORDER[lead]][0]] / x[1] for x in sL])
                h = len(lp) // 2
                # batch means for the lead share: 10 batches per quarter
                bm = lambda a: np.array([b.mean() for b in np.array_split(a, 10)])
                b3, b4 = bm(lp[:h]), bm(lp[h:])
                # 111 state from the samples (exact blocks in W): smoothed over 5 samples, > 1/2 = friendly rivals hold
                f111 = np.array([x[2][7] / x[1] for x in ts]); tt = np.array([x[0] for x in ts])
                if tb[0][3] >= 0:                                   # wedge W: exact per block
                    f111 = np.array([x[3] / x[4] for x in tb]); tt = np.array([x[1] for x in tb])
                sm = np.convolve(f111, np.ones(5) / 5, mode="same") > 0.5
                keep = tt > itend // 40                            # skip the first 2.5 % (burn-in)
                st = sm[keep]
                up = int(((~st[:-1]) & st[1:]).sum()); dn = int((st[:-1] & (~st[1:])).sum())
                dt = (tt[1] - tt[0]) if len(tt) > 1 else 1
                rows.append(dict(eP=np.mean(eP), eL=np.mean(eL), z3=(e4.mean() - e3.mean()) / math.sqrt(e3.var(ddof=1) / len(e3) + e4.var(ddof=1) / len(e4) + 1e-300),
                                 lP=lp.mean(), lL=ll.mean(), zl=(b4.mean() - b3.mean()) / math.sqrt(b3.var(ddof=1) / 10 + b4.var(ddof=1) / 10 + 1e-300),
                                 ess=len(lp) / geyer_iat(lp), up=up, dn=dn, t0=(~st).sum() * dt, t1=st.sum() * dt,
                                 aP=sum(x[2] for x in sP) / max(1, sum(x[1] for x in sP)),
                                 aL=sum(x[2] for x in sL) / max(1, sum(x[1] for x in sL))))
        if not rows:
            continue
        n = len(rows)
        a = lambda k: np.array([r[k] for r in rows])
        zP = a("eP"); z_eff = (a("z3").mean() * math.sqrt(n)) if n > 1 else a("z3")[0]
        # pooled Q4 - Q3 over replicates: mean of the per-replicate z times sqrt(n) is ~N(0,1) under stationarity
        z_lead = a("zl").mean() * math.sqrt(n)
        pe = np.mean([r["eff"] for r in prod[i]]) if i in prod else np.nan
        pl = np.mean([r["s"][lead] for r in prod[i]]) if i in prod else np.nan
        up, dn, t0, t1 = a("up").sum(), a("dn").sum(), a("t0").sum(), a("t1").sum()
        from scipy.stats import chi2
        rate = lambda k, T: "%.2g [%.2g, %.2g]" % (k / T * 1e7, (chi2.ppf(0.025, 2 * k) / 2 if k else 0) / T * 1e7, chi2.ppf(0.975, 2 * k + 2) / 2 / T * 1e7) if T > 0 else "-"
        say("  %-5d %-4s %-5s %-4s %8.4f %8.4f %8.4f %6.1f %7.3f %7.3f %7.3f %6.1f %6.0f %6.0f up %s (%.1e) down %s (%.1e)" %
            (i, TAGS[i][-3:], "yes" if i in split else "no", ORDER[lead], zP.mean(), a("eL").mean(), pe, z_eff,
             a("lP").mean(), a("lL").mean(), pl, z_lead, a("ess").min(), np.median(a("ess")), rate(up, t0), t0, rate(dn, t1), t1))
        if i not in split and (abs(z_eff) > 2 or abs(z_lead) > 2 or a("ess").min() < 50):
            fails.append(i)
        # run length (report 6): the leading atom and the efficiency of the window P against the second half L
        ldP = int(np.argmax(np.mean([r["aP"] for r in rows], 0))); ldL = int(np.argmax(np.mean([r["aL"] for r in rows], 0)))
        tPL = welch(a("eL"), a("eP"))
        if ldP != ldL or (np.isfinite(tPL) and abs(tPL) > 2.26) or ldP != [k for k in range(8) if BITS[k] == ORDER[lead]][0]:
            say("        ipt %d: window P leads %s, L leads %s (production %s); efficiency L vs P Welch t %.2f" %
                (i, BITS[ldP], BITS[ldL], ORDER[lead], tPL))
        if i not in split and (ldP != ldL or (np.isfinite(tPL) and abs(tPL) > 2.26)):
            fails_pl.append(i)
    say("  CRITERION at every unsplit game the quarters agree within 2 SE and ESS >= 50 per replicate  %s" %
        ("PASS" if not fails else "FAIL at ipt " + " ".join(map(str, fails))))
    say("  CRITERION at every unsplit game the leading atom and the efficiency agree between P and L (|Welch t| <= 2.26)  %s" %
        ("PASS" if not fails_pl else "FAIL at ipt " + " ".join(map(str, fails_pl))))


# ------------------------------------------------------------------------------------------------ the probe
def probe(data):
    """DEF, DEF2, TSV must reproduce the archived production lines exactly (timings aside); then the costs"""
    strip = lambda l: " ".join(x for k, x in enumerate(l.split()) if k not in (24, 25))
    ref = [("DEF", PROD["F1"] + "/r_1o3_30_1.txt", "r_1o3_30_1.txt"), ("DEF2", PROD["F2"] + "/r_1o3_511_1.txt", "r_1o3_511_1.txt")]
    # TSV: every game of tasks_probe.txt's TSV lines (the 8 validation games), against production F1 replicate 1
    tsv = [int(l.split()[3]) for l in open(os.path.join(HERE, "tasks_probe.txt")) if l.startswith("probe/TSV ")]
    ref += [("TSV", PROD["F1"] + "/r_1o3_%d_1.txt" % g, "r_1o3_%d_1.txt" % g) for g in tsv]
    ok = True
    for arm, a, b in ref:
        fb = os.path.join(data, "probe", arm, b)
        if not os.path.exists(fb):
            say("PROBE %s: missing %s" % (arm, fb)); ok = False; continue
        ra = [l for l in open(a) if l.startswith("REP ") and l.split()[9] == "1"][0]
        rb = [l for l in open(fb) if l.startswith("REP ")]
        same = len(rb) == 1 and strip(ra) == strip(rb[0])
        ok &= same
        say("PROBE %-4s REP line %s the archived production line (%s)" % (arm, "IDENTICAL to" if same else "DIFFERS from", a))
    # the strategy dumps: replicate 1's records against the archived production dumps (ref/archived_rep1.sha256)
    import hashlib, struct
    href = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, "ref", "archived_rep1.sha256"))
            if l.strip() and not l.startswith("#")}
    for arm, fn, key in (("DEF", "S_1o3_30_1.bin", "F1"), ("TSV", "S_1o3_30_1.bin", "F1"), ("DEF2", "S_1o3_511_1.bin", "F2")):
        f = os.path.join(data, "probe", arm, fn)
        if not os.path.exists(f):
            continue
        b = open(f, "rb").read()
        r1 = b"".join(b[k:k + 140] for k in range(0, len(b), 140) if struct.unpack("<i", b[k:k + 4])[0] == 1)
        same = hashlib.sha256(r1).hexdigest() == href[key]; ok &= same
        say("PROBE %-4s strategy dump (replicate 1) %s the archived one" % (arm, "IDENTICAL to" if same else "DIFFERS from"))
    na, nb = PROD["F1"] + "/n_1o3_30.txt", os.path.join(data, "probe", "DEF", "n_1o3_30.txt")
    if os.path.exists(nb):
        s = lambda fn: open(fn).read().split()[:-1]
        same = s(na) == s(nb); ok &= same
        say("PROBE DEF  NEU line %s the archived one (%s)" % ("IDENTICAL to" if same else "DIFFERS from", na))
    for d in sorted(glob.glob(os.path.join(data, "probe", "*"))):
        arm = os.path.basename(d); at = ARMS.get(arm.rsplit("_", 1)[0] if arm[-2:] in ("_s", "_f") else arm, ("1o3",))[0]
        R = load_arm(d, at)
        for i in sorted(R):
            for r in R[i]:
                gen = int([l for l in open(glob.glob(os.path.join(d, "r_%s_%d_*.txt" % (at, i + 1)))[0]) if l.startswith("# wfc ")][0].split("itend=")[1].split()[0])
                mu = r["mut"]
                say("PROBE %-5s ipt %3d %-4s itend %9d  t_sim %7.1f  t_class %6.1f  ndist %5d  eff %.4f  lead %s%s" %
                    (arm, i, TAGS[i][-3:], gen, r["tsim"], r["tcl"], r["ndist"], r["eff"], ORDER[int(np.argmax(r["s"]))],
                     ("  MUT nfrscan %d ex111u %d" % (mu["nfrscan"], mu["ex111u"])) if mu else ""))
    say("PROBE bit-identical checks: %s" % ("ALL PASS" if ok else "FAILURES ABOVE"))
    return ok


def selftest():
    """production F1 read as an arm: the comparison machinery must report no change at all"""
    PR = load_arm(PROD["F1"], "1o3"); ipts = sorted(PR)
    M, H, MD = load_mass([PROD["F1"]], "1o3", ipts)
    G = {i: game_stats(PR[i], M[i], H[i], MD[i]) for i in ipts}
    lead = [ORDER[G[i]["lead"]] for i in ipts]
    und = sum(not G[i]["det"] for i in ipts)
    say("SELFTEST production F1 through analyze.py: %d games; most abundant %s; undetermined %d (paper: 61); "
        "median efficiency %.3f (paper 0.987); S games led by 111: %d of %d (paper 126)" %
        (len(ipts), {c: lead.count(c) for c in set(lead)}, und, np.median([G[i]["eff"] for i in ipts]),
         sum(G[i]["lead"] == IX["111"] for i in ipts if WEDGE[i] == "S"), sum(WEDGE[i] == "S" for i in ipts)))
    for w, paper in (("W", 8.0e-7), ("N", 1.25e-4), ("E", 2.25e-4)):
        say("SELFTEST wedge %s median 111 share %.3g (paper %.2g), competitive %.4f" %
            (w, np.median([G[i]["fr"] for i in ipts if WEDGE[i] == w]), paper, np.median([G[i]["comp"] for i in ipts if WEDGE[i] == w])))
    ok = und == 61 and abs(np.median([G[i]["eff"] for i in ipts]) - 0.987) < 5e-4
    say("SELFTEST %s" % ("PASS (reproduces the paper's numbers)" if ok else "FAIL"))
    return ok


def main(argv):
    data = os.path.join(HERE, "out")
    if "--data" in argv:
        k = argv.index("--data"); data = argv[k + 1]; argv = argv[:k] + argv[k + 2:]
    if "--selftest" in argv:
        ok = selftest()
    elif "--probe" in argv:
        ok = probe(data)
    else:
        arms = [a for a in argv if not a.startswith("-")] or [a for a in ARMS if os.path.isdir(os.path.join(data, a))]
        for a in arms:
            report_arm(a, data)
        ok = True
    os.makedirs(os.path.join(HERE, "analysis"), exist_ok=True)
    import datetime
    open(os.path.join(HERE, "analysis", "summary.txt"), "a").write(
        "#### analyze.py %s  (%s)\n" % (" ".join(sys.argv[1:]), datetime.datetime.now().isoformat(timespec="seconds"))
        + "\n".join(OUT) + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
