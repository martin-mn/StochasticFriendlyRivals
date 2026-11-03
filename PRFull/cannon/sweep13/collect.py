#!/usr/bin/env python3
r"""collect.py -- SI Table 3 of PRFull/Text/ms.tex (and the rarity numbers of the main text and Methods) from the
task outputs of a rarity sweep, `census psweep`.

    python3 collect.py SOURCE [SOURCE ...] [options]

SOURCE is a directory of task outputs (a_<tag>.tsv, with a_<tag>.log or a_<tag>.err, pt_<tag>.txt, and the Slurm
j_*.err), or a file of ALPHA lines such as PartnersRivalsM2f/check/sweep_a333.tsv.  One ALPHA line per task:

    ALPHA  alpha  n  effCC effALT effZ rivP rivM  FR_W FR_S FR_N FR_E  rayP rayM

(alpha is the tag given to census: "333" permille, or "1/3" exact.)

Options
  --alpha TAG       every line must carry this alpha tag (default: the lines must agree among themselves)
  --tasks N         expected number of tasks (default 512)
  --n N             expected draws per task (default 3190000)
  --seedbase S      directory mode, from the task logs: task i must have used seed S + i
  --glob G          directory mode: the task files (default a_*.tsv; the archived 0.333 directory: 'a_333_*.tsv')
  --errglob G       directory mode: the Slurm stderr files that must be empty (default j_*.err and probe_*.err;
                    with --jobid, only that array's)
  --jobid FILE      first, over the open ssh session: `sacct -j <id> -X` must show every task COMPLETED (and
                    `squeue -j <id>` nothing), or collect.py refuses to report (~/.claude/CLAUDE.md protocol)
  --fetch REMOTE    with --jobid, after that check: rsync cannon:REMOTE/ (a remote output directory, e.g.
                    CL1/PRFull/cannon/sweep13/out) into the first SOURCE directory
  --probe           SOURCE is the probe's output directory (out_probe): the bit-identity checks against the
                    archived probe 45830542 task 4 and the cost at alpha = 1/3, sized for the full run
  --boot B          bootstrap resamples over tasks for the ratios (default 20000, 0 = none)
  --compare SOURCE  a second sweep (e.g. the alpha = 0.333 file): each rate here / there, with a 95% interval
  --json FILE       every number, machine-readable
  --reproduce       run on the archived alpha = 0.333 sweep and check every published SI Table 3 number

Statistics (as in the Ref. 1 triage, scratchpad rarity_ci.py): exact Clopper-Pearson 95% intervals for rates;
the ratio observed / independent R = p11 / (p1 p2), where a friendly rival is exactly (efficient AND competitive),
with the multinomial delta method on the log scale,
    Var(log R) = [1/p11 - 1/p1 - 1/p2 + 2 p11/(p1 p2) - 1] / n,
and a nonparametric bootstrap over the tasks (numpy default_rng(1)); the index of dispersion of every count over
the tasks (chi-square, df = tasks - 1); the N-E and rivP-rivM mirror checks.
"""
import argparse, glob, json, math, os, re, subprocess, sys
import numpy as np
from scipy.stats import beta, chi2, norm

HERE = os.path.dirname(os.path.abspath(__file__))
A333 = "/Users/martin/Documents/CL1/PartnersRivalsM2f/check/sweep_a333.tsv"
COLS = ["effCC", "effALT", "effZ", "rivP", "rivM", "FR_W", "FR_S", "FR_N", "FR_E", "rayP", "rayM"]
ORDER = ["W", "S", "N", "E"]                             # the row order of SI Table 3
MARG = {"W": ("effCC", "rivP"), "S": ("effCC", "rivM"), "N": ("effALT", "rivP"), "E": ("effALT", "rivM")}
Z = norm.ppf(0.975)
# Methods, 'Rarity': the census-predicted change of every rate from alpha = 0.333 to exactly 1/3
PRED_SHIFT = {"W": 0.0088, "S": 0.0012, "N": 0.0069, "E": 0.0069}


def die(msg):
    sys.exit("collect.py: " + msg)


# ------------------------------------------------------------------ reading
def parse_alpha_line(line, where):
    f = line.rstrip("\n").split("\t")
    if len(f) != 14 or f[0] != "ALPHA":
        die("%s: not an ALPHA line: %r" % (where, line))
    return f[1], int(f[2]), [int(x) for x in f[3:]]


SUMMARY = re.compile(r"psweep alpha=(\S+) n=(\d+): effCC (\d+) effALT (\d+) effZ (\d+) rivP (\d+) rivM (\d+) \| "
                     r"FR_W (\d+) FR_S (\d+) FR_N (\d+) FR_E (\d+) \| rayP (\d+) rayM (\d+) \| meanpin ([0-9.]+)$")


def read_dir(d, problems, seedbase, tglob="a_*.tsv", eglob=("j_*.err", "probe_*.err"), unique_seeds=True):
    """One task per a_<tag>.tsv.  Checks each against its census summary line, its W-point file and (if present)
    the task header written by task.sh; returns [(tag, alpha, n, counts, meta)]."""
    rows = []
    def natural(p):             # a_333_2 before a_333_10: tasks in index order, as in the concatenated files
        return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", os.path.basename(p))]
    for tsv in sorted(glob.glob(os.path.join(d, tglob)), key=natural):
        tag = os.path.basename(tsv)[2:-4]
        lines = [l for l in open(tsv) if l.strip()]
        if len(lines) != 1:
            problems.append("%s: %d lines, expected 1" % (tsv, len(lines)))
            continue
        al, n, c = parse_alpha_line(lines[0], tsv)
        meta = {"tag": tag}
        log = next((p for p in (os.path.join(d, "a_%s.log" % tag), os.path.join(d, "a_%s.err" % tag))
                    if os.path.exists(p)), None)
        if log is None:
            problems.append("%s: no census log" % tag)
        else:
            summ = [l.rstrip("\n") for l in open(log)]
            other = []
            for l in summ:
                m = SUMMARY.match(l)
                if m:
                    if [int(x) for x in m.groups()[2:13]] != c or int(m.group(2)) != n:
                        problems.append("%s: census summary disagrees with the ALPHA line" % tag)
                    meta["meanpin"] = float(m.group(14))
                elif l.startswith("task "):
                    t = l.split()
                    kv = dict(zip(t[2::2], t[3::2]))
                    meta.update(index=int(t[1]), seed=int(kv["seed"]), start=int(kv["start"]), host=kv.get("host"))
                    if kv["alpha"] != al or int(kv["n"]) != n:
                        problems.append("%s: task header %r disagrees with the ALPHA line" % (tag, l))
                elif l.startswith("end "):
                    meta["end"] = int(l.split()[1])
                elif l.startswith("psweep exact alpha="):
                    meta["redraws"] = int(l.split()[-2])
                else:
                    other.append(l)
            if "meanpin" not in meta:
                problems.append("%s: no census summary line in %s" % (tag, log))
            if other:
                problems.append("%s: unexpected lines in %s: %r" % (tag, log, other[:3]))
            if seedbase is not None and "seed" in meta and meta["seed"] != seedbase + meta["index"]:
                problems.append("%s: seed %d, expected %d + %d" % (tag, meta["seed"], seedbase, meta["index"]))
        pt = os.path.join(d, "pt_%s.txt" % tag)
        if os.path.exists(pt):
            nw = sum(1 for l in open(pt) if l.startswith("W"))
            if nw != c[COLS.index("FR_W")]:
                problems.append("%s: %d W points in %s, FR_W = %d" % (tag, nw, pt, c[COLS.index("FR_W")]))
        else:
            problems.append("%s: no W-point file" % tag)
        rows.append((tag, al, n, c, meta))
    errs = sorted(sum((glob.glob(os.path.join(d, g)) for g in eglob), []))
    nonempty = [e for e in errs if os.path.getsize(e) > 0]
    for e in nonempty:
        problems.append("Slurm stderr not empty: %s" % e)
    seeds = [r[4]["seed"] for r in rows if "seed" in r[4]]
    if unique_seeds and len(seeds) != len(set(seeds)):
        problems.append("repeated seeds among the tasks")
    return rows, len(errs)


def read_source(src, problems, seedbase=None, tglob="a_*.tsv", eglob=("j_*.err", "probe_*.err")):
    if os.path.isdir(src):
        rows, nerr = read_dir(src, problems, seedbase, tglob, eglob)
        return rows, {"kind": "dir", "slurm_err_files": nerr}
    rows = []
    for k, l in enumerate(open(src)):
        if l.startswith("ALPHA"):
            al, n, c = parse_alpha_line(l, "%s:%d" % (src, k + 1))
            rows.append(("line%d" % (k + 1), al, n, c, {}))
    return rows, {"kind": "file"}


# ------------------------------------------------------------------ statistics
def cp(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def analyse(rows, boot):
    T = np.array([r[3] for r in rows], dtype=np.int64)          # tasks x 11
    ntask = np.array([r[2] for r in rows], dtype=np.int64)
    N = int(ntask.sum())
    tot = {k: int(v) for k, v in zip(COLS, T.sum(0))}
    out = {"tasks": len(rows), "n": N, "counts": tot, "rates": {}, "wedges": {}, "dispersion": {}}
    for k in COLS:
        lo, hi = cp(tot[k], N)
        out["rates"][k] = {"p": tot[k] / N, "cp": [lo, hi], "relhalf": (hi - lo) / 2 / (tot[k] / N) if tot[k] else None}
    if boot and len(rows) > 1:
        idx = np.random.default_rng(1).integers(0, len(rows), (boot, len(rows)))
        NB = ntask[idx].sum(1).astype(float)
    for w in ORDER:
        f, e, c = "FR_" + w, MARG[w][0], MARG[w][1]
        p11, p1, p2 = tot[f] / N, tot[e] / N, tot[c] / N
        R = p11 / (p1 * p2)
        se = math.sqrt((1 / p11 - 1 / p1 - 1 / p2 + 2 * p11 / (p1 * p2) - 1) / N) if tot[f] else float("nan")
        d = {"count": tot[f], "rate": p11, "rate_cp": out["rates"][f]["cp"], "indep": p1 * p2, "ratio": R,
             "ratio_delta": [R * math.exp(-Z * se), R * math.exp(Z * se)], "ratio_relse": se,
             "expected_if_indep": p1 * p2 * N}
        a, b = out["rates"][f]["cp"]
        d["ratio_cp_fr_only"] = [a / (p1 * p2), b / (p1 * p2)]
        if boot and len(rows) > 1:
            i1, i2, i3 = COLS.index(f), COLS.index(e), COLS.index(c)
            Rb = (T[:, i1][idx].sum(1) / NB) / ((T[:, i2][idx].sum(1) / NB) * (T[:, i3][idx].sum(1) / NB))
            d["ratio_boot"] = [float(x) for x in np.percentile(Rb, [2.5, 97.5])]
        out["wedges"][w] = d
    if len(rows) > 1:
        m = ntask.astype(float)
        for i, k in enumerate(COLS):
            x = T[:, i].astype(float)
            p = x.sum() / m.sum()
            if 0 < p < 1:
                X2 = float(((x - m * p) ** 2 / (m * p * (1 - p))).sum()); df = len(rows) - 1
                out["dispersion"][k] = {"index": X2 / df,
                                        "p_two_sided": float(2 * min(chi2.cdf(X2, df), chi2.sf(X2, df)))}
    dNE = tot["FR_N"] - tot["FR_E"]
    out["mirror"] = {"N_minus_E": dNE, "z_NE": dNE / math.sqrt(tot["FR_N"] + tot["FR_E"]) if tot["FR_N"] + tot["FR_E"] else None,
                     "rivP_minus_rivM": tot["rivP"] - tot["rivM"],
                     "rel_rivP_rivM": abs(tot["rivP"] - tot["rivM"]) / (0.5 * (tot["rivP"] + tot["rivM"]))}
    return out


# ------------------------------------------------------------------ the table, in the format of ms.tex
def sci(x, sig):
    e = int(math.floor(math.log10(abs(x))))
    return x / 10 ** e, e


def fmt_rate(w, p, lo, hi):
    sig = 5 if w == "S" else 3
    m, e = sci(p, sig)
    s = 10 ** e
    f = "%%.%df" % (sig - 1)
    return (f + r"\ [" + f + r",\," + f + r"]\times10^{%d}") % (m, lo / s, hi / s, e)


def fmt_ratio(R, lo, hi):
    if R >= 1:
        return r"%.3f\ [%.3f,\,%.3f]" % (R, lo, hi)
    inv = 1 / R
    if inv >= 100:
        return r"1/%.0f\ [1/%.0f,\,1/%.0f]" % (inv, 1 / lo, 1 / hi)
    return r"1/%.2f\ [1/%.2f,\,1/%.2f]" % (inv, 1 / lo, 1 / hi)


def fmt_count(k):
    return "{:,}".format(k).replace(",", r"\,") if k >= 10000 else str(k)


def table_rows(res):
    rows = {}
    for w in ORDER:
        d = res["wedges"][w]
        m, e = sci(d["indep"], 3)
        if d["count"] == 0:          # nothing found: the exact upper bounds
            ub = d["rate_cp"][1]; mu, eu = sci(ub, 3)
            rate, ratio = r"0\ [0,\,%.2f\times10^{%d}]" % (mu, eu), r"0\ [0,\,%.3g]" % (ub / d["indep"])
        else:
            rate, ratio = fmt_rate(w, d["rate"], *d["rate_cp"]), fmt_ratio(d["ratio"], *d["ratio_delta"])
        rows[w] = r"$%s$ & $%s$ & $%s$ & $%.2f\times10^{%d}$ & $%s$\\" % (w, fmt_count(d["count"]), rate, m, e, ratio)
    return rows


def inv(x):
    return 1 / x if x else float("inf")


def report(res, label):
    r, n = res["rates"], res["n"]
    print("=" * 110)
    print("%s: %d tasks, %s draws" % (label, res["tasks"], "{:,}".format(n)))
    print("-" * 110)
    print("counts   " + "  ".join("%s %d" % (k, res["counts"][k]) for k in COLS))
    print("\nmarginals (Clopper-Pearson 95%)")
    for k in ("effCC", "effALT", "effZ", "rivP", "rivM"):
        print("  %-6s %.7f [%.7f, %.7f]" % (k, r[k]["p"], *r[k]["cp"]))
    print("\n  %-3s %12s %28s %12s %12s %30s %24s %24s" % ("", "count", "rate [CP 95%]", "rel.half", "indep",
          "obs/indep [delta 95%]", "1/ratio [delta]", "task bootstrap"))
    for w in ORDER:
        d = res["wedges"][w]
        lo, hi = d["ratio_delta"]
        bs = d.get("ratio_boot")
        rh = r["FR_" + w]["relhalf"]
        print("  %-3s %12d %10.5e [%.5e, %.5e] %12s %12.5e %10.6g [%.6g, %.6g] %8.2f [%.2f, %.2f] %s"
              % (w, d["count"], d["rate"], *d["rate_cp"], "%.3f%%" % (100 * rh) if rh else "-", d["indep"], d["ratio"],
                 lo, hi, inv(d["ratio"]), inv(hi), inv(lo), "[%.6g, %.6g]" % tuple(bs) if bs else "-"))
    print("\n  ratio rel. s.e. (delta):  " + "  ".join("%s %.4f%%" % (w, 100 * res["wedges"][w]["ratio_relse"]) for w in ORDER))
    print("  CP on the FR count alone, marginals fixed:  " +
          "  ".join("%s [%.6g, %.6g]" % (w, *res["wedges"][w]["ratio_cp_fr_only"]) for w in ORDER))
    print("  expected if independent:  " + "  ".join("%s %.1f" % (w, res["wedges"][w]["expected_if_indep"]) for w in ORDER))
    if res["dispersion"]:
        ds = res["dispersion"]
        print("\nindex of dispersion over tasks (df %d): " % (res["tasks"] - 1) +
              "  ".join("%s %.3f (p %.2f)" % (k, v["index"], v["p_two_sided"]) for k, v in ds.items()))
        print("  range %.2f to %.2f" % (min(v["index"] for v in ds.values()), max(v["index"] for v in ds.values())))
    mi = res["mirror"]
    print("\nmirror: N - E = %d, z = %s;  rivP - rivM = %d (%.3f%% apart)"
          % (mi["N_minus_E"], "%.2f" % mi["z_NE"] if mi["z_NE"] is not None else "-", mi["rivP_minus_rivM"],
             100 * mi["rel_rivP_rivM"]))
    W, S = res["wedges"]["W"], res["wedges"]["S"]
    if not all(res["wedges"][w]["count"] for w in ORDER):      # a short test run: no text numbers
        print("\nSI Table 3 rows (ms.tex format; a wedge with no friendly rival shows the exact upper bounds):")
        for w, t in table_rows(res).items():
            print("  " + t)
        return
    print("one draw in %s is a friendly rival in W; in S one in %.0f" % ("{:,}".format(int(round(1 / W["rate"], -2))), 1 / S["rate"]))
    print("\nSI Table 3 rows (ms.tex format):")
    for w, t in table_rows(res).items():
        print("  " + t)
    print("legend marginals: efficiency below %.4f, above %.4f; competitiveness %.5f (T>S) and %.5f (T<S)"
          % (r["effCC"]["p"], r["effALT"]["p"], r["rivP"]["p"], r["rivM"]["p"]))
    print("main text: %.1e (W), %.1e (N), %.1e (E), %.1e (S); efficiency %.1f%% / %.1f%%; competitiveness %.1f%%;"
          % (W["rate"], res["wedges"]["N"]["rate"], res["wedges"]["E"]["rate"], S["rate"], 100 * r["effCC"]["p"],
             100 * r["effALT"]["p"], 50 * (r["rivP"]["p"] + r["rivM"]["p"])))
    print("           %.0f x less in W (%.0f to %.0f), %.1f x less in N and E (%.1f to %.1f), %.1f x more in S (%.3f to %.3f)"
          % (1 / W["ratio"], 1 / W["ratio_delta"][1], 1 / W["ratio_delta"][0],
             1 / res["wedges"]["N"]["ratio"], 1 / max(res["wedges"]["N"]["ratio_delta"][1], res["wedges"]["E"]["ratio_delta"][1]),
             1 / min(res["wedges"]["N"]["ratio_delta"][0], res["wedges"]["E"]["ratio_delta"][0]),
             S["ratio"], *S["ratio_delta"]))


def compare(res, ref, label, rlabel):
    print("\n" + "=" * 110)
    print("%s  /  %s   (rate ratio, 95%% interval for two independent binomials on the log scale)" % (label, rlabel))
    out = {}
    for k in COLS + ["ratio_" + w for w in ORDER]:
        if k.startswith("ratio_"):
            w = k[6:]
            a, b = res["wedges"][w], ref["wedges"][w]
            if not (a["count"] and b["count"]):
                continue
            q = a["ratio"] / b["ratio"]
            se = math.hypot(a["ratio_relse"], b["ratio_relse"])
        else:
            k1, n1, k2, n2 = res["counts"][k], res["n"], ref["counts"][k], ref["n"]
            if not (k1 and k2):
                continue
            q = (k1 / n1) / (k2 / n2)
            se = math.sqrt(1 / k1 - 1 / n1 + 1 / k2 - 1 / n2)
        lo, hi = q * math.exp(-Z * se), q * math.exp(Z * se)
        pred = PRED_SHIFT.get(k[3:]) if k.startswith("FR_") else None
        out[k] = {"ratio": q, "ci": [lo, hi], "census_predicted": pred}
        print("  %-9s %+8.3f%%  [%+.3f%%, %+.3f%%]%s" % (k, 100 * (q - 1), 100 * (lo - 1), 100 * (hi - 1),
              "   census-predicted %+.2f%%" % (100 * pred) if pred is not None else ""))
    return out


# ------------------------------------------------------------------ the published alpha = 0.333 numbers
PUBLISHED = {   # PRFull/Text/ms.tex SI Table 3 (rows), its legend, the main text and Methods
    "W": r"$W$ & $6388$ & $3.91\ [3.82,\,4.01]\times10^{-6}$ & $2.63\times10^{-3}$ & $1/673\ [1/690,\,1/657]$\\",
    "S": r"$S$ & $21\,928\,869$ & $1.3426\ [1.3421,\,1.3432]\times10^{-2}$ & $2.63\times10^{-3}$ & $5.096\ [5.095,\,5.097]$\\",
    "N": r"$N$ & $162\,398$ & $9.94\ [9.89,\,9.99]\times10^{-5}$ & $1.18\times10^{-3}$ & $1/11.90\ [1/11.96,\,1/11.84]$\\",
    "E": r"$E$ & $161\,960$ & $9.92\ [9.87,\,9.96]\times10^{-5}$ & $1.18\times10^{-3}$ & $1/11.94\ [1/11.99,\,1/11.88]$\\",
}


def reproduce(res, boot_exact=True):
    """Every number of SI Table 3 and of the rarity sentences, against the text; returns the number of failures."""
    r, w = res["rates"], res["wedges"]
    checks = [("n", res["n"], 1633280000), ("tasks", res["tasks"], 512)]
    rows = table_rows(res)
    for k in ORDER:
        checks.append(("SI Table 3 row " + k, rows[k], PUBLISHED[k]))
    checks += [
        ("legend: efficiency below", "%.4f" % r["effCC"]["p"], "0.1442"),
        ("legend: efficiency above", "%.4f" % r["effALT"]["p"], "0.0648"),
        ("legend: competitiveness T>S", "%.5f" % r["rivP"]["p"], "0.01827"),
        ("legend: competitiveness T<S", "%.5f" % r["rivM"]["p"], "0.01827"),
        ("main: W rate", "%.1e" % w["W"]["rate"], "3.9e-06"), ("main: N rate", "%.1e" % w["N"]["rate"], "9.9e-05"),
        ("main: E rate", "%.1e" % w["E"]["rate"], "9.9e-05"), ("main: S rate", "%.1e" % w["S"]["rate"], "1.3e-02"),
        ("main: W 673x less", "%.0f" % (1 / w["W"]["ratio"]), "673"),
        ("main: W interval 657-690", "%.0f-%.0f" % (1 / w["W"]["ratio_delta"][1], 1 / w["W"]["ratio_delta"][0]), "657-690"),
        ("main: N, E 12x less", "%.0f %.0f" % (1 / w["N"]["ratio"], 1 / w["E"]["ratio"]), "12 12"),
        ("main: N, E interval 11.8-12.0", "%.1f-%.1f" % (1 / max(w["N"]["ratio_delta"][1], w["E"]["ratio_delta"][1]),
                                                    1 / min(w["N"]["ratio_delta"][0], w["E"]["ratio_delta"][0])), "11.8-12.0"),
        ("main: S 5.1x more", "%.1f" % w["S"]["ratio"], "5.1"),
        ("main: S interval 5.095-5.097", "%.3f-%.3f" % tuple(w["S"]["ratio_delta"]), "5.095-5.097"),
        ("Methods: N-E in sd", "%.2f" % res["mirror"]["z_NE"], "0.77"),
        ("Methods: dispersion range", "%.2f-%.2f" % (min(v["index"] for v in res["dispersion"].values()),
                                                      max(v["index"] for v in res["dispersion"].values())), "0.95-1.12"),
        ("old Table 4: rivP vs rivM apart", "%.3f%%" % (100 * res["mirror"]["rel_rivP_rivM"]), "0.023%"),
        ("old Table 4: one in 255,700 (W)", "{:,}".format(int(round(1 / w["W"]["rate"], -2))), "255,700"),
        ("old Table 4: one in 74 (S)", "%.0f" % (1 / w["S"]["rate"]), "74"),
    ]
    # the Ref. 1 triage (report 6, scratchpad rarity_ci.py), as printed there: each number is compared at the
    # precision it was printed with (agreement to half a unit of its last digit)
    tri = {"W": ("0.00148489", "0.00144893", "0.00152175"), "S": ("5.09621", "5.09505", "5.09737"),
           "N": ("0.084023", "0.083617", "0.084431"), "E": ("0.083777", "0.083372", "0.084185")}
    # (bootstrap: as rarity_ci.py prints it; the triage text's S "[5.0950, 5.0975]" is a hand rounding of its 5.09745)
    tboot = {"W": ("0.00145", "0.00152"), "S": ("5.09496", "5.09745"), "N": ("0.08359", "0.08445"), "E": ("0.08338", "0.08417")}
    tcp = {"FR_W": ("3.9111e-6", "3.8158e-6", "4.0083e-6"), "FR_S": ("1.342628e-2", "1.342070e-2", "1.343186e-2"),
           "FR_N": ("9.9431e-5", "9.8948e-5", "9.9915e-5"), "FR_E": ("9.9162e-5", "9.8680e-5", "9.9647e-5"),
           "effCC": ("0.1441852", "0.1441682", "0.1442023"), "effALT": ("0.0647789", "0.0647670", "0.0647909"),
           "rivP": ("0.0182679", "0.0182614", "0.0182744"), "rivM": ("0.0182721", "0.0182656", "0.0182786")}

    def at(vals, pubs, rel=None):   # got, printed at the published precision; passes if within half a unit
        from decimal import Decimal     # (or, for a Monte Carlo number, within the relative tolerance rel)
        out, ok = [], True
        for v, p in zip(vals, pubs):
            q = Decimal(p); unit = Decimal(1).scaleb(q.as_tuple().exponent)
            out.append(str(Decimal(repr(v)).quantize(unit)) if abs(q) >= 1e-3 else "%.*e" % (len(q.as_tuple().digits) - 1, v))
            ok &= (abs(Decimal(repr(v)) - q) <= unit / 2 * Decimal("1.000001") if rel is None
                   else abs(v - float(q)) <= rel * abs(float(q)))
        return " ".join(out), ok
    for k, pub in tri.items():
        g, ok = at([w[k]["ratio"], *w[k]["ratio_delta"]], pub)
        checks.append(("triage: ratio [delta] " + k, g, g if ok else " ".join(pub)))
    # the bootstrap resamples tasks by position, so it is exact only in the file's task order (sweep_a333.tsv is
    # `cat out/a_*.tsv` in Cannon's en_US glob order: 0, 100..109, 10, 110, ...); in any other order its
    # quantiles agree to Monte Carlo precision, about 1e-4 relative at 20 000 resamples
    for k, pub in tboot.items():
        g, ok = at(w[k]["ratio_boot"], pub, None if boot_exact else 3e-4)
        checks.append(("triage: ratio bootstrap " + k + ("" if boot_exact else " (MC, 3e-4)"), g, g if ok else " ".join(pub)))
    for k, pub in tcp.items():
        g, ok = at([r[k]["p"], *r[k]["cp"]], pub)
        checks.append(("triage: rate [CP] " + k, g, g if ok else " ".join(pub)))
    fails = 0
    print("\nREPRODUCTION of the published alpha = 0.333 numbers")
    for name, got, want in checks:
        ok = got == want
        fails += not ok
        print("  %-4s %-34s %s%s" % ("ok" if ok else "FAIL", name, got, "" if ok else "   (published: %s)" % want))
    print("  %d of %d agree" % (len(checks) - fails, len(checks)))
    return fails


# ------------------------------------------------------------------ the probe
OLD_MS_PER_DRAW = 1935.2 * 3600e3 / 1633280000    # 45831015: 1935.2 core-hours (sacct) for 1.633e9 draws, 4.27 ms
OLD_MAX_OVER_MEAN = (3 + 55 / 60) / (1935.2 / 512)  # its longest task, 3:55, over its mean task


def probe(d, jid):
    problems = []
    rows, nerr = read_dir(d, problems, None, unique_seeds=False)   # tasks 0 and 1 share seed 90004 on purpose
    if jid is not None and len(glob.glob(os.path.join(d, "probe_%s_*.err" % jid))) != 4:
        problems.append("expected 4 Slurm stderr files of array %s in %s" % (jid, d))
    by = {r[0]: r for r in rows}
    ref = open(os.path.join(HERE, "ref", "probe45830542_task4.tsv")).read().rstrip("\n").split("\t")
    fails = 0
    print("PROBE %s" % d)
    for tag, want_alpha in (("p333_s90004", "333"), ("p333o1000_s90004", "333/1000")):
        if tag not in by:
            problems.append("no output for %s" % tag); continue
        _, al, n, c, meta = by[tag]
        got = [str(n)] + [str(x) for x in c]
        same = al == want_alpha and got == ref[2:]
        fails += not same
        print("  %-18s alpha %-8s %s the archived probe 45830542 task 4 (seed 90004, 1e5 draws)%s"
              % (tag, al, "IDENTICAL to" if same else "DIFFERS from", "" if same else ": %s vs %s" % (got, ref[2:])))
    rates = []
    for tag, (_, al, n, c, meta) in sorted(by.items()):
        if "end" in meta:
            ms = 1000.0 * (meta["end"] - meta["start"]) / n
            print("  %-18s alpha %-8s n %6d  %5d s  %.2f ms per draw  (host %s, seed %d)"
                  % (tag, al, n, meta["end"] - meta["start"], ms, meta.get("host"), meta["seed"]))
            if al == "1/3":
                rates.append(ms)
    if rates:
        ms = sum(rates) / len(rates)
        mean_h = 3190000 * ms / 3.6e6
        print("  alpha = 1/3: %.2f ms per draw (45831015 at 0.333: %.2f ms) -> a 3 190 000-draw task %.2f h on average,"
              " about %.2f h for the longest (x %.3f as in 45831015); %.0f core-hours in all.  --time=8:00:00 is %s."
              % (ms, OLD_MS_PER_DRAW, mean_h, mean_h * OLD_MAX_OVER_MEAN, OLD_MAX_OVER_MEAN, 512 * mean_h,
                 "ample (x %.1f)" % (8 / (mean_h * OLD_MAX_OVER_MEAN)) if mean_h * OLD_MAX_OVER_MEAN * 1.5 < 8
                 else "NOT enough: raise it (and use intermediate beyond 3 days)"))
    for p in problems:
        print("  PROBLEM: " + p)
    return 1 if (fails or problems) else 0


# ------------------------------------------------------------------ cluster side
def sacct_clean(jobfile):
    jid = open(jobfile).read().strip()
    if not jid.isdigit():
        die("%s does not hold a job ID" % jobfile)
    run = lambda c: subprocess.run(["ssh", "-o", "BatchMode=yes", "cannon", c], capture_output=True, text=True, check=True).stdout
    q = run("squeue -j %s -h -o '%%i %%T'" % jid).strip()
    if q:
        die("array %s is still queued or running:\n%s" % (jid, q))
    s = run("sacct -j %s --format=JobID,State,ExitCode -X --noheader" % jid).splitlines()
    bad = [l for l in s if l.split() and l.split()[1] != "COMPLETED"]
    if not s or bad:
        die("array %s is NOT clean per sacct -- refusing to report:\n%s" % (jid, "\n".join(bad) or "(no sacct rows)"))
    print("sacct: array %s, %d tasks, all COMPLETED" % (jid, len(s)))
    return jid, len(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("source", nargs="*")
    ap.add_argument("--alpha"); ap.add_argument("--tasks", type=int, default=512)
    ap.add_argument("--n", type=int, default=3190000); ap.add_argument("--seedbase", type=int)
    ap.add_argument("--jobid"); ap.add_argument("--fetch"); ap.add_argument("--boot", type=int, default=20000)
    ap.add_argument("--compare"); ap.add_argument("--json"); ap.add_argument("--reproduce", action="store_true")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--glob", default="a_*.tsv"); ap.add_argument("--errglob", action="append")
    a = ap.parse_args()
    if a.probe:
        a.tasks = 4 if a.tasks == 512 else a.tasks
    if a.reproduce and not a.source:
        a.source = [A333]
        a.alpha = a.alpha or "333"
    if not a.source:
        ap.error("no SOURCE")
    if a.jobid:
        jid, ns = sacct_clean(a.jobid)
        if ns != a.tasks:
            die("sacct lists %d tasks, expected %d" % (ns, a.tasks))
        if a.fetch:
            os.makedirs(a.source[0], exist_ok=True)
            subprocess.run(["rsync", "-a", "cannon:%s/" % a.fetch.rstrip("/"), a.source[0].rstrip("/") + "/"], check=True)
    if a.probe:
        sys.exit(probe(a.source[0], jid if a.jobid else None))
    problems, rows, meta = [], [], {}
    eglob = tuple(a.errglob) if a.errglob else (("j_*_%s_*.err" % jid, "probe_%s_*.err" % jid) if a.jobid
                                                 else ("j_*.err", "probe_*.err"))
    for s in a.source:
        r, m = read_source(s, problems, a.seedbase, a.glob, eglob)
        rows += r; meta[s] = m
    if a.jobid:        # only this array's Slurm stderr counts
        for s in a.source:
            if meta[s]["kind"] == "dir":
                n_err = len(glob.glob(os.path.join(s, "j_*_%s_*.err" % jid)))
                if n_err != a.tasks:
                    problems.append("%d Slurm stderr files of array %s in %s, expected %d" % (n_err, jid, s, a.tasks))
    alphas = sorted({r[1] for r in rows})
    if len(alphas) != 1:
        problems.append("more than one alpha among the lines: %r" % alphas)
    if a.alpha and alphas != [a.alpha]:
        problems.append("alpha %r, expected %r" % (alphas, a.alpha))
    if len(rows) != a.tasks:
        problems.append("%d task lines, expected %d" % (len(rows), a.tasks))
    badn = [r[0] for r in rows if r[2] != a.n]
    if badn:
        problems.append("%d tasks with n != %d: %s" % (len(badn), a.n, badn[:5]))
    idx = sorted(r[4]["index"] for r in rows if "index" in r[4])
    if idx and idx != list(range(len(rows))):
        problems.append("task indices are not 0..%d" % (len(rows) - 1))
    if problems:
        print("NOT CLEAN -- refusing to report:")
        for p in problems:
            print("  " + p)
        sys.exit(1)
    res = analyse(rows, a.boot)
    res["alpha"] = alphas[0]
    label = "alpha = %s (%s)" % (alphas[0], ", ".join(a.source))
    print("checks passed: %d tasks, alpha %s, n %d each%s" % (len(rows), alphas[0], a.n,
          "; logs, W-point files and Slurm stderr consistent" if any(m["kind"] == "dir" for m in meta.values()) else ""))
    times = [r[4]["end"] - r[4]["start"] for r in rows if "end" in r[4] and "start" in r[4]]
    if times:
        res["walltime_s"] = {"min": min(times), "max": max(times), "sum": sum(times)}
        print("wall time per task %.2f to %.2f h, %.1f core-hours, %.2f ms per draw"
              % (min(times) / 3600, max(times) / 3600, sum(times) / 3600, 1000 * sum(times) / res["n"]))
    red = [r[4]["redraws"] for r in rows if "redraws" in r[4]]
    if red:
        res["redraws"] = sum(red)
        print("exact-rational redraws: %d in all" % sum(red))
    report(res, label)
    fails = 0
    if a.reproduce:
        fails = reproduce(res, boot_exact=all(m["kind"] == "file" for m in meta.values()))
    if a.compare:
        p2 = []
        r2, _ = read_source(a.compare, p2)
        if p2:
            die("comparison source not clean: %r" % p2)
        ref = analyse(r2, 0)
        res["compare"] = {"source": a.compare, "alpha": sorted({r[1] for r in r2}),
                          "ratios": compare(res, ref, "alpha " + alphas[0], "alpha " + ",".join(sorted({r[1] for r in r2})))}
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
