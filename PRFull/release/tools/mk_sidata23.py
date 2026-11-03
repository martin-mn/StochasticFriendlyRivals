#!/usr/bin/env python3
"""mk_sidata23.py -- SI Data 2 (the sensitivity runs) and SI Data 3 (the reruns at four times the run length), from
the kits' per-game tables, as deposited in PRFull/Text/SIData2_sensitivity.csv and SIData3_runlength.csv.

    python3 mk_sidata23.py [--prf PRFULL] [--wfc2 DIR] [--cw2 DIR] [--sidata1 CSV] [--wfc2-summary FILE] [--out DIR]
                                                                                    (seconds; standard library only)
Defaults: the deposited tables under PRFULL (this tree's PRFull): cannon/wfc2, cannon/cw2, Text/SIData1_games512.csv,
cannon/wfc2/analysis/summary.txt.  --out writes the two files there; the verdict compares them with the deposited
PRFull/Text files byte for byte.

Written after the fact (2026-09-27, referee 2, major 6): the two CSV files had been assembled from the tables below
without a script; this one rebuilds them, and reproduce_figures.sh (step sidata23) checks that it gives the deposited
files byte for byte.
  SI Data 2, one row per (arm, game): memory two, arms A1 A2 A3 K0 E3 E5 (N = 1000) and E3F2 E5F2 (N = 100), from
    cannon/wfc2/analysis/<arm>_games.tsv (analyze.py), with the production run's most enriched atom at the game
    from SI Data 1; memory one, arms M1E3 M1E5, from cannon/cw2/results/E3_games.tsv, E5_games.tsv (collect.py).
  SI Data 3, one row per game of the 24: memory one from cannon/cw2/results/TS_games.tsv; memory two from the
    time-series block that analyze.py TS prints (cannon/wfc2/analysis/summary.txt, the last "analyze.py ... TS" run,
    or --wfc2-summary: the output of a fresh `analyze.py ... TS`), since analysis/TS_games.tsv keeps only the arm
    table.
"""
import argparse
import csv
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PRF0 = os.path.dirname(os.path.dirname(HERE))
ORDER = ["000", "100", "010", "001", "110", "011", "111"]
M2ARMS = [("A1", 1000, "alpha=1/10"), ("A2", 1000, "alpha=1/5"), ("A3", 1000, "alpha=2/5"), ("K0", 1000, "local kernel"),
          ("E3", 1000, "eps=1e-3"), ("E5", 1000, "eps=1e-5"), ("E3F2", 100, "eps=1e-3"), ("E5F2", 100, "eps=1e-5")]
M1ARMS = [("M1E3", "E3", "eps=1e-3"), ("M1E5", "E5", "eps=1e-5")]
H2 = ["run", "memory", "N", "change", "game", "tag", "wedge", "efficiency", "most_abundant", "replicates_naming_it",
      "most_enriched", "enriched_determined"] + ["share_" + a for a in ORDER] + ["share_competitive", "ref_efficiency",
      "ref_most_abundant", "ref_most_enriched"]
H3 = ["memory", "game", "tag", "wedge", "production_replicates_agree", "production_leader", "leader_window_P",
      "leader_window_L", "eff_P", "eff_L", "eff_production", "z_quarters_eff", "z_quarters_leader_share", "leader_share_P",
      "leader_share_L", "leader_share_production", "ESS_min", "ESS_median", "takeovers_111", "takeover_rate_per_1e7",
      "escapes_111", "escape_rate_per_1e7"]


def tsv(path):
    lines = [l.rstrip("\n") for l in open(path) if l.strip() and not l.startswith("#")]
    hdr = lines[0].split("\t")
    out = []
    for l in lines[1:]:
        f = l.split("\t")
        out.append(dict(zip(hdr, f)) if len(set(hdr)) == len(hdr) else f)
    return hdr, out


def sidata2(wfc2, cw2, sid1):
    rows = []
    for arm, n, change in M2ARMS:
        _, R = tsv(os.path.join(wfc2, "analysis", arm + "_games.tsv"))
        ref = "m2_N%d_most_enriched" % n
        for r in R:
            rows.append([arm, "2", str(n), change, r["ipt"], r["tag"], r["wedge"], r["eff"], r["lead"], r["agree"],
                         r["most_enriched"], r["determined"]] + [r["s" + a] for a in ORDER] +
                        [r["comp_share"], r["prod_eff"], r["prod_lead"], sid1[int(r["ipt"])][ref]])
    for run, arm, change in M1ARMS:
        _, R = tsv(os.path.join(cw2, "results", arm + "_games.tsv"))
        for r in R:
            rows.append([run, "1", "1000", change, r["ipt"], r["tag"], r["wedge"], r["E_arm"], r["lead_arm"], r["agree_arm"],
                         r["enr_arm"], ""] + [r["a" + a] for a in ORDER] +
                        [r["comp_arm"], r["E_prod"], r["lead_prod"], r["enr_prod"]])
    return rows


TSLINE = re.compile(r"^\s+(\d+)\s+([A-Z]{3})\s+(yes|no)\s+(\d{3})\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+([\d.]+)\s+"
                    r"([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+(\d+)\s+(\d+)\s+up\s+(\S+)\s.*?down\s+(\S+)\s")
LEADS = re.compile(r"^\s+ipt (\d+): window P leads (\d{3}), L leads (\d{3})")


def ts_block(summary):
    """the per-game lines of the last TS block of analyze.py's printed output"""
    txt = open(summary).read()
    k = [m.start() for m in re.finditer(r"^TS: 24 games", txt, re.M)]
    if not k:
        sys.exit("no TS block in %s" % summary)
    blk = txt[k[-1]:]
    blk = blk[:blk.index("CRITERION")]
    games, leads = [], {}
    for l in blk.split("\n"):
        m = TSLINE.match(l)
        if m:
            games.append(m.groups())
        m = LEADS.match(l)
        if m:
            leads[m.group(1)] = (m.group(2), m.group(3))
    return games, leads


def sidata3(wfc2, cw2, summary):
    rows = []
    hdr, R = tsv(os.path.join(cw2, "results", "TS_games.tsv"))
    ix = {h: i for i, h in reversed(list(enumerate(hdr)))}          # first occurrence of a repeated name
    rate = [i for i, h in enumerate(hdr) if h == "rate_per1e7"]      # takeovers' rate, then escapes'
    for f in R:
        g = lambda k: f[ix[k]]
        rows.append(["1", g("ipt"), g("tag"), g("wedge"), "yes" if g("prod_reps") == "agree" else "no", g("lead_prod"),
                     g("lead_P"), g("lead_L"), g("E_P"), g("E_L"), g("E_prod"), g("zmax_eff"), g("zmax_lead"), "", "", "",
                     g("ESSlead_P_min"), g("ESSlead_P_med"), g("takeovers"), f[rate[0]], g("escapes"), f[rate[1]]])
    _, T = tsv(os.path.join(wfc2, "analysis", "TS_games.tsv"))
    tag = {r["ipt"]: (r["tag"], r["wedge"]) for r in T}
    games, leads = ts_block(summary)
    for (ipt, q, split, lead, eP, eL, eprd, zq, lP, lL, lpr, zql, essmin, essmed, up, down) in games:
        wP, wL = leads.get(ipt, (lead, lead))
        rows.append(["2", ipt, tag[ipt][0], tag[ipt][1], "no" if split == "yes" else "yes", lead, wP, wL, eP, eL, eprd, zq,
                     zql, lP, lL, lpr, essmin, essmed, "", up, "", down])
    return rows


def write(rows, hdr):
    b = io.StringIO(newline="")
    w = csv.writer(b)                                                # the csv module's own line ends, CR LF
    w.writerow(hdr)
    w.writerows(rows)
    return b.getvalue()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prf", default=PRF0)
    ap.add_argument("--wfc2", default=None, help="the wfc2 kit (default PRFULL/cannon/wfc2)")
    ap.add_argument("--cw2", default=None, help="the cw2 kit (default PRFULL/cannon/cw2)")
    ap.add_argument("--sidata1", default=None, help="SI Data 1 (default PRFULL/Text/SIData1_games512.csv)")
    ap.add_argument("--wfc2-summary", default=None)
    ap.add_argument("--out", default=None, help="directory for the two CSV files (default: print their sizes only)")
    a = ap.parse_args()
    wfc2 = a.wfc2 or os.path.join(a.prf, "cannon", "wfc2")
    cw2 = a.cw2 or os.path.join(a.prf, "cannon", "cw2")
    sid1 = list(csv.DictReader(open(a.sidata1 or os.path.join(a.prf, "Text", "SIData1_games512.csv"))))
    assert [int(r["game"]) for r in sid1] == list(range(512))
    s2 = write(sidata2(wfc2, cw2, sid1), H2)
    s3 = write(sidata3(wfc2, cw2, a.wfc2_summary or os.path.join(wfc2, "analysis", "summary.txt")), H3)
    for name, s in (("SIData2_sensitivity.csv", s2), ("SIData3_runlength.csv", s3)):
        dep = os.path.join(PRF0, "Text", name)
        same = os.path.exists(dep) and open(dep, "rb").read() == s.encode()
        print("%s: %d rows, %s the deposited file" % (name, s.count("\n") - 1, "byte-identical to" if same else "DIFFERS from"))
        if a.out:
            with open(os.path.join(a.out, name), "w", newline="") as fh:
                fh.write(s)


if __name__ == "__main__":
    main()
