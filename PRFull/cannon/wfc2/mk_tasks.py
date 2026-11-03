#!/usr/bin/env python3
"""mk_tasks.py -- writes every task file, game list and Slurm script of the wfc2 kit, so they cannot disagree.

    python3 mk_tasks.py          (local, instant; re-run after changing a list, a seed or a --time below)

Game numbers are TASK numbers = lines of games.dat (1-based); the sunflower index is ipt = task - 1 (REP field 3).
Seeds: replicate r of game g is seeded seedbase + 10 (g - 1) + r (production's rule).  Production used the bases
9000000 (F1), 9100000 (F2), 9500000 and 9600000 (probes), neutral draws 777000000 + g and 8800134xx; every base below
is at least 3.1e7, and the arms' ranges (base + 1 .. base + 5120) are pairwise disjoint.  The probe's replicates use
base + PROBEOFF (250000): no arm's base is 250000 off another's, so every probe seed is disjoint from every arm's and
from every other probe arm's (main() asserts it); only DEF, DEF2 and TSV use production's seeds, on purpose.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TAGS = [l.split()[7] for l in open(os.path.join(HERE, "games.dat")) if l.strip() and not l.startswith("#")]
assert len(TAGS) == 512 and all(TAGS[i].startswith("V_p%03d_" % i) for i in range(512))
WEDGE = {t + 1: TAGS[t][-1] for t in range(512)}           # by task number
QUAD = {t + 1: TAGS[t][-3:-1] for t in range(512)}

# ---- the triage's 64 games (report 5): 16 per wedge, evenly spaced in sunflower index, with the 8 validation games
G64 = {"W": [1, 13, 30, 38, 56, 82, 114, 148, 187, 224, 263, 302, 342, 378, 420, 465],
       "S": [18, 47, 83, 115, 146, 172, 180, 209, 243, 277, 311, 345, 379, 413, 447, 481],
       "N": [24, 58, 92, 126, 160, 184, 192, 223, 257, 291, 325, 354, 388, 417, 448, 480],
       "E": [55, 110, 144, 165, 178, 212, 238, 262, 267, 296, 322, 351, 382, 406, 432, 461]}
VALID8 = [30, 38, 47, 58, 165, 172, 184, 262]
# ---- the special games of the critic (ipt): 466 411 301 (W friendly-rival takeovers), 237 326 415 (HAE) and 500 (PDN)
# where friendly rivals lead outside S, 246 (the one PD without a partner lead) and 337 (the one game above the switch
# line without a partner lead).  301 and 237 are already among the 64 (tasks 302 and 238), so the subset has 71 games.
SPECIAL_IPT = [466, 411, 301, 237, 326, 415, 500, 246, 337]
G71 = sorted(set(sum(G64.values(), [])) | {i + 1 for i in SPECIAL_IPT})
# ---- the 24 time-series games (report 6: 6 per wedge), ipt.  Split = replicates of production F1 disagree on the
# leading atom (35 games); typical = unsplit, of the 64, efficiency nearest its wedge-and-quadrant median.
TS_IPT = {"W": [466, 411, 301, 186, 377, 29],     # 2 split PDW takeovers, 301 (111 up to 6%), 3 typical PDW
          "S": [184, 417, 451, 208, 17, 171],     # 3 split SHS, typical HAS 208 17, typical SHS 171 (validation 172)
          "N": [500, 390, 445, 479, 337, 256],    # the 4 split PDN (500 special), 337 special split SDN, typical SDN
          "E": [237, 326, 415, 109, 164, 117]}    # the 3 special split HAE, 3 split SDE
G24TS = sorted(i + 1 for v in TS_IPT.values() for i in v)
# ---- the 24-game F2 subset for E3F2/E5F2: representative, 6 per wedge of the 64 = the 2 validation games + 4 evenly
# spaced others (positions 0, 4, 9, 13 of the remaining 14)
G24F2 = []
for w in "WSNE":
    val = [g for g in G64[w] if g in VALID8]; rest = [g for g in G64[w] if g not in VALID8]
    G24F2 += val + [rest[i] for i in (0, 4, 9, 13)]
G24F2 = sorted(G24F2)
for w in "WSNE":
    assert all(WEDGE[g] == w for g in G64[w]) and all(WEDGE[i + 1] == w for i in TS_IPT[w])
assert len(G71) == 71 and len(G24TS) == 24 and len(G24F2) == 24 and set(VALID8) <= set(G71)

F1 = dict(itend=10000000, NT=1000, NI=4)
F2 = dict(itend=100000000, NT=4000, NI=1)
TS = dict(itend=40000000, NT=1000, NI=4)      # 4x F1; production-style sampling over [2e7, 4e7] (TS_* below for the rest)
TSENV = "TS_NT=2000 TS_NI=10 TS_NB=400"

# ---- the arms: name -> (binary, alpha, games, regime, seedbase, neutral draws, neutral seed base, --time, what)
ARMS = [
    ("A1", "wfc", "1/10", G71, F1, 31000000, 20000, 38000000, "0-06:00:00", "alpha = 1/10, global kernel, F1"),
    ("A2", "wfc", "1/5", G71, F1, 32000000, 20000, 38100000, "0-04:00:00", "alpha = 1/5, global kernel, F1"),
    ("A3", "wfc", "2/5", G71, F1, 30000000, 20000, 38200000, "0-04:00:00", "alpha = 2/5, global kernel, F1"),
    ("K0", "wfc_k0", "1/3", G71, F1, 33000000, 0, 0, "0-08:00:00", "local kernel K0 at alpha = 1/3, F1"),
    ("E3", "wfc_e3", "1/3", G71, F1, 35000000, 0, 0, "0-04:00:00", "eps = 1e-3, F1"),
    ("E5", "wfc_e5", "1/3", G71, F1, 36000000, 0, 0, "0-04:00:00", "eps = 1e-5, F1"),
    ("E3F2", "wfc2_e3", "1/3", G24F2, F2, 35500000, 0, 0, "0-04:00:00", "eps = 1e-3, F2 (N=100, beta=3, mu=1e-4)"),
    ("E5F2", "wfc2_e5", "1/3", G24F2, F2, 36500000, 0, 0, "0-04:00:00", "eps = 1e-5, F2 (N=100, beta=3, mu=1e-4)"),
]
B0 = [("B0", "wfc_b0", 34000000), ("B0K0", "wfc_b0k0", 34500000)]        # beta = 0, the 8 validation games
B0TIME = "0-08:00:00"
TSSEED, TSTIME = 37000000, "0-02:00:00"
PROBEOFF = 250000        # probe seeds: base + PROBEOFF (with +500000 the E3/E5/E3F2/E5F2 probes reused E3F2/E5F2/E5/TS seeds)

HDR = "# outdir bin alpha game rep0 nrep itend NT NI seedbase neutralM neutralseed -- "


def line(od, b, a, g, r0, nr, reg, sb, nm, ns):
    return "%s %s %s %d %d %d %d %d %d %d %d %d" % (od, b, a, g, r0, nr, reg["itend"], reg["NT"], reg["NI"], sb, nm, ns)


def sbatch(name, ntask, time, taskfile, part="sapphire", env=""):
    return """#!/bin/bash
#SBATCH --job-name=PRFull.wfc2.{name}
#SBATCH --comment=PRFull
#SBATCH --partition={part}
#SBATCH --array=1-{ntask}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time={time}
#SBATCH --output=out/{name}_%A_%a.out
#SBATCH --error=out/{name}_%A_%a.err
export TASKFILE={taskfile}{env}
srun ./task.sh
""".format(name=name, ntask=ntask, time=time, taskfile=taskfile, part=part, env=(" " + env if env else ""))


def write(fn, txt):
    open(os.path.join(HERE, fn), "w").write(txt)


def main():
    summary = []
    for name, b, a, games, reg, sb, nm, nsb, time, what in ARMS:
        L = [HDR + "%s: %d games x 10 replicates" % (what, len(games))]
        L += [line(name, b, a, g, 1, 10, reg, sb, nm, (nsb + g) if nm else 0) for g in games]
        write("tasks_%s.txt" % name, "\n".join(L) + "\n")
        write("run_%s.sbatch" % name, sbatch(name, len(games), time, "tasks_%s.txt" % name))
        summary.append((name, len(games), len(games) * 10, reg["itend"]))
    L = [HDR + "beta = 0 (neutral), global kernel (B0) and local kernel K0 (B0K0): 8 validation games x 10 replicates"]
    for od, b, sb in B0:
        L += [line(od, b, "1/3", g, 1, 10, F1, sb, 0, 0) for g in VALID8]
    write("tasks_B0.txt", "\n".join(L) + "\n")
    write("run_B0.sbatch", sbatch("B0", 2 * len(VALID8), B0TIME, "tasks_B0.txt"))
    summary.append(("B0", 16, 160, F1["itend"]))
    L = [HDR + "time series, 4x length, 24 games x 10 replicates, one replicate per task (%s)" % TSENV]
    L += [line("TS", "wfc_ts", "1/3", g, r, 1, TS, TSSEED, 0, 0) for g in G24TS for r in range(1, 11)]
    write("tasks_TS.txt", "\n".join(L) + "\n")
    write("run_TS.sbatch", sbatch("TS", 240, TSTIME, "tasks_TS.txt", env=TSENV))
    summary.append(("TS", 240, 240, TS["itend"]))

    # ---- the probe: partition test, <= 8 tasks per arm.  DEF/DEF2/TSV repeat archived production replicates exactly
    short = {k: dict(v, itend=v["itend"] // 10) for k, v in (("F1", F1), ("F2", F2), ("TS", TS))}
    P = [HDR + "PROBE (partition test).  DEF, DEF2, TSV: production settings and seeds, must equal the archived REP/NEU "
               "lines; the rest: 1/10-length replicates (cost, sanity) and one or two full replicates (per-replicate rate)"]
    P += [line("probe/DEF", "wfc", "1/3", 30, 1, 1, F1, 9000000, 20000, 777000030),       # = a333 r_1o3_30_1, n_1o3_30
          line("probe/DEF2", "wfc2", "1/3", 511, 1, 1, F2, 9100000, 0, 0),                 # = f2 r_1o3_511_1
          ] + [line("probe/TSV", "wfc_ts", "1/3", g, 1, 1, F1, 9000000, 0, 0) for g in VALID8]  # = a333 r_1o3_<g>_1, TS on
    four = [30, 47, 58, 165]
    for name, b, a, games, reg, sb, nm, nsb, time, what in ARMS:
        rs = short["F2"] if reg is F2 else short["F1"]
        # <arm>_s: 1/10-length replicates, <arm>_f: full length (separate directories: the file names carry only
        # alpha, game and rep0, so a short and a full run at one game must not share a directory)
        P += [line("probe/%s_s" % name, b, a, g, 1, 1, rs, sb + PROBEOFF, 0, 0) for g in four]
        full = [30] if reg is F2 else ([467, 47] if name == "K0" else [467])
        P += [line("probe/%s_f" % name, b, a, g, 1, 1, reg, sb + PROBEOFF, nm, (nsb + 500000 + g) if nm else 0) for g in full]
    for od, b, sb in B0:
        P += [line("probe/%s_s" % od, b, "1/3", g, 1, 1, short["F1"], sb + PROBEOFF, 0, 0) for g in (30, 47)]
        P += [line("probe/%s_f" % od, b, "1/3", 30, 1, 1, F1, sb + PROBEOFF, 0, 0)]
    P += [line("probe/TS_s", "wfc_ts", "1/3", g, 1, 1, short["TS"], TSSEED + PROBEOFF, 0, 0) for g in (467, 185)]
    P += [line("probe/TS_f", "wfc_ts", "1/3", 467, 1, 1, TS, TSSEED + PROBEOFF, 0, 0)]
    outs = [l.split()[0] + " " + l.split()[3] + " " + l.split()[4] + " " + l.split()[2] for l in P[1:]]
    assert len(set(outs)) == len(outs), "two probe tasks would write the same files"
    # replicate seeds: disjoint between arms (probe arms counted as their own arm), and none in production's ranges
    owner = {}
    for fn in ["tasks_%s.txt" % a[0] for a in ARMS] + ["tasks_B0.txt", "tasks_TS.txt"]:
        for l in open(os.path.join(HERE, fn)):
            if not l.startswith("#"):
                f = l.split(); g, r0, nr, sb = int(f[3]), int(f[4]), int(f[5]), int(f[9])
                for r in range(r0, r0 + nr):
                    owner.setdefault(sb + 10 * (g - 1) + r, set()).add(f[0])
    for l in P[1:]:
        f = l.split(); g, r0, nr, sb = int(f[3]), int(f[4]), int(f[5]), int(f[9])
        arm = f[0].split("/")[1]
        if arm in ("DEF", "DEF2", "TSV"):
            assert sb in (9000000, 9100000), "the validation tasks must use production's seeds"
            continue
        for r in range(r0, r0 + nr):
            owner.setdefault(sb + 10 * (g - 1) + r, set()).add("probe_" + arm.rsplit("_", 1)[0])
    shared = {k: v for k, v in owner.items() if len(v) > 1}
    assert not shared, "seeds shared between arms: %s" % sorted(shared.items())[:5]
    assert min(owner) > 9600000 + 5120, "a seed falls in production's range"
    write("tasks_probe.txt", "\n".join(P) + "\n")
    write("probe.sbatch", sbatch("probe", len(P) - 1, "0-03:00:00", "tasks_probe.txt", part="test", env=TSENV))

    # ---- game lists (for the M1 kit and the analysis)
    G = ["# task ipt tag wedge in71 inTS24 inF2sub valid8"]
    for t in range(1, 513):
        if t in G71 or t in G24TS or t in G24F2:
            G.append("%d %d %s %s %d %d %d %d" % (t, t - 1, TAGS[t - 1], WEDGE[t], t in G71, t in G24TS, t in G24F2, t in VALID8))
    write("games_subsets.txt", "\n".join(G) + "\n")
    print("arm     tasks  replicates  generations")
    for s in summary:
        print("%-6s %6d %11d %12d" % s)
    print("probe  %6d" % (len(P) - 1))
    print("71 games:", " ".join(map(str, G71)))
    print("TS 24 (task):", " ".join(map(str, G24TS)))
    print("F2 24 (task):", " ".join(map(str, G24F2)))


if __name__ == "__main__":
    main()
