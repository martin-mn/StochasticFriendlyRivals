#!/usr/bin/env python3
"""mktasks.py -- write the task files of PRFull cannon/cw2 from games71.txt and games24.txt (deterministic).

    python3 mktasks.py            -> tasks_e3.txt tasks_e5.txt tasks_ts.txt tasks_probe.txt tasks_local.txt
                                     tasks_tsx.txt (the OPTIONAL arm TSX)
    python3 mktasks.py --check    -> exit 1 if any task file on disk differs from what this script writes

Seeds (cw.f: replicate irep of task isl uses seedbase + (isl-1)*10 + irep, MT19937 init_genrand):
  production F1  isl = 3*ipt+3, base 40000000 -> 40000000 + (3*ipt+2)*10 + irep       (PartnersRivalsM1f run1)
  E3             isl = 3*ipt+2, base 40000000 -> 40000000 + (3*ipt+1)*10 + irep       disjoint (other isl)
  E5             isl = 3*ipt+1, base 40000000 -> 40000000 + (3*ipt  )*10 + irep       disjoint (other isl)
  TS             isl = 3*ipt+3, base 43000000 -> 43000000 + (3*ipt+2)*10 + irep       disjoint (other base)
  taken elsewhere: 40M (F1 production, and the E3/E5 slots above), 41M+isl (kappa), 42M (F2), 6M..39M (other
  campaigns, cw/README.md).  All E3/E5/TS seeds lie in [40000001, 40015360] or [43000001, 43015360], and the
  (isl, irep) pairs of E3, E5 and production never coincide, so no two replicates of any campaign share a seed.
The probe rows VAL and VALTS deliberately use the PRODUCTION seed (isl 3, replicate 1): they must reproduce the
archived production line ref/prod/3 line 1 bit for bit.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def games(fn):
    out = []
    for l in open(os.path.join(HERE, fn)):
        if l.startswith("#") or not l.strip():
            continue
        f = l.split()
        out.append((int(f[0]), f[3]))
    return out


G71 = games("games71.txt")
G24 = games("games24.txt")
GTX = games("games_tsx.txt")
assert len(G71) == 71 and len(G24) == 24
assert len({i for i, _ in G71}) == 71 and len({i for i, _ in G24}) == 24

HDR = "# arm  exe  isl  ipt  irep0  irep1  seedbase  eps  itend  outdir     (task.sh; k = data line number)\n"
TSBASE = 43000000


def e3():
    s = ("# PRFull cw2 arm E3: eps = 1e-3 (eps slot ie = 2, NO code change: cw.x = the production program), F1 regime\n"
         "# N = 1000, beta = 100, u = 1e-2, itend = 1e7, 10 replicates per task (the production layout), %d games\n"
         "# (games71.txt).  Seeds 40000000 + (isl-1)*10 + irep with isl = 3*ipt+2: disjoint from production.\n" % len(G71))
    s += HDR
    for i, t in G71:
        s += "E3 cw.x %d %d 1 10 - 1e-3 10000000 out/E3     # %s\n" % (3 * i + 2, i, t)
    return s


def e5():
    s = ("# PRFull cw2 arm E5: eps = 1e-5 (eps slot ie = 1 of the -DEPS1=1.0d-5 build cw_e5.x), F1 regime, 10\n"
         "# replicates per task, %d games (games71.txt).  Seeds 40000000 + (isl-1)*10 + irep, isl = 3*ipt+1: disjoint.\n"
         % len(G71))
    s += HDR
    for i, t in G71:
        s += "E5 cw_e5.x %d %d 1 10 - 1e-5 10000000 out/E5     # %s\n" % (3 * i + 1, i, t)
    return s


def ts():
    s = ("# PRFull cw2 arm TS: time series at 4x production length (itend 4e7, 2000 blocks of 2e4 generations),\n"
         "# eps = 1e-4 (slot 3), F1 regime, %d games (games24.txt) x 10 replicates, ONE replicate per task.\n"
         "# Seeds %d + (isl-1)*10 + irep with isl = 3*ipt+3: disjoint from production (base 40000000).\n"
         "# (Using base 40000000 instead would make every TS run the exact continuation of the published replicate.)\n"
         % (len(G24), TSBASE))
    s += HDR
    for i, t in G24:
        for r in range(1, 11):
            s += "TS cw_ts.x %d %d %d %d %d 1e-4 40000000 out/TS/p%03d_r%02d     # %s\n" % (3 * i + 3, i, r, r, TSBASE, i, r, t)
    return s


def tsx():
    s = ("# PRFull cw2 OPTIONAL arm TSX: the TS build and design at the %d memory-one split games of games_tsx.txt,\n"
         "# 10 replicates, one per task.  Seeds %d + (isl-1)*10 + irep (isl = 3*ipt+3), disjoint from production and\n"
         "# from TS (no TSX game is a TS game).\n" % (len(GTX), TSBASE))
    assert not ({i for i, _ in GTX} & {i for i, _ in G24})
    s += HDR
    for i, t in GTX:
        for r in range(1, 11):
            s += "TSX cw_ts.x %d %d %d %d %d 1e-4 40000000 out/TSX/p%03d_r%02d     # %s\n" % (3 * i + 3, i, r, r, TSBASE, i, r, t)
    return s


def probe():
    rows = [
        ("VAL", "cw.x", 3, 0, 1, 1, "-", "1e-4", 10000000, "out/probe/VAL",
         "default build, PRODUCTION seed: must equal ref/prod/3 line 1 and ref/prod/w3 line 3 bit for bit"),
        ("VALTS", "cw_ts1.x", 3, 0, 1, 1, "-", "1e-4", 10000000, "out/probe/VALTS",
         "-DTSERIES at itend 1e7, production seed: the same two lines; blocks asserted = atom counts"),
        ("PE3", "cw.x", 3 * 466 + 2, 466, 1, 1, "-", "1e-3", 10000000, "out/probe/E3",
         "E3 full replicate at the slowest production game (ipt 466, 3094 s per task)"),
        ("PE3s", "cwp.x", 3 * 0 + 2, 0, 1, 1, "-", "1e-3", 1000000, "out/probe/E3s", "E3 short, PD centre"),
        ("PE3s", "cwp.x", 3 * 480 + 2, 480, 1, 1, "-", "1e-3", 1000000, "out/probe/E3s", "E3 short, deep SHS"),
        ("PE5", "cw_e5.x", 3 * 466 + 1, 466, 1, 1, "-", "1e-5", 10000000, "out/probe/E5", "E5 full replicate, ipt 466"),
        ("PE5s", "cwp_e5.x", 3 * 0 + 1, 0, 1, 1, "-", "1e-5", 1000000, "out/probe/E5s", "E5 short, PD centre"),
        ("PE5s", "cwp_e5.x", 3 * 480 + 1, 480, 1, 1, "-", "1e-5", 1000000, "out/probe/E5s", "E5 short, deep SHS"),
        ("PTS", "cw_ts.x", 3 * 466 + 3, 466, 1, 1, TSBASE, "1e-4", 40000000, "out/probe/TS",
         "TS full task (4e7) at ipt 466 = the TS arm's task (466, rep 1)"),
        ("PTSs", "cwp_ts.x", 3 * 0 + 3, 0, 1, 1, TSBASE, "1e-4", 4000000, "out/probe/TSs", "TS short (4e6), PD centre"),
        ("PTSs", "cwp_ts.x", 3 * 184 + 3, 184, 1, 1, TSBASE, "1e-4", 4000000, "out/probe/TSs", "TS short, split SHS"),
    ]
    s = ("# PRFull cw2 PROBE (partition test): validates each switch on Cannon's compiler and measures the rates.\n"
         "# Rows 1-2 validate the default build and -DTSERIES against the archived production line (isl 3, rep 1);\n"
         "# rows 3-11 give the per-replicate cost of each arm (a full replicate at the slowest game + short runs);\n"
         "# rows 12-13 compute the kernel masses at eps 1e-3 and 1e-5 at the 71 games (M = 1e7, 2 s per game).\n")
    s += HDR
    for r in rows:
        s += "%s %s %d %d %d %d %s %s %d %s     # %s\n" % r
    s += "KA kappa.x 2 games71.txt 10000000 - - 1e-3 - out/probe/KA3     # kernel masses, eps 1e-3 (slot 2)\n"
    s += "KA kappa_e5.x 1 games71.txt 10000000 - - 1e-5 - out/probe/KA5     # kernel masses, eps 1e-5 (slot 1, -DEPS1)\n"
    return s


def local(its=20000):
    s = ("# PRFull cw2 LOCAL smoke test of every arm (test builds at ITS = %d: bash compile.sh test), a few seconds each.\n"
         % its)
    s += HDR
    s += "E3 cw_t.x %d %d 1 10 - 1e-3 %d out/local/E3     # E3 at ipt 466\n" % (3 * 466 + 2, 466, its)
    s += "E5 cw_te5.x %d %d 1 10 - 1e-5 %d out/local/E5     # E5 at ipt 466\n" % (3 * 466 + 1, 466, its)
    for r in (1, 2):
        s += "TS cw_tts.x %d %d %d %d %d 1e-4 %d out/local/TS/p466_r%02d     # TS at ipt 466\n" % (
            3 * 466 + 3, 466, r, r, TSBASE, its, r)
    s += "KA kappa_e5.x 1 games24.txt 100000 - - 1e-5 - out/local/KA5     # kernel masses, 1e5 draws\n"
    return s


FILES = {"tasks_e3.txt": e3(), "tasks_e5.txt": e5(), "tasks_ts.txt": ts(), "tasks_probe.txt": probe(),
         "tasks_local.txt": local(), "tasks_tsx.txt": tsx()}

if __name__ == "__main__":
    bad = 0
    for fn, txt in FILES.items():
        p = os.path.join(HERE, fn)
        if "--check" in sys.argv:
            if not os.path.exists(p) or open(p).read() != txt:
                print("DIFFERS: %s" % fn)
                bad = 1
        else:
            open(p, "w").write(txt)
            n = sum(1 for l in txt.splitlines() if l.strip() and not l.startswith("#"))
            print("wrote %-16s %4d tasks" % (fn, n))
    sys.exit(bad)
