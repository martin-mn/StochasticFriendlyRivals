#!/usr/bin/env python3
"""
Pack the cw run into the two files Figures/cwdata.py reads (formats defined in its docstring), with no numpy,
so it runs on Cannon as it is.

    python3 pack.py run1 kappa packed           -> packed/cw1_e4.tsv, packed/cw1_e4.rep.tsv, packed/kappa_e4.tsv
    python3 pack.py run2 kappa packed --run=F2  -> packed/cw2_e4.tsv, packed/cw2_e4.rep.tsv, packed/cw2_e4.diag.tsv
                                                  (run F2: N = 100, beta = 3, u = 1e-4, itend = 1e8, nper 10).
                                                  kappa_e4.tsv is SHARED with F1 (the kernel and the games do not
                                                  depend on N): F2 reads and checks the same kappa directory, and
                                                  writes kappa_e4.tsv only if <out> has none; if <out> has one, its
                                                  data lines must equal what this kappa directory gives, byte for
                                                  byte, or the pack is refused
    python3 pack.py <rundir> <kappadir> <out> --dev
                                                  the same for a SHORT LOCAL run: itend/nper are taken from each
                                                  h-file's own header and kappa's M from its lines, and the first
                                                  line of every file is '# SHORT LOCAL RUN ...', which cwdata.py
                                                  refuses as production data (and requires in --short mode)

    cw1_e4.tsv       one row per cell, ipt 0..511:
                     ipt gid u v emax pay ef s000 s100 s010 s001 s110 s011 s111 s101 meanlive
                     s<ENR> = h-file atom count / '# total' (all 10 replicates summed), %.17g
                     pay = mean of the 10 replicate pays, ef = pay/emax (checked against the mean of the
                     replicate efs), meanlive = mean of the 10 replicate means
    cw1_e4.rep.tsv   one row per replicate: ipt irep pay ef meanlive solves_per_gen s000 .. s101
    kappa_e4.tsv     one row per cell: ipt k000 k100 k010 k001 k110 k011 k111 k101 z00 .. z80
                     (kappa.x's atom fractions, reordered from index 4E+2N+R to the named codes; its 81
                     stratum fractions unchanged)

The atom index in every cw/kappa output is 4E+2N+R.  The column for code 'ENR' is index int('ENR', 2):
    000->0  100->4  010->2  001->1  110->6  011->3  111->7  101->5

CHECKED, CELL BY CELL, not assumed:
  * TASK DECODE: the h-file's '# cell:' line names isl, gid and tag; isl must be 3*ipt+3 and gid/tag must
    equal games.dat row (isl-1)/3+1.  kappa's line names isl and gid; the same test.  ie must be 3 (eps 1e-4).
  * the game: (u, v) = (cS, cT-1) of that games.dat row, which must also equal its display columns
  * the run parameters in the h-file header (umut 1e-2, beta 100, eps 1e-4; nper 10, itend 1e7 unless --dev)
  * the h-file: exactly 8 atom, 81 stratum, 16 corner count lines; '# total' present; atoms, strata and
    total all equal nper*(itend/2)*N; corners equal the corner strata; the 101 count exactly 0
  * the summary file: exactly nper replicate lines, each naming this gid, each with 8 shares summing to 1
  * kappa: 98 fields, both groups summing to 1, k101 exactly 0
Any failure is listed and the exit status is 1; the files are still written for inspection.
"""
import os
import sys

N = 1000
NPER, ITEND = 10, 10 ** 7
UMUT, BETA, EPS = 1e-2, 100.0, 1e-4
MPROD = 10 ** 7
# the two runs.  F1 = the module constants above (and its outputs are byte-identical to the pre-F2 pack.py);
# F2 = DiskM1WF's sf2.f parameters (cw2.f = mkcw2.py(cw.f)).  eps, nper, the games and kappa's M are shared.
RUNS = {
    "F1": dict(N=N, NPER=NPER, ITEND=ITEND, UMUT=UMUT, BETA=BETA, name="cw1", sf="sf1.f", ulab="1e-2",
               seed=40000000, kappa_own=True),
    "F2": dict(N=100, NPER=10, ITEND=10 ** 8, UMUT=1e-4, BETA=3.0, name="cw2", sf="sf2.f", ulab="1e-4",
               seed=42000000, kappa_own=False),
}
ORDER8 = ["000", "100", "010", "001", "110", "011", "111", "101"]
IDX8 = [int(c, 2) for c in ORDER8]
CORNER_STRATUM = [sum(((k >> j) & 1) * 3 ** j for j in range(4)) for k in range(16)]


def games(path):
    """row (1-based) -> (gid, cR, cS, cT, u, v, tag)"""
    out, n = {}, 0
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        n += 1
        out[n] = (int(f[0]), float(f[2]), float(f[3]), float(f[4]), float(f[5]), float(f[6]), f[7])
    return out


def read_h(hp):
    H = dict(sec={}, cell=None, nper=None, itend=None, ie=None, total=None)
    sec = None
    for line in open(hp):
        if line.startswith("# cell:"):
            f = line.split()
            # '# cell: isl <isl>  gid <gid> <tag>  eps <e>  umut <u>  beta <b>'
            H["cell"] = (int(f[3]), int(f[5]), f[6], float(f[8]), float(f[10]), float(f[12]))
        elif line.startswith("# run:"):
            f = line.split()
            H["nper"], H["itend"], H["ie"] = int(f[3]), int(f[5]), int(f[7])
        elif line.startswith("# atoms"):
            sec = "a"; H["sec"][sec] = []
        elif line.startswith("# strata"):
            sec = "s"; H["sec"][sec] = []
        elif line.startswith("# corners"):
            sec = "c"; H["sec"][sec] = []
        elif line.startswith("# total"):
            H["total"] = int(line.split()[2]); sec = None
        elif line.startswith("#") or not line.strip():
            continue
        else:
            if sec is None:
                raise ValueError("%s: count line outside a section" % hp)
            f = line.split()
            H["sec"][sec].append((int(f[0]), int(f[1])))
    return H


NEARTIE_MAX = 1e-6      # share of sampled abundance on flag-1 lineages above which a cell is refused


def main(run, kdir, outdir, dev, allow_neartie=False, which="F1"):
    R = RUNS[which]
    N, NPER, ITEND, UMUT, BETA = R["N"], R["NPER"], R["ITEND"], R["UMUT"], R["BETA"]
    G = games(os.path.join(run, "games.dat"))
    assert len(G) == 512, "games.dat in %s has %d rows, expected 512" % (run, len(G))
    cells, reps, kap, bad, diag = [], [], [], [], []
    seen = set()
    for ipt in range(512):
        isl, row = 3 * ipt + 3, ipt + 1
        gid, cR, cS, cT, ud, vd, tag = G[row]
        u, v = cS, cT - 1.0
        if cR != 1.0 or (u, v) != (ud, vd):
            bad.append("%d: games.dat row %d: (cR,cS,cT) and (u,v) columns disagree" % (isl, row)); continue
        hp, sp = os.path.join(run, "h%d" % isl), os.path.join(run, "%d" % isl)
        if not (os.path.exists(hp) and os.path.exists(sp)):
            bad.append("%d: missing product" % isl); continue
        H = read_h(hp)
        if H["cell"] is None or H["nper"] is None:
            bad.append("%d: h-file has no header" % isl); continue
        hisl, hgid, htag, heps, humut, hbeta = H["cell"]
        if (hisl, hgid, htag) != (isl, gid, tag):
            bad.append("%d: h-file says isl %d gid %d %s, games.dat row %d is gid %d %s -- TASK DECODE"
                       % (isl, hisl, hgid, htag, row, gid, tag)); continue
        if H["ie"] != 3 or abs(heps - EPS) > 1e-12 or abs(humut - UMUT) > 1e-12 or abs(hbeta - BETA) > 1e-9:
            bad.append("%d: ran with ie %d eps %g umut %g beta %g" % (isl, H["ie"], heps, humut, hbeta)); continue
        nper, itend = H["nper"], H["itend"]
        if not dev and (nper, itend) != (NPER, ITEND):
            bad.append("%d: nper %d itend %d, production is %d / %d" % (isl, nper, itend, NPER, ITEND)); continue
        seen.add((nper, itend))
        if H["total"] is None:
            bad.append("%d: no '# total' line -- h-file TRUNCATED" % isl); continue
        a, s, c = (H["sec"].get(k, []) for k in "asc")
        if [k for k, _ in a] != list(range(8)) or [k for k, _ in s] != list(range(81)) \
                or [k for k, _ in c] != list(range(16)):
            bad.append("%d: h-file sections have %d/%d/%d lines, expected 8/81/16" % (isl, len(a), len(s), len(c)))
            continue
        a, s, c = [x for _, x in a], [x for _, x in s], [x for _, x in c]
        expect = nper * (itend // 2) * N
        if not (sum(a) == sum(s) == H["total"] == expect):
            bad.append("%d: atoms %d strata %d total %s expected %d" % (isl, sum(a), sum(s), H["total"], expect))
            continue
        if any(c[k] != s[CORNER_STRATUM[k]] for k in range(16)):
            bad.append("%d: corner counts differ from the corner strata" % isl); continue
        if a[5] != 0:
            bad.append("%d: atom 101 holds %d -- the theorem or the classifier fails" % (isl, a[5])); continue
        rl = [l.split() for l in open(sp) if l.strip()]
        if len(rl) != nper:
            bad.append("%d: %d replicate lines, expected %d" % (isl, len(rl), nper)); continue
        emax = max(cR, 0.5 * (cS + cT))
        ok = True
        for irep, f in enumerate(rl, 1):
            if len(f) != 17 or int(f[0]) != gid or int(f[1]) != 3:
                bad.append("%d: replicate %d malformed or names gid %s" % (isl, irep, f[0])); ok = False; break
            sh = [float(x) for x in f[9:17]]
            if abs(sum(sh) - 1) > 1e-9:
                bad.append("%d: replicate %d shares sum to %.12g" % (isl, irep, sum(sh))); ok = False; break
            if abs(float(f[4]) - emax) > 1e-9:
                bad.append("%d: replicate %d Emax %s, games.dat gives %.10f" % (isl, irep, f[4], emax)); ok = False
                break
            reps.append((ipt, irep, float(f[5]), float(f[6]), float(f[7]), float(f[8])) + tuple(sh[i] for i in IDX8))
        if not ok:
            continue
        pay = sum(float(f[5]) for f in rl) / nper
        efm = sum(float(f[6]) for f in rl) / nper
        ef = pay / emax
        if abs(ef - efm) > 1e-9:
            bad.append("%d: mean ef %.12g but mean pay / emax %.12g" % (isl, efm, ef)); continue
        live = sum(float(f[7]) for f in rl) / nper
        # the classifier's own diagnostics, one '# ' line per replicate in w<isl> (cw.f [cw 6]):
        #   isl irun ncls ncall nneartie neartie_abund nhw hit meanlive solves gcount nexacttie exacttie_abund
        # neartie = m1cls flag 1, a Nash margin within 1e-9*scale of 0 decided by the tau band: the only
        # classifications that rounding could flip.  They must be visible, not silent.
        wp = os.path.join(run, "w%d" % isl)
        if not os.path.exists(wp):
            bad.append("%d: w-file missing" % isl); continue
        wl = [l.split()[1:] for l in open(wp) if l.startswith("# ") and len(l.split()) == 14
              and l.split()[1].isdigit()]
        if len(wl) != nper or any(int(w[0]) != isl for w in wl):
            bad.append("%d: w-file has %d diagnostic lines, expected %d" % (isl, len(wl), nper)); continue
        nnt, ntab = sum(int(w[4]) for w in wl), sum(int(w[5]) for w in wl)
        diag.append((ipt, sum(int(w[2]) for w in wl), sum(int(w[3]) for w in wl), nnt, ntab,
                     ntab / H["total"], max(int(w[6]) for w in wl), sum(int(w[12]) for w in wl) / H["total"]))
        if ntab / H["total"] > NEARTIE_MAX and not allow_neartie:
            bad.append("%d: near-tie (flag 1) abundance share %.3g exceeds %g -- inspect before plotting"
                       % (isl, ntab / H["total"], NEARTIE_MAX)); continue
        cells.append((ipt, gid, u, v, emax, pay, ef, [a[i] / H["total"] for i in IDX8], live))

        kp = os.path.join(kdir, "k%d" % isl)
        if not os.path.exists(kp):
            bad.append("%d: kappa file missing" % isl); continue
        kl = [l.split() for l in open(kp) if l.strip()]
        if len(kl) != 1 or len(kl[0]) != 98:
            bad.append("%d: kappa file has %d lines / %s fields, expected 1 / 98"
                       % (isl, len(kl), len(kl[0]) if kl else 0)); continue
        f = kl[0]
        if (int(f[0]), int(f[1]), int(f[2])) != (isl, gid, 3) or (float(f[4]), float(f[5])) != (u, v):
            bad.append("%d: kappa line says isl %s gid %s ie %s (u,v) (%s,%s) -- TASK DECODE" % ((isl,) + tuple(f[:6])))
            continue
        M = int(f[6])
        if not dev and M != MPROD:
            bad.append("%d: kappa M = %d, production is %d" % (isl, M, MPROD)); continue
        ka = [float(x) for x in f[7:15]]
        kz = [float(x) for x in f[15:96]]
        if abs(sum(ka) - 1) > 1e-9 or abs(sum(kz) - 1) > 1e-9 or ka[5] != 0:
            bad.append("%d: kappa sums %.12g %.12g, k101 %g" % (isl, sum(ka), sum(kz), ka[5])); continue
        kap.append((ipt, [ka[i] for i in IDX8], kz, M))

    Ms = sorted({k[3] for k in kap})
    if len(seen) > 1:
        bad.append("cells ran with different (nper, itend): %s" % sorted(seen))
    if len(Ms) > 1:
        bad.append("kappa lines have different M: %s" % Ms)
    nper, itend = sorted(seen)[0] if seen else (NPER, ITEND)
    M = Ms[0] if Ms else MPROD
    mark = []
    if dev:
        mark = ["# SHORT LOCAL RUN -- itend = %d, M = %d, gfortran on the Mac: a test of the chain run -> pack -> "
                "figure, NOT a result" % (itend, M)]
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, "%s_e4.tsv" % R["name"])
    with open(p, "w") as fh:
        fh.write("\n".join(mark + [
            "# PartnersRivalsM1f %s: the WF process (%s's) over the CONTINUOUS memory-one strategy space,"
            % (R["name"], R["sf"]),
            "# exact-alpha kernel (alpha = 1/3 per face), 512 Voronoi cells of BinM2Ev ca6's sunflower; kit cannon/cw",
            "# s<ENR> = time-averaged share of the population in atom (E,N,R) of the exact eps->0 classification",
            "# N = %d, beta = %d, u = %s, eps = 1e-4, itend = %d, %d replicates, iseedb = %d"
            % (N, BETA, R["ulab"], itend, nper, R["seed"]),
            "# " + "\t".join(["ipt", "gid", "u", "v", "emax", "pay", "ef"] + ["s" + c for c in ORDER8] + ["meanlive"])])
            + "\n")
        for ipt, gid, u, v, emax, pay, ef, sh, live in cells:
            fh.write("\t".join(["%d" % ipt, "%d" % gid, "%.17g" % u, "%.17g" % v, "%.17g" % emax, "%.17g" % pay,
                                "%.17g" % ef] + ["%.17g" % x for x in sh] + ["%.6f" % live]) + "\n")
    q = os.path.join(outdir, "%s_e4.rep.tsv" % R["name"])
    with open(q, "w") as fh:
        fh.write("\n".join(mark + ["# " + "\t".join(["ipt", "irep", "pay", "ef", "meanlive", "solves_per_gen"]
                                                   + ["s" + c for c in ORDER8]) + "    (one line per replicate)"])
                 + "\n")
        for r in reps:
            fh.write("%d\t%d\t%.10f\t%.10f\t%.4f\t%.4f\t" % r[:6] + "\t".join("%.12f" % x for x in r[6:]) + "\n")
    k = os.path.join(outdir, "kappa_e4.tsv")
    ktxt = "\n".join(mark + [
        "# PartnersRivalsM1f kappa: the exact-alpha kernel's mass on each atom and stratum per game (kappa.x)",
        "# M = %d, seed = 41000000 + isl" % M,
        "# " + "\t".join(["ipt"] + ["k" + c for c in ORDER8] + ["z%02d" % i for i in range(81)])]) + "\n"
    ktxt += "".join("\t".join(["%d" % ipt] + ["%.12g" % x for x in ka] + ["%.12g" % x for x in kz]) + "\n"
                    for ipt, ka, kz, _ in kap)
    kmsg = None
    if R["kappa_own"]:
        with open(k, "w") as fh:
            fh.write(ktxt)
    elif not os.path.exists(k):
        # F2 may create the SHARED file only from a complete, consistent pack: a partial one would later
        # make a clean pack look like a 'kappa mismatch' (audit F2, nonblocking 3)
        if bad:
            kmsg = "%s: NOT written -- this F2 pack is incomplete, and kappa_e4.tsv is shared with F1" % k
        else:
            with open(k, "w") as fh:
                fh.write(ktxt)
    else:
        # F2 shares F1's kappa_e4.tsv: leave it untouched, but only if it is exactly what this kappa directory
        # gives (the '# SHORT LOCAL RUN' mark names the run's itend, so it is excluded from the comparison)
        def body(t):
            return [l for l in t.split("\n") if not l.startswith("# SHORT LOCAL RUN")]
        if body(open(k).read()) != body(ktxt):
            bad.append("%s exists and differs from what %s gives -- kappa_e4.tsv is shared by F1 and F2 and "
                       "must be one kernel measurement" % (k, kdir))
        else:
            kmsg = "%s: already present, and identical (data lines) to %s -- shared with F1, left untouched" % (k, kdir)
    d = os.path.join(outdir, "%s_e4.diag.tsv" % R["name"])
    with open(d, "w") as fh:
        fh.write("\n".join(mark + [
            "# per-cell classifier diagnostics, summed over replicates (from the w-files)",
            "# neartie = m1cls flag 1 (Nash margin within 1e-9*scale of 0, decided by the tau band);"
            " exacttie = flag 2 (margin computed as exactly 0.0, in practice a structural tie)",
            "# " + "\t".join(["ipt", "ncls", "ncall", "nneartie", "neartie_abund", "neartie_share", "nhw_max",
                               "exacttie_share"])]) + "\n")
        for r in diag:
            fh.write("%d\t%d\t%d\t%d\t%d\t%.6g\t%d\t%.6g\n" % r)
    nt = [r for r in diag if r[3]]
    print("%s: %d cells; flag-1 near ties in %d cells, total %d lineages, max abundance share %.3g"
          % (d, len(diag), len(nt), sum(r[3] for r in diag), max([r[5] for r in diag] or [0])))
    print("%s: %d cells of 512; %s: %d replicate lines; %s: %d cells of 512"
          % (p, len(cells), q, len(reps), k, len(kap)))
    if kmsg:
        print(kmsg)
    if bad:
        print("INCOMPLETE OR INCONSISTENT (%d):" % len(bad))
        for b in bad[:40]:
            print("   ", b)
        sys.exit(1)
    print("all 512 cells complete and consistent (task decode checked against games.dat at every cell)%s%s"
          % ("  [--dev: itend %d, M %d]" % (itend, M) if dev else "", "" if which == "F1" else "  [run %s]" % which))


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    runs = [x.split("=", 1)[1] for x in sys.argv[1:] if x.startswith("--run=")]
    which = runs[-1] if runs else "F1"
    if which not in RUNS:
        sys.exit("--run=%s: the runs are %s" % (which, ", ".join(sorted(RUNS))))
    main(a[0] if a else ("run1" if which == "F1" else "run2"), a[1] if len(a) > 1 else "kappa",
         a[2] if len(a) > 2 else "packed", "--dev" in sys.argv, "--allow-neartie" in sys.argv, which)
