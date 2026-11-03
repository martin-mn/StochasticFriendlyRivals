#!/usr/bin/env python3
"""collect.py -- PRFull cannon/cw2: tables for the E3, E5 and TS arms, and their comparison with production.

    python3 collect.py                      every arm with output under out/ (E3, E5, TS); tables to results/
    python3 collect.py E3 E5 TS TSX         just these (TSX: the optional memory-one split games)
    python3 collect.py --local              the local smoke test out/local (short itend; NOT a result)
    python3 collect.py --kappa-check D3 D5  (pure python) kernel masses at eps 1e-3 / 1e-5 vs production kappa_e4
    python3 collect.py --probe-times LOGDIR PREFIX   (pure python) probe rates and the --time each arm needs

Production comparator (read only, copies in ref/packed/, sha256 in ref/SHA256SUMS): PartnersRivalsM1f
cannon/cw/packed/cw1_e4.rep.tsv (the ten replicates of run F1, eps = 1e-4, at all 512 games) and kappa_e4.tsv
(kernel masses, M = 1e7 draws per game).  The atoms are written in the paper's ORDER 000 100 010 001 110 011 111
(E N R: efficient, stable, competitive; 101 is empty).  Paper conventions (PRFull Figures/fig45.py census):
  efficiency  = mean per-round payoff / Emax (the mean of the ten replicate values; each replicate's own line)
  share       = time-averaged fraction of the population in the atom over the sampled half, ten replicates pooled
  most abundant = argmax share;  most enriched = argmax share/kappa over atoms with kappa > 0 (kappa = 0 but share
                  > 0 floored at 3/M, as M1f's cwdata);  friendly rivals = 111;  competitive = 001 + 011 + 111
  replicate agreement = number of the ten replicates whose own most abundant atom is the pooled one
Medians per wedge get a 95 % bootstrap interval over replicates within games (4000 resamples, fixed seed), as
the triage's rare-event intervals (ref1_triage report 6).

Success criteria of the eps arms (triage, issue 7, after Revised2:801), evaluated and printed:
  (1) the most abundant atom unchanged at >= 95 % of the games;  (2) efficiency within replicate spread
  (|Welch t| <= 2 at >= 95 % of games, none beyond 3);  (3) friendly rivals still lead every S game;
  (4) the W/N/E friendly-rival and competitive medians inside production's bootstrap intervals.
Time-series diagnostics (triage, issue 8): Q3 vs Q4 of the production window (5e6,1e7] with batch-means SE
(20 batches), Geyer's initial-positive-sequence integrated autocorrelation time and ESS, the run-length check
(5e6,1e7] vs (2e7,4e7], friendly-rival takeover and escape counts with exact intervals.  Success: at every game
whose production replicates agree on the leader, the quarters agree (at most 2 of 10 replicates beyond 2 SE),
the ESS of the leading-atom share over the production window is >= 50 in every replicate, and the leader and efficiency agree
between the two windows within replicate spread.
"""
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.join(HERE, "ref", "packed")
ORDER = ["000", "100", "010", "001", "110", "011", "111"]
RAW = [0, 4, 2, 1, 6, 3, 7]               # raw index 4E+2N+R of each ORDER atom (raw 5 = 101, empty)
EFF, STAB, COMP = ("100", "110", "111"), ("010", "110", "011", "111"), ("001", "011", "111")
NPOP = 1000
WEDGES = ["W", "S", "N", "E"]


# ------------------------------------------------------------------------------------------------ plain helpers
def games():
    g = {}
    for l in open(os.path.join(HERE, "games.dat")):
        if l.lstrip().startswith("#") or not l.strip():
            continue
        f = l.split()
        tag = f[7]
        g[int(tag.split("_")[1][1:])] = dict(gid=int(f[0]), tag=tag, wedge=tag[-1], quad=tag[-3:-1],
                                             u=float(f[5]), v=float(f[6]))
    assert len(g) == 512
    return g


def subset(fn):
    return [int(l.split()[0]) for l in open(os.path.join(HERE, fn)) if l.strip() and not l.startswith("#")]


def kappa_prod():
    K = {}
    for l in open(os.path.join(PACK, "kappa_e4.tsv")):
        if l.startswith("#") or not l.strip():
            continue
        f = l.split("\t")
        K[int(f[0])] = [float(x) for x in f[1:9]]          # ORDER + 101
    return K, 10000000


def kappa_dir(d):
    """k<isl> files of kappa.x: isl gid ie eps u v M f_raw[0..7] ... -> {ipt: ORDER+101 fractions}, M"""
    K, M = {}, None
    for fn in glob.glob(os.path.join(d, "k*")):
        f = open(fn).read().split()
        if len(f) != 98:
            continue
        isl, m = int(f[0]), int(f[6])
        raw = [float(x) for x in f[7:15]]
        K[(isl - 1) // 3] = [raw[i] for i in RAW] + [raw[5]]
        M = m if M is None else M
        assert M == m
    return K, M


# ------------------------------------------------------------------------------------------------ --kappa-check
def kappa_check(dirs):
    Kp, Mp = kappa_prod()
    worst, nbad, nall = 0.0, 0, 0
    for d in dirs:
        K, M = kappa_dir(d)
        if not K:
            print("   %s: no kappa files" % d)
            nbad += 1
            continue
        zmax, where = 0.0, None
        for ipt, k in sorted(K.items()):
            for j in range(8):
                p1, p2 = k[j], Kp[ipt][j]
                pb = (p1 * M + p2 * Mp) / (M + Mp)
                se = math.sqrt(max(pb * (1 - pb), 0.0) * (1.0 / M + 1.0 / Mp))
                z = abs(p1 - p2) / se if se > 0 else (0.0 if p1 == p2 else float("inf"))
                nall += 1
                if z > zmax:
                    zmax, where = z, (ipt, (ORDER + ["101"])[j], p1, p2)
                if j == 7 and p1 != 0.0:
                    print("   %s ipt %d: atom 101 has mass %g (must be 0)" % (d, ipt, p1))
                    nbad += 1
        print("   %s: %d games, M = %d; max |z| = %.2f at ipt %d atom %s (%.6g vs production %.6g)"
              % (d, len(K), M, zmax, where[0], where[1], where[2], where[3]))
        worst = max(worst, zmax)
    ok = worst < 5.0 and nbad == 0
    print("   kernel masses %s eps-independent within sampling (max |z| %.2f over %d comparisons; threshold 5)"
          % ("ARE" if ok else "are NOT", worst, nall))
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------ --probe-times
def probe_times(logdir, prefix):
    rows = [l.split("#")[0].split() for l in open(os.path.join(HERE, "tasks_probe.txt"))
            if l.strip() and not l.startswith("#")]
    rate = {}
    logs = glob.glob(os.path.join(logdir, prefix + "*.out"))
    for fn in sorted(logs, key=lambda x: int(x.rsplit("_", 1)[1].split(".")[0])):
        k = int(fn.rsplit("_", 1)[1].split(".")[0])
        r = rows[k - 1]
        for l in open(fn):
            if "complete in" in l and "s wall" in l and r[0] != "KA":
                w = int(l.split(" s wall")[0].split()[-1])
                itend = int(r[8])
                per1e7 = w / (itend / 1e7) / (int(r[5]) - int(r[4]) + 1)
                rate.setdefault(r[0], []).append((w, itend, per1e7, r[3]))
                print("      %-6s ipt %4s itend %9d: %5d s wall = %6.1f s per 1e7 generations per replicate"
                      % (r[0], r[3], itend, w, per1e7))
    def worst(keys):
        v = [x[2] for k in keys for x in rate.get(k, [])]
        return max(v) if v else None
    for arm, keys, gens, reps, cur in (("E3", ("PE3", "PE3s"), 1e7, 10, 105), ("E5", ("PE5", "PE5s"), 1e7, 10, 105),
                                       ("TS", ("PTS", "PTSs"), 4e7, 1, 50)):
        w = worst(keys)
        if w is None:
            print("   %s: no probe timing" % arm)
            continue
        task = w * gens / 1e7 * reps
        need = 2 * task / 60
        print("   %s: worst %.1f s per 1e7 per replicate -> task %.0f s; 2x margin = %.0f min (sbatch has %d min)%s"
              % (arm, w, task, need, cur, "  ** RAISE --time **" if need > cur else ""))
    return 0


# ------------------------------------------------------------------------------------------------ the analysis (numpy)
def analysis(arms, local=False):
    import numpy as np
    rng = np.random.default_rng(20260927)
    G = games()
    root = os.path.join(HERE, "out", "local") if local else os.path.join(HERE, "out")
    resd = os.path.join(HERE, "results_local" if local else "results")
    os.makedirs(resd, exist_ok=True)
    B = 4000

    # production, per replicate
    Rp, Ep = {}, {}
    for l in open(os.path.join(PACK, "cw1_e4.rep.tsv")):
        if l.startswith("#") or not l.strip():
            continue
        f = l.split("\t")
        i = int(f[0])
        s = [float(x) for x in f[6:14]]
        assert s[7] == 0.0
        Rp.setdefault(i, []).append(s[:7])
        Ep.setdefault(i, []).append(float(f[3]))
    Rp = {i: np.array(v) for i, v in Rp.items()}
    Ep = {i: np.array(v) for i, v in Ep.items()}
    assert all(v.shape == (10, 7) for v in Rp.values()) and len(Rp) == 512
    Kp, Mp = kappa_prod()
    iprod_lead = {i: int(Rp[i].mean(0).argmax()) for i in Rp}
    prod_split = {i: bool((Rp[i].argmax(1) != iprod_lead[i]).any()) for i in Rp}

    def per_game(R, E, K, M):
        S = R.mean(0)
        Kf = np.where((K == 0) & (S > 0), 3.0 / M, K)
        enr = np.where(Kf > 0, S / np.where(Kf > 0, Kf, 1), -np.inf)
        lr = R.argmax(1)
        wa = int(S.argmax())
        return dict(E=float(E.mean()), sdE=float(E.std(ddof=1)), S=S, wa=wa, we=int(enr.argmax()),
                    nagree=int((lr == wa).sum()), eff=float(S[[ORDER.index(a) for a in EFF]].sum()),
                    stab=float(S[[ORDER.index(a) for a in STAB]].sum()),
                    comp=float(S[[ORDER.index(a) for a in COMP]].sum()), s111=float(S[6]))

    def boot_medians(Rs, idxs, fn):
        """95 % bootstrap interval of the median over the games idxs of fn(game mean shares), resampling the ten
        replicates within each game"""
        A = np.stack([fn(Rs[i]) for i in idxs])                   # (G, 10): the statistic per replicate
        pick = rng.integers(0, A.shape[1], size=(B, A.shape[0], A.shape[1]))
        M = np.take_along_axis(A[None, :, :], pick, axis=2).mean(2)          # (B, G): resampled game means
        v = np.median(M, axis=1)
        return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))

    f111 = lambda M: M[..., 6]
    fcomp = lambda M: M[..., [ORDER.index(a) for a in COMP]].sum(-1)

    report = []
    P = report.append

    def census(name, D, Rs, idxs):
        P("  %s" % name)
        P("    %-6s %5s %8s %8s %8s | most abundant %s | most enriched %s | %s | %s | all-10-agree"
          % ("wedge", "games", "E min", "E med", "E max", " ".join("%4s" % a for a in ORDER),
             " ".join("%4s" % a for a in ORDER), "111 share median [95%]", "competitive median [95%]"))
        for w in WEDGES + ["all"]:
            ii = [i for i in idxs if w == "all" or G[i]["wedge"] == w]
            if not ii:
                continue
            e = np.array([D[i]["E"] for i in ii])
            wa = np.array([D[i]["wa"] for i in ii])
            we = np.array([D[i]["we"] for i in ii])
            m111 = np.median([D[i]["s111"] for i in ii])
            mcomp = np.median([D[i]["comp"] for i in ii])
            b1 = boot_medians(Rs, ii, f111)
            b2 = boot_medians(Rs, ii, fcomp)
            P("    %-6s %5d %8.4f %8.4f %8.4f | %s | %s | %.3g [%.3g, %.3g] | %.4g [%.4g, %.4g] | %d"
              % (w, len(ii), e.min(), np.median(e), e.max(), " ".join("%4d" % (wa == j).sum() for j in range(7)),
                 " ".join("%4d" % (we == j).sum() for j in range(7)), m111, b1[0], b1[1], mcomp, b2[0], b2[1],
                 sum(D[i]["nagree"] == 10 for i in ii)))
        return

    # ---------------------------------------------------------------------------- the eps arms
    for arm in [a for a in arms if a in ("E3", "E5")]:
        ie, eps = (2, "1e-3") if arm == "E3" else (1, "1e-5")
        d = os.path.join(root, arm)
        if not os.path.isdir(d):
            print("%s: no output directory %s" % (arm, d))
            continue
        idxs = subset("games71.txt") if not local else sorted({(int(os.path.basename(p)) - 1) // 3 for p in
                                                               glob.glob(os.path.join(d, "[0-9]*"))})
        Ra, Ea, missing = {}, {}, []
        for i in idxs:
            fn = os.path.join(d, str(3 * i + ie))
            if not os.path.exists(fn):
                missing.append(i)
                continue
            L = [l.split() for l in open(fn) if l.strip()]
            if len(L) != 10:
                missing.append(i)
                continue
            for f in L:
                assert int(f[0]) == G[i]["gid"] and int(f[1]) == ie, (fn, f[:2])
                s = [float(x) for x in f[9:17]]
                assert s[5] == 0.0 and abs(sum(s) - 1) < 1e-9, fn
            Ra[i] = np.array([[float(f[9 + RAW[j]]) for j in range(7)] for f in L])
            Ea[i] = np.array([float(f[6]) for f in L])
        if missing:
            print("%s: %d games without 10 replicate lines: %s" % (arm, len(missing), missing[:20]))
        idxs = [i for i in idxs if i in Ra]
        if not idxs:
            continue
        # [verifier fix] ONE denominator for the arm and for production.  The kernel's law does not depend on eps
        # (atoms are read on p = (s-eps)/(1-2eps), uniform whatever eps), so production's kappa_e4 is the arm's
        # kernel mass too.  Dividing the arm by its own independent 1e7-draw estimate (KA) and production by
        # kappa_e4 would let Monte-Carlo noise in kappa -- worst for atoms of a few draws, where the 3/M floor
        # applies to one estimate and not the other -- flip "most enriched" with no change in any share.  KA is
        # the eps-independence check (verify.sh probe, --kappa-check); here it only reports how many games' most
        # enriched atom would change if the arm's own estimate were used instead.
        kd = os.path.join(HERE, "out", "probe", "KA3" if arm == "E3" else "KA5")
        Ka, Ma = kappa_dir(kd) if os.path.isdir(kd) else ({}, None)
        have_ka = bool(idxs) and all(i in Ka for i in idxs)
        ksrc = ("production's kappa_e4 (eps 1e-4, M = 1e7) for the arm AND for production: the kernel's law does not "
                "depend on eps")
        Da = {i: per_game(Ra[i], Ea[i], np.array(Kp[i][:7]), Mp) for i in idxs}
        Dp = {i: per_game(Rp[i], Ep[i], np.array(Kp[i][:7]), Mp) for i in idxs}
        if have_ka:
            nka = sum(per_game(Ra[i], Ea[i], np.array(Ka[i][:7]), Ma)["we"] != Da[i]["we"] for i in idxs)
            ksrc += "; with the arm's own KA masses (%s, M = %s) the most enriched atom would differ at %d of %d games" \
                    % (os.path.relpath(kd, HERE), Ma, nka, len(idxs))
        P("=" * 118)
        P("%s: eps = %s, memory-one cube, F1 regime (N = 1000, beta = 100, u = 1e-2, 10 x 1e7 generations), %d games%s"
          % (arm, eps, len(idxs), "  [LOCAL SMOKE TEST, NOT A RESULT]" if local else ""))
        P("  kernel masses: %s" % ksrc)
        census("%s (eps = %s)" % (arm, eps), Da, Ra, idxs)
        census("production at the same games (eps = 1e-4)", Dp, Rp, idxs)
        # comparison
        tE, dE, same_a, same_e = [], [], 0, 0
        rows = []
        for i in idxs:
            a, p = Da[i], Dp[i]
            se = math.sqrt(a["sdE"] ** 2 / 10 + p["sdE"] ** 2 / 10)
            t = (a["E"] - p["E"]) / se if se > 0 else (0.0 if a["E"] == p["E"] else float("inf"))
            tE.append(t)
            dE.append(a["E"] - p["E"])
            same_a += a["wa"] == p["wa"]
            same_e += a["we"] == p["we"]
            rows.append((i, G[i]["tag"], G[i]["wedge"], a["E"], p["E"], a["E"] - p["E"], t, ORDER[a["wa"]], ORDER[p["wa"]],
                         ORDER[a["we"]], ORDER[p["we"]], a["nagree"], p["nagree"], a["s111"], p["s111"], a["comp"],
                         p["comp"], a["eff"], p["eff"]) + tuple(a["S"]) + tuple(p["S"]))
        tE, dE = np.array(tE), np.array(dE)
        n = len(idxs)
        P("  comparison with production at the same %d games" % n)
        P("    most abundant atom the same at %d of %d games (%.1f %%)%s; most enriched the same at %d"
          % (same_a, n, 100.0 * same_a / n, "" if same_a == n else ": differs at " + ", ".join(
              "%d (%s -> %s)" % (i, ORDER[Dp[i]["wa"]], ORDER[Da[i]["wa"]]) for i in idxs if Da[i]["wa"] != Dp[i]["wa"]),
             same_e))
        P("    efficiency: median dE %+.2e, max |dE| %.2e; |Welch t| > 2 at %d, > 3 at %d of %d games (5 %% expected by "
          "chance: %.1f)" % (np.median(dE), np.abs(dE).max(), (np.abs(tE) > 2).sum(), (np.abs(tE) > 3).sum(), n, 0.05 * n))
        sS = [i for i in idxs if G[i]["wedge"] == "S"]
        c1 = same_a >= 0.95 * n
        c2 = (np.abs(tE) <= 2).sum() >= 0.95 * n and (np.abs(tE) > 3).sum() == 0
        c3 = all(Da[i]["wa"] == 6 for i in sS)
        c4, c4txt = True, []
        for w in ("W", "N", "E"):
            ii = [i for i in idxs if G[i]["wedge"] == w]
            if not ii:
                continue
            for lab, key, fn in (("111", "s111", f111), ("competitive", "comp", fcomp)):
                lo, hi = boot_medians(Rp, ii, fn)
                m = float(np.median([Da[i][key] for i in ii]))
                inside = lo <= m <= hi
                c4 &= inside
                c4txt.append("%s %s %.3g %s [%.3g, %.3g]" % (w, lab, m, "in" if inside else "OUTSIDE", lo, hi))
        P("    (1) most abundant unchanged at >= 95 %% of games: %s" % ("PASS" if c1 else "FAIL"))
        P("    (2) efficiency within replicate spread: %s" % ("PASS" if c2 else "FAIL"))
        P("    (3) friendly rivals lead all %d S games: %s (%d of %d)" % (len(sS), "PASS" if c3 else "FAIL",
                                                                        sum(Da[i]["wa"] == 6 for i in sS), len(sS)))
        P("    (4) W/N/E medians inside production's bootstrap intervals: %s -- %s" % ("PASS" if c4 else "FAIL",
                                                                                   "; ".join(c4txt)))
        hdr = ("ipt tag wedge E_arm E_prod dE welch_t lead_arm lead_prod enr_arm enr_prod agree_arm agree_prod "
               "s111_arm s111_prod comp_arm comp_prod effsh_arm effsh_prod " + " ".join("a%s" % a for a in ORDER) + " " +
               " ".join("p%s" % a for a in ORDER))
        with open(os.path.join(resd, "%s_games.tsv" % arm), "w") as fo:
            fo.write("# %s (eps %s) vs production (eps 1e-4), memory-one cube F1; shares in ORDER; kappa: %s\n"
                     % (arm, eps, ksrc))
            fo.write("\t".join(hdr.split()) + "\n")
            for r in rows:
                fo.write("\t".join(str(x) if not isinstance(x, float) else "%.10g" % x for x in r) + "\n")
        P("  per-game table: %s" % os.path.relpath(os.path.join(resd, "%s_games.tsv" % arm), HERE))

    # ---------------------------------------------------------------------------- the time series
    for arm, gf in (("TS", "games24.txt"), ("TSX", "games_tsx.txt")):
        if arm not in arms:
            continue
        d = os.path.join(root, arm)
        if os.path.isdir(d):
            ts_analysis(np, rng, G, d, resd, P, Rp, Ep, iprod_lead, prod_split, local, arm, gf)
        else:
            print("%s: no output directory %s" % (arm, d))

    txt = "\n".join(report)
    print(txt)
    with open(os.path.join(resd, "summary%s.txt" % ("_" + "_".join(arms))), "w") as fo:
        fo.write(txt + "\n")


def iat_geyer(np, x):
    """Geyer's initial positive sequence (with the monotone correction): tau in units of the series step;
    nan for a constant series"""
    x = np.asarray(x, float)
    n = len(x)
    xc = x - x.mean()
    g0 = float(xc @ xc) / n
    if g0 <= 0:
        return float("nan")
    f = np.fft.rfft(xc, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n] / n
    s, prev, k = 0.0, float("inf"), 0
    while 2 * k + 1 < n:
        G = ac[2 * k] + ac[2 * k + 1]
        if G <= 0:
            break
        G = min(G, prev)
        s += G
        prev = G
        k += 1
    return max((-g0 + 2 * s) / g0, 1.0)


def bm_se(np, x, nb=20):
    m = np.array([b.mean() for b in np.array_split(np.asarray(x, float), nb)])
    return float(m.std(ddof=1) / math.sqrt(nb))


def ts_analysis(np, rng, G, d, resd, P, Rp, Ep, iprod_lead, prod_split, local, arm="TS", gf="games24.txt"):
    from scipy.stats import beta as betad, chi2
    ipts = subset(gf)
    P("=" * 118)
    P("%s: time series, memory-one cube, F1 regime, eps = 1e-4, itend 4x production, games of %s%s"
      % (arm, gf, "  [LOCAL SMOKE TEST, NOT A RESULT]" if local else ""))
    rows, series = [], {}
    for i in ipts:
        dirs = sorted(glob.glob(os.path.join(d, "p%03d_r*" % i)))
        if not dirs:
            continue
        isl = 3 * i + 3
        effs, shs, reps, hdr = [], [], [], None
        for dd in dirs:
            fn = os.path.join(dd, "t%d" % isl)
            if not os.path.exists(fn):
                continue
            h, eff, C = None, [], []
            for l in open(fn):
                if l.startswith("# ts:"):
                    f = l.split()
                    h = dict(itend=int(f[3]), nblk=int(f[5]), lblk=int(f[7]))
                elif l.startswith("TS"):
                    f = l.split()
                    eff.append(float(f[4]))
                    C.append([int(x) for x in f[7:15]])
            C = np.array(C)
            if h is None or len(eff) != h["nblk"]:
                continue
            assert (C.sum(1) == NPOP * h["lblk"]).all(), fn
            hdr = h
            effs.append(eff)
            shs.append(C[:, RAW] / float(NPOP * h["lblk"]))
            reps.append(dd)
        if not effs:
            continue
        eff, sh = np.array(effs), np.array(shs)                     # (r, nblk), (r, nblk, 7)
        nb, L, it = hdr["nblk"], hdr["lblk"], hdr["itend"]
        itp = it // 4                                                # the production length inside the TS run
        blk = lambda a, b: slice(a // L, b // L)                     # generations (a, b] -> 0-based block slice
        Pw, Q3, Q4, Lw = blk(itp // 2, itp), blk(itp // 2, 3 * itp // 4), blk(3 * itp // 4, itp), blk(it // 2, it)
        lead = iprod_lead[i]
        series[i] = (eff, sh)
        nr = eff.shape[0]
        zE, zL, essL_P, essL_L, essE_P, iatL_P = [], [], [], [], [], []
        for r in range(nr):
            for x, zs in ((eff[r], zE), (sh[r, :, lead], zL)):
                a, b = x[Q3], x[Q4]
                se = math.hypot(bm_se(np, a), bm_se(np, b))
                zs.append((b.mean() - a.mean()) / se if se > 0 else (0.0 if a.mean() == b.mean() else float("inf")))
            tP = iat_geyer(np, sh[r, Pw, lead])
            tL = iat_geyer(np, sh[r, Lw, lead])
            tE = iat_geyer(np, eff[r, Pw])
            nP, nL = Pw.stop - Pw.start, Lw.stop - Lw.start
            essL_P.append(nP / tP if tP == tP else nP)
            essL_L.append(nL / tL if tL == tL else nL)
            essE_P.append(nP / tE if tE == tE else nP)
            iatL_P.append(tP * L if tP == tP else 0.0)
        zE, zL = np.array(zE), np.array(zL)
        # windows, pooled over replicates
        SP, SL = sh[:, Pw].mean(1), sh[:, Lw].mean(1)                # (r, 7) replicate means
        EP, EL = eff[:, Pw].mean(1), eff[:, Lw].mean(1)
        leadP, leadL = int(SP.mean(0).argmax()), int(SL.mean(0).argmax())
        seW = math.sqrt(EP.var(ddof=1) / nr + EL.var(ddof=1) / nr) if nr > 1 else float("nan")
        tW = (EL.mean() - EP.mean()) / seW if seW > 0 else 0.0
        seX = math.sqrt(EP.var(ddof=1) / nr + Ep[i].var(ddof=1) / 10) if nr > 1 else float("nan")
        tX = (EP.mean() - Ep[i].mean()) / seX if seX > 0 else 0.0
        # friendly rivals: blocks with 111 above one half; takeovers (0 -> 1) and escapes (1 -> 0)
        st = sh[:, :, 6] > 0.5
        up = int((~st[:, :-1] & st[:, 1:]).sum())
        dn = int((st[:, :-1] & ~st[:, 1:]).sum())
        exp0 = float((~st[:, :-1]).sum() * L) / 1e7
        exp1 = float((st[:, :-1]).sum() * L) / 1e7
        ever = int(st.any(1).sum())

        def garwood(k, e):
            if e <= 0:
                return (float("nan"),) * 3
            lo = chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0
            return k / e, lo / e, chi2.ppf(0.975, 2 * k + 2) / 2 / e
        ru, rd = garwood(up, exp0), garwood(dn, exp1)
        ev_lo = betad.ppf(0.025, ever, nr - ever + 1) if ever > 0 else 0.0
        ev_hi = betad.ppf(0.975, ever + 1, nr - ever) if ever < nr else 1.0
        split = prod_split[i]
        ok_q = int((np.abs(zE) > 2).sum()) <= 2 and int((np.abs(zL) > 2).sum()) <= 2
        # [verifier fix] the triage's criterion is 'ESS at least 50 PER REPLICATE': every replicate, not the median
        ok_ess = float(np.min(essL_P)) >= 50
        ok_w = leadP == leadL and abs(tW) <= 2
        rows.append((i, G[i]["tag"], G[i]["wedge"], "split" if split else "agree", nr, ORDER[lead], ORDER[leadP],
                     ORDER[leadL], float(Ep[i].mean()), float(EP.mean()), float(EL.mean()), tW, tX,
                     float(np.median(essL_P)), float(np.min(essL_P)), float(np.median(essL_L)), float(np.median(essE_P)),
                     float(np.median(iatL_P)), float(np.abs(zE).max()), int((np.abs(zE) > 2).sum()),
                     float(np.abs(zL).max()), int((np.abs(zL) > 2).sum()), ever, ev_lo, ev_hi, up, ru[0], ru[1], ru[2],
                     dn, rd[0], rd[1], rd[2], ok_q, ok_ess, ok_w))
    if not rows:
        P("  no complete time series")
        return
    hdr = ("ipt tag wedge prod_reps nrep lead_prod lead_P lead_L E_prod E_P E_L t_P_vs_L t_P_vs_prod ESSlead_P_med "
           "ESSlead_P_min ESSlead_L_med ESSeff_P_med IATlead_P_gen zmax_eff n_z2_eff zmax_lead n_z2_lead "
           "reps_111_above_half ci_lo ci_hi takeovers rate_per1e7 lo hi escapes rate_per1e7 lo hi "
           "ok_quarters ok_ess ok_windows").split()
    with open(os.path.join(resd, "%s_games.tsv" % arm), "w") as fo:
        fo.write("# TS diagnostics per game (memory one, F1, eps 1e-4, 4e7 generations): P = production window "
                 "(itend/8, itend/4], L = (itend/2, itend]; lead share = the production pooled leader's share\n")
        fo.write("\t".join(hdr) + "\n")
        for r in rows:
            fo.write("\t".join("%.6g" % x if isinstance(x, float) else str(x) for x in r) + "\n")
    np.savez_compressed(os.path.join(resd, "%s_series.npz" % arm),
                        **{"eff_%d" % i: v[0] for i, v in series.items()}, **{"share_%d" % i: v[1] for i, v in series.items()},
                        order=np.array(ORDER))
    P("  %-4s %-11s %-5s %4s %4s %4s %4s  %7s %7s %7s %6s  %7s %7s %9s  %5s %5s  %-9s %s"
      % ("ipt", "tag", "prod", "lead", "P", "L", "nrep", "E_prod", "E_P", "E_L", "t(PL)", "ESS_P", "ESS_L", "IAT(gen)",
         "z>2E", "z>2L", "111>1/2", "verdict"))
    for r in rows:
        (i, tag, w, sp, nr, lp, lP, lL, E0, EPm, ELm, tW, tX, eP, ePmin, eL, eEP, iat, zE, nzE, zL, nzL, ev, elo, ehi,
         up, *_rest) = r
        okq, oke, okw = r[-3], r[-2], r[-1]
        verdict = ("PASS" if (okq and oke and okw) else "FAIL") if sp == "agree" else "split: report"
        P("  %-4d %-11s %-5s %4s %4s %4s %4d  %7.4f %7.4f %7.4f %6.2f  %7.0f %7.0f %9.3g  %5d %5d  %2d/%-2d    %s"
          % (i, tag, sp, lp, lP, lL, nr, E0, EPm, ELm, tW, eP, eL, iat, nzE, nzL, ev, nr, verdict))
    na = [r for r in rows if r[3] == "agree"]
    P("  games whose production replicates agree: %d; all three TS criteria met at %d of them"
      % (len(na), sum(r[-3] and r[-2] and r[-1] for r in na)))
    P("  per-game table: %s; the block series: %s" % (os.path.relpath(os.path.join(resd, "%s_games.tsv" % arm), HERE),
                                                    os.path.relpath(os.path.join(resd, "%s_series.npz" % arm), HERE)))


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["--kappa-check"]:
        sys.exit(kappa_check(a[1:]))
    if a[:1] == ["--probe-times"]:
        sys.exit(probe_times(a[1], a[2]))
    local = "--local" in a
    arms = [x for x in a if not x.startswith("--")] or ["E3", "E5", "TS"]
    analysis(arms, local=local)
