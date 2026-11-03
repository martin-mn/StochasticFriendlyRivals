#!/usr/bin/env python3
"""
cwdata.py -- the CONTINUOUS memory-one Wright-Fisher run (kit ../cannon/cw) read for fig4c.py and fig5c.py,
the continuum counterparts of ../../PartnersRivals1/FinalFigures/fig4.py and fig5.py.

    load()                  the real pack:      ../cannon/cw/packed/cw1_e4.tsv   + ../cannon/cw/packed/kappa_e4.tsv
    load(synthetic=True)    fabricated stand-in: synthetic/cw1_e4.SYNTHETIC.tsv + synthetic/kappa_e4.SYNTHETIC.tsv
                            (written by make_synthetic(); NEVER into ../cannon/cw/packed; a file whose first line
                            says SYNTHETIC is refused by load() unless synthetic=True, and vice versa)
    load(short=True)        a SHORT LOCAL RUN of the real chain: ../cannon/cw/packed_local/ (pack.py --dev), whose
                            first line is '# SHORT LOCAL RUN ...'; real data from a run far too short to be a result.
                            Refused unless short=True, and vice versa; any itend and any M are accepted.
    load(..., run="F2")     the same three kinds for run F2, the N = 100 run of SI Figures 4 and 5 (sf2.f's parameters:
                            N = 100, beta = 3, u = 1e-4, eps = 1e-4, itend = 1e8, 10 replicates; kit program cw2.x):
                            real ../cannon/cw/packed/cw2_e4.tsv, short ../cannon/cw/packed_local2/cw2_e4.tsv,
                            synthetic synthetic/cw2_e4.SYNTHETIC.tsv.  The KERNEL FILE IS F1's -- kappa_e4.tsv of
                            packed/, packed_local/, synthetic/ respectively: the atoms of a strategy depend only on the
                            game and the kernel is the same exact-alpha kernel, so kappa does not depend on N, beta or u
                            (see KAPPA_F2).  run="F1" is the default, and every F1 call is exactly what it was.

returns, per cell in ipt order 0..511 (the Voronoi cells of wfdata.cells(), the same games as the binary runs):
    S     (512, 7)  share of the population in each atom, in ORDER = 000 100 010 001 110 011 111
    K     (512, 7)  the mutation kernel's mass on each atom at that game (kappa; Monte Carlo, M draws per game)
    E     (512,)    efficiency ef = pay / Emax
    UV    (512, 2)  the game
    plus  emax, pay, meanlive (512,), Z (512, 81) the kernel's stratum masses, M, D (the run line), synthetic.

THE PACK FORMATS (defined here; the kit's packer must write exactly this).  Tab-separated, one row per cell,
ipt order 0..511, '#' lines are header.  The column-header line is REQUIRED -- the reader maps columns by name
and refuses a file without it, so the atom order is never a matter of convention:

  cw1_e4.tsv
    # <free description lines>
    # N = 1000, beta = 100, u = 1e-2, eps = 1e-4, itend = 10000000, 10 replicates, iseedb = 40000000
    # ipt  gid  u  v  emax  pay  ef  s000  s100  s010  s001  s110  s011  s111  s101  meanlive
    - u, v      exactly as games.dat (%.17g); checked against the binary run's (u, v) at the same ipt
    - emax      max(R, (S+T)/2) = max(1, (1+u+v)/2)
    - pay, ef   mean over the 10 replicates of sf1's per-replicate pay and ef; ef = pay/emax is asserted
    - s<code>   the 8 atom shares, the 7 of ORDER then the empty 101: atomcnt(atom)/sum(atomcnt) summed over
                the 10 replicates, atom = 4E + 2N + R from m1cls; written %.17g; they must sum to 1 (1e-9)
                and s101 must be exactly 0 (the theorem: efficient + competitive => Nash)
    - meanlive  mean number of live lineage classes per sampled generation (diagnostic only)

  cw2_e4.tsv   (run F2) the same columns and rules; its run line is
    # N = 100, beta = 3, u = 1e-4, eps = 1e-4, itend = 100000000, 10 replicates, iseedb = 42000000
    (the reader checks the line up to and including 'itend = 100000000,'; in short mode up to 'itend = ')

  kappa_e4.tsv
    # <free description lines>
    # M = 10000000, seed = 41000000 + isl
    # ipt  k000  k100  k010  k001  k110  k011  k111  k101  z00  z01 ... z80
    - k<code>   fraction of the M kernel draws in each atom at that game (the same 8 codes, the same order)
    - z<nn>     fraction in stratum nn = sum_k code_k 3^k, code 0: p_k = 0, 1: p_k = 1, 2: interior (two digits)
    both groups sum to 1; k101 must be exactly 0; each z is 1/81 for the exact-alpha kernel (checked at 6 sigma).
"""
import contextlib
import os
import sys

import numpy as np

sys.dont_write_bytecode = True                                     # write no __pycache__ into the read-only trees imported below
HERE = os.path.dirname(os.path.abspath(__file__))
FINAL = "/Users/martin/Documents/CL1/PartnersRivals1/FinalFigures"  # fig4.py, fig5.py, wfdata.py: imported, READ-ONLY
if FINAL not in sys.path:
    sys.path.insert(0, FINAL)
PACKDIR = os.path.join(HERE, os.pardir, "cannon", "cw", "packed")
SYNDIR = os.path.join(HERE, "synthetic")
SHORTDIR = os.path.join(HERE, os.pardir, "cannon", "cw", "packed_local")
ORDER = ["000", "100", "010", "001", "110", "011", "111"]
ATOMS8 = ORDER + ["101"]
SCOLS = ["s" + c for c in ATOMS8]
KCOLS = ["k" + c for c in ATOMS8]
ZCOLS = ["z%02d" % i for i in range(81)]
PCOLS = ["ipt", "gid", "u", "v", "emax", "pay", "ef"] + SCOLS + ["meanlive"]
QCOLS = ["ipt"] + KCOLS + ZCOLS
RUNS = {"F1": dict(N=1000, beta=100, muplain="1e-2", eps="1e-4", itn="10000000", nrep=10),     # sf1.f's parameters
        "F2": dict(N=100, beta=3, muplain="1e-4", eps="1e-4", itn="100000000", nrep=10)}      # sf2.f's parameters
STEM = {"F1": "cw1_e4", "F2": "cw2_e4"}                            # the population pack of each run
SEEDB = {"F1": 40000000, "F2": 42000000}                           # the continuum seed base of cw.x, cw2.x
SHORTDIRS = {"F1": SHORTDIR, "F2": os.path.join(HERE, os.pardir, "cannon", "cw", "packed_local2")}
KAPPA_F2 = ("the F1 kernel file: kappa depends on the game and the kernel only, and the atoms of a strategy on the game only; "
            "N, beta and u enter neither")
MARK = "# SYNTHETIC"
MARK_SHORT = "# SHORT LOCAL RUN"


def runline(run="F1", short=False):
    """the prefix the run line of the population pack must start with"""
    R = RUNS[run]
    full = "# N = %d, beta = %d, u = %s, eps = %s, itend = %s," % (R["N"], R["beta"], R["muplain"], R["eps"], R["itn"])
    return full.split("itend = ")[0] + "itend = " if short else full


RUN = RUNS["F1"]                                                   # the F1 names, kept for every existing caller
RUNLINE = runline("F1")
RUNLINE_SHORT = runline("F1", short=True)


def paths(synthetic=False, short=False, run="F1"):
    """(population pack, kernel pack); the kernel pack is F1's for both runs (KAPPA_F2)"""
    st = STEM[run]
    if synthetic:
        return os.path.join(SYNDIR, st + ".SYNTHETIC.tsv"), os.path.join(SYNDIR, "kappa_e4.SYNTHETIC.tsv")
    if short:
        return os.path.join(SHORTDIRS[run], st + ".tsv"), os.path.join(SHORTDIR, "kappa_e4.tsv")
    return os.path.join(PACKDIR, st + ".tsv"), os.path.join(PACKDIR, "kappa_e4.tsv")


def _kind(hdr):
    if hdr and hdr[0].startswith(MARK):
        return "SYNTHETIC"
    if hdr and hdr[0].startswith(MARK_SHORT):
        return "SHORT LOCAL RUN"
    return "real"


def _read(path, cols, synthetic, short=False):
    hdr, rows, names = [], [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                hdr.append(line.rstrip("\n"))
                f = line[1:].split()
                if f and f[0] == "ipt":
                    names = f
                continue
            if line.strip():
                rows.append(line.rstrip("\n").split("\t"))
    kind = _kind(hdr)
    want = "SYNTHETIC" if synthetic else "SHORT LOCAL RUN" if short else "real"
    if kind != want:
        raise ValueError("%s is %s data but was asked for as %s" % (path, kind, want))
    if names is None:
        raise ValueError("%s has no '# ipt ...' column-header line; the reader maps columns by name" % path)
    miss = [c for c in cols if c not in names]
    if miss:
        raise ValueError("%s lacks the columns %s" % (path, miss))
    assert len(rows) == 512, "%s has %d cells, expected 512" % (path, len(rows))
    assert all(len(r) == len(names) for r in rows), "%s: a row has the wrong number of fields" % path
    a = np.array([[float(x) for x in r] for r in rows])
    col = {n: a[:, i] for i, n in enumerate(names)}
    assert (col["ipt"] == np.arange(512)).all(), "%s is not in ipt order 0..511" % path
    return col, hdr


def load(synthetic=False, uv_ref=None, short=False, run="F1"):
    """the continuous run; uv_ref (512, 2), if given, is the binary run's games, which must be the same"""
    assert not (synthetic and short)
    assert run in RUNS, run
    pp, kp = paths(synthetic, short, run)
    for p in (pp, kp):
        if not os.path.exists(p):
            raise FileNotFoundError("%s does not exist%s" % (p, " -- run cwdata.make_synthetic()" if synthetic else
                                                               " -- run the kit's local short run and pack.py --dev"
                                                               if short else
                                                               " -- run with --synthetic, or pack the cw run first"))
    P, ph = _read(pp, PCOLS, synthetic, short)
    Q, qh = _read(kp, QCOLS, synthetic, short)
    rline = [h for h in ph if h.startswith("# N =")]
    rl = runline(run, short)
    assert len(rline) == 1 and rline[0].startswith(rl), "%s: run line %r, expected %r..." % (pp, rline, rl)
    mline = [h for h in qh if h.startswith("# M =")]
    assert len(mline) == 1, "%s has no single '# M = ...' line" % kp
    M = int(float(mline[0].split("=")[1].split(",")[0]))

    S8 = np.stack([P[c] for c in SCOLS], 1)
    K8 = np.stack([Q[c] for c in KCOLS], 1)
    Z = np.stack([Q[c] for c in ZCOLS], 1)
    UV = np.stack([P["u"], P["v"]], 1)
    emax, pay, E = P["emax"], P["pay"], P["ef"]
    assert np.allclose(emax, np.maximum(1.0, 0.5 * (1.0 + UV[:, 0] + UV[:, 1])), rtol=1e-12, atol=0), "emax != max(R,(S+T)/2)"
    assert np.allclose(E, pay / emax, rtol=1e-9, atol=1e-12), "ef != pay/emax"
    assert np.allclose(S8.sum(1), 1.0, atol=1e-9), "atom shares do not sum to 1 (worst %.3e)" % abs(S8.sum(1) - 1).max()
    assert np.allclose(K8.sum(1), 1.0, atol=1e-9) and np.allclose(Z.sum(1), 1.0, atol=1e-9), "kappa does not sum to 1"
    assert (S8 >= 0).all() and (K8 >= 0).all() and (Z >= 0).all()
    assert (S8[:, 7] == 0).all(), "the 101 atom holds population at %d cells -- the theorem fails; do not draw" % (S8[:, 7] > 0).sum()
    assert (K8[:, 7] == 0).all(), "the kernel puts mass on 101 at %d cells -- the classifier is wrong" % (K8[:, 7] > 0).sum()
    sd = np.sqrt((1 / 81) * (80 / 81) / M)
    assert (abs(Z - 1 / 81) < 6 * sd).all(), "a stratum mass is %.1f sigma from 1/81" % (abs(Z - 1 / 81) / sd).max()
    if uv_ref is not None:
        assert np.allclose(UV, uv_ref, rtol=0, atol=1e-12), "the cw games are not the binary run's games"
    return dict(S=S8[:, :7], K=K8[:, :7], E=E, UV=UV, emax=emax, pay=pay, meanlive=P["meanlive"], Z=Z, M=M,
                gid=P["gid"].astype(int), D=dict(RUNS[run], run=rline[0], mark=ph[0] if short else None), synthetic=synthetic,
                short=short, which=run, kmark=qh[0] if short else None)


def kfloor(M):
    """where the estimate of kappa is 0 but the population holds the atom: the 95% upper bound of a zero count (rule of three)"""
    return 3.0 / M


def enrichment(S, K, M):
    """share / kappa per (cell, atom); -inf where kappa = 0 and the share is 0 (no strategy -- the original's N = 0);
    where kappa_hat = 0 but share > 0 the kernel mass is below the Monte Carlo resolution: kappa is set to its
    95% upper bound 3/M, so that entry is a conservative LOWER bound on the enrichment. Returns (enr, nfloored)."""
    Kx = np.where((K == 0) & (S > 0), kfloor(M), K)
    enr = np.where(Kx > 0, S / np.where(Kx > 0, Kx, 1.0), -np.inf)
    return enr, int(((K == 0) & (S > 0)).sum())


# ----------------------------------------------------------------------------------------------- the watermark
FOOT_SYN = "SYNTHETIC DATA (cwdata.make_synthetic) -- fabricated to test the layout, NOT a result"


def short_footer(run="F1"):
    """the footer for a short local run, from the pack's own first line"""
    with open(paths(short=True, run=run)[0]) as fh:
        first = fh.readline()
    return "Memory-1 continuous row: " + first.lstrip("# ").strip()


@contextlib.contextmanager
def watermark(on, text="SYNTHETIC", footer=FOOT_SYN, size=110):
    """while active, every pyplot figure is stamped with `text` just before it is saved (fig4.draw and fig5.draw
    build and save their own figures; this reaches in without touching their code)"""
    if not on:
        yield
        return
    import matplotlib.pyplot as plt
    orig = plt.figure

    def figure(*a, **k):
        f = orig(*a, **k)
        save = f.savefig

        def savefig(*aa, **kk):
            if not getattr(f, "_cw_wm", False):
                f.text(0.5, 0.5, text, rotation=28, fontsize=size, fontweight="bold", color="red", alpha=0.20,
                       ha="center", va="center", zorder=1000)
                f.text(0.5, 0.004, footer, ha="center", va="bottom", fontsize=11, color="red", zorder=1000)
                f._cw_wm = True
            return save(*aa, **kk)
        f.savefig = savefig
        return f
    plt.figure = figure
    try:
        yield
    finally:
        plt.figure = orig


# ----------------------------------------------------------------------------------------------- the synthetic stand-in
def make_synthetic(M=10_000_000, seed=20260924, run="F1"):
    """fabricate plausible files in exactly the pack formats above, seeded from the binary memory-one run F1:
    kappa ~ multinomial(M, Dirichlet(0.2 + 3 |a|/16)) with 101 = 0; shares = kappa tilted by the binary enrichment;
    efficiency = the binary one, shrunk and jittered.  Watermarked SYNTHETIC in the files and on the page.
    run="F2": only synthetic/cw2_e4.SYNTHETIC.tsv, from the binary run F2 and the SAME synthetic kernel file
    (written first, by the F1 call, if it is missing), as the real F2 pack shares F1's kappa_e4.tsv."""
    if run == "F2":
        return _make_synthetic_f2(seed=20260925 if seed == 20260924 else seed)
    assert run == "F1", run
    import wfdata                                                   # noqa: E402  (FinalFigures on sys.path)
    import m1load                                                   # noqa: E402
    rng = np.random.default_rng(seed)
    Sb, Nb, Eb, UVb, _ = wfdata.load("m1", "F1")
    rows, _ = m1load.read_pack("F1")
    gid = np.array([r[1] for r in rows])
    os.makedirs(SYNDIR, exist_ok=True)
    pp, kp = paths(True)
    K7 = np.zeros((512, 7)); S7 = np.zeros((512, 7))
    for i in range(512):
        p = rng.dirichlet(0.2 + 3.0 * Nb[i] / 16.0)
        K7[i] = rng.multinomial(M, p) / M
        tilt = np.where(Nb[i] > 0, Sb[i] / np.maximum(Nb[i] / 16.0, 1e-300), 0.3 * rng.lognormal(0, 1, 7))
        w = K7[i] * (tilt + 0.02) ** 0.8 * rng.lognormal(0, 0.3, 7)
        if w.sum() == 0:
            w = K7[i].copy()
        S7[i] = w / w.sum()
    Z = rng.multinomial(M, np.full(81, 1 / 81), size=512) / M
    u, v = UVb[:, 0], UVb[:, 1]
    emax = np.maximum(1.0, 0.5 * (1.0 + u + v))
    ef = np.clip(0.85 * Eb + 0.12 + rng.normal(0, 0.03, 512), -0.3, 1.0)
    pay = ef * emax
    meanlive = rng.uniform(40, 400, 512)
    desc = [MARK + " -- FABRICATED by PartnersRivalsM1f/Figures/cwdata.py make_synthetic(seed=%d); NOT DATA" % seed]
    with open(pp, "w") as fh:
        fh.write("\n".join(desc + ["# PartnersRivalsM1f cw1 (stand-in): the continuous memory-one WF run, exact-alpha kernel",
                                   RUNLINE + " 10 replicates, iseedb = 40000000",
                                   "# " + "\t".join(PCOLS)]) + "\n")
        for i in range(512):
            f = [str(i), str(gid[i]), repr(float(u[i])), repr(float(v[i])), repr(float(emax[i])), repr(float(pay[i])),
                 repr(float(pay[i] / emax[i]))] + [repr(float(x)) for x in S7[i]] + ["0.0", "%.3f" % meanlive[i]]
            fh.write("\t".join(f) + "\n")
    with open(kp, "w") as fh:
        fh.write("\n".join(desc + ["# PartnersRivalsM1f kappa (stand-in): kernel masses per atom and per stratum",
                                   "# M = %d, seed = 41000000 + isl" % M, "# " + "\t".join(QCOLS)]) + "\n")
        for i in range(512):
            fh.write("\t".join([str(i)] + [repr(float(x)) for x in K7[i]] + ["0.0"] + [repr(float(x)) for x in Z[i]]) + "\n")
    print("wrote %s and %s (SYNTHETIC)" % (os.path.relpath(pp, HERE), os.path.relpath(kp, HERE)))
    return pp, kp


def _make_synthetic_f2(seed):
    import wfdata                                                   # noqa: E402  (FinalFigures on sys.path)
    import m1load                                                   # noqa: E402
    pp, kp = paths(True, run="F2")
    if not os.path.exists(kp):
        make_synthetic()
    Q, _ = _read(kp, QCOLS, True)
    K7 = np.stack([Q["k" + c] for c in ORDER], 1)
    rng = np.random.default_rng(seed)
    Sb, Nb, Eb, UVb, _ = wfdata.load("m1", "F2")
    rows, _ = m1load.read_pack("F2")
    gid = np.array([r[1] for r in rows])
    S7 = np.zeros((512, 7))
    for i in range(512):
        tilt = np.where(Nb[i] > 0, Sb[i] / np.maximum(Nb[i] / 16.0, 1e-300), 0.3 * rng.lognormal(0, 1, 7))
        w = K7[i] * (tilt + 0.02) ** 0.8 * rng.lognormal(0, 0.3, 7)
        if w.sum() == 0:
            w = K7[i].copy()
        S7[i] = w / w.sum()
    u, v = UVb[:, 0], UVb[:, 1]
    emax = np.maximum(1.0, 0.5 * (1.0 + u + v))
    ef = np.clip(0.85 * Eb + 0.12 + rng.normal(0, 0.03, 512), -0.3, 1.0)
    pay = ef * emax
    meanlive = rng.uniform(1, 20, 512)
    desc = [MARK + " -- FABRICATED by PartnersRivalsM1f/Figures/cwdata.py make_synthetic(run='F2', seed=%d); NOT DATA" % seed]
    with open(pp, "w") as fh:
        fh.write("\n".join(desc + ["# PartnersRivalsM1f cw2 (stand-in): the continuous memory-one WF run F2, exact-alpha kernel",
                                   runline("F2") + " 10 replicates, iseedb = %d" % SEEDB["F2"],
                                   "# " + "\t".join(PCOLS)]) + "\n")
        for i in range(512):
            f = [str(i), str(gid[i]), repr(float(u[i])), repr(float(v[i])), repr(float(emax[i])), repr(float(pay[i])),
                 repr(float(pay[i] / emax[i]))] + [repr(float(x)) for x in S7[i]] + ["0.0", "%.3f" % meanlive[i]]
            fh.write("\t".join(f) + "\n")
    print("wrote %s (SYNTHETIC; kernel file %s, shared with F1)" % (os.path.relpath(pp, HERE), os.path.relpath(kp, HERE)))
    return pp, kp
