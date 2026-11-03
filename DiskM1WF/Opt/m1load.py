#!/usr/bin/env python3
"""
The data layer of this project's figures: one packed WF run plus the exact atom classification, per
game, in the shape PartnersRivals/FinalFigures expects.

    python3 m1load.py            print both runs' per-atom summary and run every check

`collect(which)` is the memory-one counterpart of DiskM2WF/Opt/figF4.collect and returns the same
four arrays, so FinalFigures/atomshares.py's machinery can draw this project's runs unchanged:

    pts   (512,)    the sunflower cell index of each game, DiskM2WF/Opt/vor's ipt
    S     (512, 8)  the time-averaged share of the population in each atom
    N     (512, 8)  how many of the 16 strategies are in each atom at that game
    UV    (512, 2)  the game

TWO THINGS ARE CALLED F-SOMETHING AND THEY ARE NOT THE SAME.  The RUNS keep DiskM2WF's names,
F1 and F2 -- that is what `which` means everywhere in this module.  The four PAGES were renamed
F1..F4 on 2026-09-18.  Each run draws two pages:

    run F1   N = 1000, beta = 100, mu = 1e-2, eps = 1e-4   Data/mw1_e4.tsv   -> pages F1, F2
    run F2   N = 100,  beta = 3,   mu = 1e-4, eps = 1e-4   Data/mw2_e4.tsv   -> pages F3, F4
    run F1s  run F1 with itend cut 100x                    Data/mw1s_e4.tsv  -> converge.py only
    run F2s  run F2 with itend cut 100x                    Data/mw2s_e4.tsv  -> converge.py only

Figures/parameters.tex tabulates every parameter behind the four pages.

The two premises of the seven-panel figures are checked here rather than assumed, cell by cell, as
figF4.collect checks them for memory two: that the eight shares are a partition of the abundance
(they sum to 1), and that the 101 atom -- efficient and competitive but not nash -- is empty.  Also
checked: that the packed file's own u, v agree with the classification's, so a run cannot be paired
with the wrong game list.

Set MW_DATA to a directory holding mw1_e4.tsv / mw2_e4.tsv to draw a different pack (Data/dev holds
the short local run used while the cluster campaign was outstanding; the figures print which pack
they drew, and a dev pack refuses to be drawn without MW_DATA set, so a development page cannot be
mistaken for a final one).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("MW_DATA") or os.path.join(HERE, os.pardir, "Data")
NSTRAT = 16
K = 4.0                        # the projection of Figure 2; asserted by atomshares.py

DATASETS = {
    "F1": dict(stem="mw1", pack="mw1_e4.tsv", kit="mw", N=1000, beta=100,
               mu="10^{-2}", muplain="1e-2", itend="10^7", itplain="1e7", itn="10000000",
               nrep=10, subs="10^5"),
    "F2": dict(stem="mw2", pack="mw2_e4.tsv", kit="mw", N=100, beta=3,
               mu="10^{-4}", muplain="1e-4", itend="10^8", itplain="1e8", itn="100000000",
               nrep=10, subs="10^4"),
    # The convergence controls (Cannon/mw prod1s / prod2s): the same runs with
    # itend cut 100x and an independent seed base.  They are NOT deliverables --
    # converge.py compares them against F1 / F2 to say what the run length buys.
    "F1s": dict(stem="mw1s", pack="mw1s_e4.tsv", kit="mw", N=1000, beta=100,
                mu="10^{-2}", muplain="1e-2", itend="10^5", itplain="1e5", itn="100000",
                nrep=10, subs="10^3"),
    "F2s": dict(stem="mw2s", pack="mw2s_e4.tsv", kit="mw", N=100, beta=3,
                mu="10^{-4}", muplain="1e-4", itend="10^6", itplain="1e6", itn="1000000",
                nrep=10, subs="10^2"),
}
CONTROL = {"F1": "F1s", "F2": "F2s"}                   # deliverable -> its control


def packpath(which):
    return os.path.join(DATA, DATASETS[which]["pack"])


def read_pack(which):
    """the packed run: per cell (ipt, gid, tag, u, v, ntot, counts[16]) and the header's own run line"""
    p = packpath(which)
    rows, hdr = [], []
    for line in open(p):
        if line.startswith("#"):
            hdr.append(line.rstrip())
            continue
        f = line.split("\t")
        rows.append((int(f[0]), int(f[1]), f[2], float(f[3]), float(f[4]), float(f[5]),
                     float(f[6]), float(f[7]), int(f[8]), [int(x) for x in f[9:9 + NSTRAT]]))
    assert len(rows) == 512, "%s has %d cells, expected 512" % (p, len(rows))
    assert [r[0] for r in rows] == list(range(512)), "%s is not in ipt order" % p
    run = [h for h in hdr if h.startswith("# N =")]
    assert len(run) == 1, "%s has no single run line" % p
    return rows, run[0]


def collect(which="F1"):
    """per cell: the eight atom shares and the eight set sizes (figF4.collect's signature)"""
    sys.path.insert(0, HERE)
    import m1disk                                                     # noqa: E402
    d = m1disk.load()
    rows, run = read_pack(which)
    D = DATASETS[which]
    want = "# N = %d, beta = %d, u = %s," % (D["N"], D["beta"], D["muplain"])
    assert run.startswith(want), "%s says %r, not %r" % (packpath(which), run, want + " ...")
    dev = ("itend = %s," % D["itn"]) not in run
    if dev and not os.environ.get("MW_DATA"):
        raise ValueError("%s does not have %s's run length (%s) -- set MW_DATA to draw it "
                         "deliberately" % (packpath(which), which, run))
    pts, S, N, UV = [], [], [], []
    codes, atoms, uv, tags = d["codes"], d["atoms"], d["uv"], d["tags"]
    for ipt, gid, tag, u, v, emax, pay, ef, ntot, cnt in rows:
        assert tag == str(tags[ipt]), "cell %d: pack says %s, classification says %s" % (ipt, tag, tags[ipt])
        assert (u, v) == tuple(uv[ipt]), "cell %d (%s): pack (u, v) = (%r, %r), classification (%r, %r)" \
            % (ipt, tag, u, v, uv[ipt][0], uv[ipt][1])
        c = codes[ipt]
        s = np.zeros(8)
        for k in range(NSTRAT):
            s[c[k]] += cnt[k]
        s /= float(ntot)
        n = atoms[ipt]
        if n[5]:
            raise ValueError("%s (u=%.4f, v=%.4f): the 101 atom is NOT empty -- %d strategies are "
                             "efficient and competitive but not nash.  The figure's premise fails; "
                             "do not draw." % (tag, u, v, n[5]))
        if abs(s.sum() - 1.0) > 1e-12:
            raise ValueError("%s: the eight shares sum to %.14f, not 1" % (tag, s.sum()))
        if s[5] != 0.0:
            raise ValueError("%s: the empty 101 atom holds %.3e of the population" % (tag, s[5]))
        pts.append(ipt); S.append(s); N.append(np.asarray(n, float)); UV.append((u, v))
    return np.array(pts), np.array(S), np.array(N), np.array(UV)


def summary(which):
    import m1disk
    pts, S, N, UV = collect(which)
    D = DATASETS[which]
    rows, run = read_pack(which)
    print("%s  (%s)" % (which, run[2:]))
    print("  %-5s %-34s %8s %8s %8s %8s %7s" % ("atom", "", "mean", "median", "min", "max", "empty"))
    for c in m1disk.ORDER:
        i = int(c, 2)
        v, gap = S[:, i], N[:, i] == 0
        print("  %-5s %-34s %8.4f %8.4f %8.5f %8.4f %7d"
              % (c, m1disk.NAME[c], v.mean(), np.median(v), v[~gap].min() if (~gap).any() else np.nan,
                 v.max(), int(gap.sum())))
    r = np.where(N[:, [int(c, 2) for c in m1disk.ORDER]] == 0, np.nan,
                 S[:, [int(c, 2) for c in m1disk.ORDER]] * NSTRAT
                 / np.maximum(N[:, [int(c, 2) for c in m1disk.ORDER]], 1))
    with np.errstate(divide="ignore"):
        lr = np.log10(np.maximum(r, 1e-300))
    print("  enrichment r = share / (|atom| / 16):  min %.3g, max %.3g   (the ceiling is 16/|atom| <= 16)"
          % (np.nanmin(r), np.nanmax(r)))
    print("  log10 r in [%.2f, %.2f];  cells outside 10^-2..10^2: %d of %d defined"
          % (np.nanmin(lr), np.nanmax(lr), int((np.abs(lr) > 2).sum()), int(np.isfinite(lr).sum())))
    win = np.nanargmax(S[:, [int(c, 2) for c in m1disk.ORDER]], 1)
    wr = np.nanargmax(np.where(np.isnan(lr), -np.inf, lr), 1)
    print("  most abundant atom: " + ", ".join("%s %d" % (m1disk.ORDER[i], int((win == i).sum()))
                                               for i in range(7) if (win == i).any()))
    print("  most enriched atom: " + ", ".join("%s %d" % (m1disk.ORDER[i], int((wr == i).sum()))
                                               for i in range(7) if (wr == i).any()))
    return S, N


if __name__ == "__main__":
    for w in ("F1", "F2"):
        if os.path.exists(packpath(w)):
            summary(w)
        else:
            print("%s: %s not present yet" % (w, packpath(w)))
