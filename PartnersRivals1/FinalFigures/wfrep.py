#!/usr/bin/env python3
"""
wfrep.py -- the Wright-Fisher runs of the referee round (Cannon kits eh, el, c1..c4; Cannon/mkkits.py), read into the frame
of Figures 4 and 5, plus their per-replicate atom shares.

    load("eh_e3")   -> dict with
        S      (512, 7)      pooled share of the population in each atom, ORDER = 000 100 010 001 110 011 111
        N      (512, 7)      strategies of the space in each atom at the game (exact eps -> 0 atoms, atoms512.bin)
        E      (512,)        efficiency: mean per-round payoff over the sampled half / Emax
        UV     (512, 2)      the game
        R      (512, 10, 7)  the share of each REPLICATE's sampled population in each atom (the .rep file)
        ER     (512, 10)     each replicate's efficiency
        LEAD   (512, 10)     the atom (index into ORDER) of each replicate's most abundant strategy (the .win file)
        nrep   (512,)        replicates found (10 when complete)
        D      dict(N=..., beta=..., mu="10^{-2}", muplain="1e-2", eps="10^{-3}", epsplain="1e-3", label=...)
    Cells in ipt order 0..511 (gid 2000 + ipt), like wfdata.load.

The packs are <kit>_<e>.{idx,bin,win,rep} in ../Data/ (written on Cannon by each kit's pack.py; the .rep is new in this
round).  atoms512.bin (Cannon/atoms512.bin) gives the atom of every strategy at every cell, so the pooled shares S are
computed here from the .bin exactly as pack.py computes the replicate shares, and the two are checked against each other.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "Data")
ATOMS = os.path.join(HERE, "..", "Cannon", "atoms512.bin")
ORDER = ["000", "100", "010", "001", "110", "011", "111"]
IDX = [int(c, 2) for c in ORDER]                       # atom code 4*eff + 2*nash + riv, as in atoms512.bin
NC = 65536

# the kits and what they are
RUNS = {
    "eh_e5": dict(N=1000, beta=100, mu="10^{-2}", muplain="1e-2", eps="10^{-5}", epsplain="1e-5"),
    "eh_e3": dict(N=1000, beta=100, mu="10^{-2}", muplain="1e-2", eps="10^{-3}", epsplain="1e-3"),
    "el_e5": dict(N=100, beta=3, mu="10^{-4}", muplain="1e-4", eps="10^{-5}", epsplain="1e-5"),
    "el_e3": dict(N=100, beta=3, mu="10^{-4}", muplain="1e-4", eps="10^{-3}", epsplain="1e-3"),
    "c1_e4": dict(N=100, beta=100, mu="10^{-2}", muplain="1e-2", eps="10^{-4}", epsplain="1e-4"),
    "c2_e4": dict(N=100, beta=3, mu="10^{-2}", muplain="1e-2", eps="10^{-4}", epsplain="1e-4"),
    "c3_e4": dict(N=100, beta=100, mu="10^{-4}", muplain="1e-4", eps="10^{-4}", epsplain="1e-4"),
    "c4_e4": dict(N=1000, beta=3, mu="10^{-4}", muplain="1e-4", eps="10^{-4}", epsplain="1e-4"),
    # the seeded re-runs of Figure 4's two runs (kits f1, f2: dw1's and dw2's parameters and seeds, plus the r-file)
    "f1_e4": dict(N=1000, beta=100, mu="10^{-2}", muplain="1e-2", eps="10^{-4}", epsplain="1e-4"),
    "f2_e4": dict(N=100, beta=3, mu="10^{-4}", muplain="1e-4", eps="10^{-4}", epsplain="1e-4"),
}
for k, d in RUNS.items():
    d["label"] = "Memory-2, $N = %d$\n$\\beta = %d$, $\\mu = %s$, $\\epsilon = %s$" % (d["N"], d["beta"], d["mu"], d["eps"])


def atoms():
    a = np.fromfile(ATOMS, dtype=np.uint8)
    assert a.size == 512 * NC, ATOMS
    return a.reshape(512, NC)


def load(run, data=DATA, strict=True):
    D = dict(RUNS[run])
    base = os.path.join(data, run)
    rows = [l.split() for l in open(base + ".idx") if not l.startswith("#")]
    blob = np.fromfile(base + ".bin", dtype="<f4")
    assert blob.size == len(rows) * NC, (base, blob.size, len(rows))
    blob = blob.reshape(len(rows), NC).astype(np.float64)
    A = atoms()
    ncell = 512 if strict else len(rows)
    S = np.full((512, 7), np.nan); N = np.zeros((512, 7), dtype=np.int64); E = np.full(512, np.nan)
    UV = np.full((512, 2), np.nan); nrep = np.zeros(512, dtype=np.int64)
    R = np.full((512, 10, 7), np.nan); ER = np.full((512, 10), np.nan); LEAD = np.full((512, 10), -1, dtype=np.int64)
    seen = []
    for k, r in enumerate(rows):
        gid, tag = int(r[0]), r[2]
        ipt = int(tag.split("_")[1][1:]); assert gid == 2000 + ipt, (gid, tag)
        cR, cS, cT = map(float, r[5:8]); assert cR == 1.0
        UV[ipt] = (cS, cT - 1.0)
        E[ipt] = float(r[10]); nrep[ipt] = int(r[13])
        p = blob[k] / blob[k].sum()
        full = np.bincount(A[ipt], weights=p, minlength=8)
        assert full[5] == 0.0
        S[ipt] = full[IDX]
        N[ipt] = np.bincount(A[ipt], minlength=8)[IDX]
        seen.append(ipt)
    assert N.sum(1)[seen].tolist() == [NC] * len(seen)
    for l in open(base + ".rep"):
        if l.startswith("#"):
            continue
        f = l.split(); ipt = int(f[0]) - 2000; irun = int(f[1]) - 1
        ER[ipt, irun] = float(f[2]); R[ipt, irun] = [float(x) for x in f[3:10]]
    for l in open(base + ".win"):
        if l.startswith("#"):
            continue
        f = l.split(); ipt = int(f[0]) - 2000; irun = int(f[1]) - 1
        LEAD[ipt, irun] = IDX.index(int(A[ipt, int(f[2])]))
    if strict:
        assert sorted(seen) == list(range(512)), "%s: %d cells" % (run, len(seen))
        assert (nrep == 10).all(), "cells short on replicates: %s" % np.nonzero(nrep != 10)[0]
    # the replicate shares average to the pooled shares (same sampled draws per replicate)
    for ipt in seen:
        if nrep[ipt] == 10:
            assert np.allclose(R[ipt].mean(0), S[ipt], atol=2e-6), (run, ipt, R[ipt].mean(0) - S[ipt])
    return dict(S=S, N=N, E=E, UV=UV, R=R, ER=ER, LEAD=LEAD, nrep=nrep, D=D, cells=np.array(sorted(seen)))


def summary(g):
    """one line of the numbers the text quotes: efficiency below 0.9, winners by atom, replicate agreement"""
    c = g["cells"]; S, E, R, LEAD = g["S"][c], g["E"][c], g["R"][c], g["LEAD"][c]
    win = S.argmax(1)
    agree_lead = (LEAD == win[:, None]).sum(1)                       # replicates whose leading strategy is in the pooled winner atom
    agree_share = (R.argmax(2) == win[:, None]).sum(1)               # replicates whose own dominant atom is the pooled winner
    out = ["%d cells; efficiency < 0.9 at %d, < 0.5 at %d, negative at %d" % (len(c), (E < 0.9).sum(), (E < 0.5).sum(), (E < 0).sum())]
    out.append("most abundant: " + ", ".join("%s %d" % (ORDER[i], (win == i).sum()) for i in range(7) if (win == i).any()))
    out.append("replicates agreeing with the pooled dominant atom: median %d/10 by their own shares, %d/10 by their leading strategy; all ten at %d cells"
               % (np.median(agree_share), np.median(agree_lead), (agree_share == 10).sum()))
    out.append("efficiency s.d. over replicates: median %.3f, max %.3f" % (np.median(np.nanstd(g["ER"][c], axis=1, ddof=1)), np.nanstd(g["ER"][c], axis=1, ddof=1).max()))
    return "\n".join(out)


if __name__ == "__main__":
    import sys
    for run in sys.argv[1:] or [r for r in RUNS if os.path.exists(os.path.join(DATA, r + ".idx"))]:
        g = load(run, strict=False)
        print("== %s: %s" % (run, RUNS[run]["label"].replace("\n", " ")))
        print(summary(g))
