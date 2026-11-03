#!/usr/bin/env python3
"""
wfdata.py -- the four Wright-Fisher runs behind Figures 4 and 5 and SI Figures 2 and 3 (fig4.py, fig5.py, sifig2.py, sifig3.py;
also other/figwfprop.py), read the same way:

    load("m2", "F1")   memory two, N = 1000, beta = 100, mu = 1e-2   (DiskM2WF run F1: Figures 3, 4)
    load("m2", "F2")   memory two, N = 100,  beta = 3,   mu = 1e-4   (DiskM2WF run F2: SI Figures 2, 3)
    load("m1", "F1")   memory one, N = 1000, the same parameters     (DiskM1WF run F1: SI Figures 6, 7)
    load("m1", "F2")   memory one, N = 100                            (DiskM1WF run F2: SI Figures 8, 9)

returns, per cell in ipt order 0..511 (the Voronoi cells of the Fermat sunflower, vor.cells()):
    S     (512, 7)  the share of the population in each atom, in the order ORDER = 000 100 010 001 110 011 111
    N     (512, 7)  how many strategies of the space are in the atom at that game (exact eps -> 0 classification)
    E     (512,)    the efficiency: mean per-round payoff over the second half of the run / Emax
    UV    (512, 2)  the game
    D     the run's parameter row (N, beta, muplain, ...)
Memory two through DiskM2WF/Opt (figF4.collect, wfload), memory one through DiskM1WF/Opt/m1load; all read-only.
"""
import sys

import numpy as np

M1OPT = "/Users/martin/Documents/CL1/DiskM1WF/Opt"
if M1OPT not in sys.path:
    sys.path.insert(0, M1OPT)
import atomshares as A                                             # noqa: E402  puts DiskM2WF/Opt on sys.path
import figF4, finalpanels, vor, wfload                             # noqa: E402
import m1load                                                      # noqa: E402

ORDER = A.ORDER
IDX = [int(c, 2) for c in ORDER]
NSTRAT = {"m1": 16.0, "m2": 65536.0}
LABEL = {("m1", "F2"): "Memory-1, $N = 100$", ("m1", "F1"): "Memory-1, $N = 1000$",
         ("m2", "F2"): "Memory-2, $N = 100$", ("m2", "F1"): "Memory-2, $N = 1000$"}
ROWS = [("m1", "F2"), ("m1", "F1"), ("m2", "F2"), ("m2", "F1")]   # the rows of Figures 4 and 5, top to bottom


def load(space, which):
    if space == "m2":
        pts, S, N, UV = figF4.collect(which)
        G = wfload.load_run(finalpanels.idxpath(which), ie=3)
        eff = {vor.index_of(g["tag"]): float(g["eff"]) for g in G.values() if vor.index_of(g["tag"]) is not None}
        D = finalpanels.DATASETS[which]
    else:
        pts, S, N, UV = m1load.collect(which)
        rows, _ = m1load.read_pack(which)
        eff = {r[0]: float(r[7]) for r in rows}
        D = m1load.DATASETS[which]
    assert sorted(pts.tolist()) == list(range(512)) and sorted(eff) == list(range(512)), (space, which)
    assert np.allclose(S.sum(1), 1.0) and not S[:, 5].any() and not N[:, 5].any()
    assert (N.sum(1) == NSTRAT[space]).all()
    o = np.argsort(pts)
    return S[o][:, IDX], N[o][:, IDX], np.array([eff[i] for i in range(512)]), UV[o], D


def cells():
    """the 512 Voronoi polygons and the indices of those with at least three vertices"""
    c = vor.cells()
    keep = [i for i in range(512) if len(c[i]) >= 3]
    return [c[i] for i in keep], keep
