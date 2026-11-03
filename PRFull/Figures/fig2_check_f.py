#!/usr/bin/env python3
"""
fig2_check_f.py -- a receipt for panel f of Figure 2: the friendly rivals of wedge S drawn from the formula
FR_S = {p1 = 1, p2 > 0, p3 > 0, p2 + p3 >= 1, (p4 > 0 or max(p2, p3) = 1)}, checked piece by piece against the RAW
limit tests of PartnersRivalsM1f/check/m1lib.py (the classifiers of fig2_check.py): competitive where T < S by the
rival margin against the 16 co-players, efficient by the efficiency deficit, both at two error rates.

    python3 fig2_check_f.py      ->  prints, and writes cache/fig2_check_f.txt

Points are drawn uniformly in the relative interior of every piece of the prism Riv- (the interior, its five faces,
its nine edges), at a wedge-S game (2.5, -4) and at a second one (0.5, -2), and the raw verdict is compared with the
formula.  Read-only imports; no bytecode is written.
"""
import os
import sys

sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M1F = "/Users/martin/Documents/CL1/PartnersRivalsM1f"
sys.path.insert(0, os.path.join(M1F, "check"))
import m1lib as L                                                  # noqa: E402

RNG = np.random.default_rng(20260927)
NPT = 25
lines = []


def say(s):
    print(s)
    lines.append(s)


def efficient(p, u, v):
    d1 = L.efficiency_deficit(p, u, v, 1e-7)
    d2 = L.efficiency_deficit(p, u, v, 1e-9)
    return d2 > -1e-5 or d2 > 0.2 * d1


def rival(p, sign):
    m1 = L.rival_margin(p, sign, 1e-7)
    m2 = L.rival_margin(p, sign, 1e-9)
    return m2 > -1e-6 or m2 > 0.2 * m1


def formula(p):
    _, p2, p3, p4 = p
    return p2 > 0 and p3 > 0 and p2 + p3 >= 1 and (p4 > 0 or max(p2, p3) == 1)


def r(n=1):
    return RNG.uniform(0.05, 0.95, n)


def tri():                                         # (p2, p3) in the open triangle p2 + p3 > 1, p2, p3 < 1
    while True:
        a, b = r(2)
        if a + b > 1.02:
            return a, b


PIECES = {                                         # name: sampler of (p2, p3, p4), p1 = 1
    "interior":                         lambda: (*tri(), r()[0]),
    "face p4 = 0 (open)":               lambda: (*tri(), 0.0),
    "face p4 = 1 (open)":               lambda: (*tri(), 1.0),
    "face p2 = 1 (open)":               lambda: (1.0, r()[0], r()[0]),
    "face p3 = 1 (open)":               lambda: (r()[0], 1.0, r()[0]),
    "face p2 + p3 = 1 (open)":          lambda: (lambda a: (a, 1 - a, r()[0]))(r()[0]),
    "edge p4 = 0, p2 = 1":              lambda: (1.0, r()[0], 0.0),
    "edge p4 = 0, p3 = 1":              lambda: (r()[0], 1.0, 0.0),
    "edge p4 = 0, p2 + p3 = 1":         lambda: (lambda a: (a, 1 - a, 0.0))(r()[0]),
    "edge p4 = 1, p2 = 1":              lambda: (1.0, r()[0], 1.0),
    "edge p4 = 1, p3 = 1":              lambda: (r()[0], 1.0, 1.0),
    "edge p4 = 1, p2 + p3 = 1":         lambda: (lambda a: (a, 1 - a, 1.0))(r()[0]),
    "edge (p2, p3) = (1, 1)":           lambda: (1.0, 1.0, r()[0]),
    "edge (p2, p3) = (0, 1), TFT-1011": lambda: (0.0, 1.0, r()[0]),
    "edge (p2, p3) = (1, 0), deleted":  lambda: (1.0, 0.0, r()[0]),
}
bad = 0
for (u, v) in ((2.5, -4.0), (0.5, -2.0)):
    assert u + v < 1 and v < u - 1, "not a wedge-S game"
    say("wedge-S game (u, v) = (%.1f, %.1f): piece, formula, raw friendly rival (competitive and efficient)" % (u, v))
    for name, f in PIECES.items():
        verd = []
        for _ in range(NPT):
            p = np.array((1.0,) + tuple(float(x) for x in f()))
            raw = rival(p, -1) and efficient(p, u, v)
            verd.append((formula(p), raw))
        fs = {a for a, _ in verd}
        agree = sum(a == b for a, b in verd)
        bad += NPT - agree
        say("  %-34s formula %-11s raw agrees at %d/%d" % (name, "/".join(sorted(map(str, fs))), agree, NPT))
say("disagreements: %d" % bad)
os.makedirs(os.path.join(HERE, "cache"), exist_ok=True)
with open(os.path.join(HERE, "cache", "fig2_check_f.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
