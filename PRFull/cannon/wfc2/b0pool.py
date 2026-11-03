#!/usr/bin/env python3
"""b0pool.py -- the beta = 0 neutral check with 20 replicates per game: B0 (array 48926421, seed base 34000000) and B0b
(array 48942211, seed base 34200000) pooled for the global kernel, B0K0 and B0K0b (34500000, 34700000) for the local
kernel K0.  Same test as analyze.py report_neutral; the B0b/B0K0b replicates are renumbered 11..20.  Local, seconds."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze as A

def pooled(a, b):
    d = os.path.join(HERE, "out")
    R = A.load_arm(os.path.join(d, a), "1o3"); Rb = A.load_arm(os.path.join(d, b), "1o3")
    assert sorted(R) == sorted(Rb)
    for i in R:
        for r in Rb[i]:
            r = dict(r); r["irep"] += 10; R[i].append(r)
        R[i].sort(key=lambda r: r["irep"])
        assert len({r["irep"] for r in R[i]}) == len(R[i]) == 20
    ipts = sorted(R)
    M, H, MD = A.load_mass([A.PROD["F1"]], "1o3", ipts)
    G = {i: A.game_stats(R[i], M.get(i), H.get(i), MD.get(i)) for i in ipts}
    A.say(""); A.say("=" * 118)
    A.say("%s + %s pooled: %d games, 20 replicates each" % (a, b, len(ipts)))
    A.report_neutral(a + "+" + b, R, G, M, MD, ipts)
    # each half alone, for comparison
    for arm in (a, b):
        Ra = A.load_arm(os.path.join(d, arm), "1o3")
        Ga = {i: A.game_stats(Ra[i], M.get(i), H.get(i), MD.get(i)) for i in ipts}
        A.say("--- %s alone" % arm)
        A.report_neutral(arm, Ra, Ga, M, MD, ipts)

pooled("B0", "B0b")
pooled("B0K0", "B0K0b")
open(os.path.join(HERE, "analysis", "b0pool.txt"), "w").write("\n".join(A.OUT) + "\n")
