#!/usr/bin/env python3
"""wfc_probe.py <dir of r_*.txt>  --  the alpha = 1/2 validation of cannon/wfc/wfc.c against the binary run of Figure 4d-f.

At alpha = 1/2 the continuum process mutates uniformly over the 65536 binary genomes, which is EXACTLY sf.f's kernel, so it
is the same stochastic process as PartnersRivals1's run F1 (DiskM2WF dw1 / PartnersRivals Cannon f1, N = 1000, beta = 100,
mu = 1e-2, eps = 1e-4, 10 x 1e7 generations), on different random numbers.  Per game its ten replicates must be
statistically indistinguishable from F1's ten, whose per-replicate efficiency and atom shares are PartnersRivals/Data/f1_e4.rep.
Printed per game: the mean and replicate s.d. of both, and Welch's t for the efficiency and for each atom share."""
import sys, glob, os, math
import numpy as np
# usage: wfc_probe.py <dir of r_500_*.txt> [F1|F2]   -- F1: N = 1000 against f1_e4.rep; F2: N = 100 against f2_e4.rep
RUN = sys.argv[2] if len(sys.argv) > 2 else "F1"
REP = "/Users/martin/Documents/CL1/PartnersRivals/Data/%s_e4.rep" % RUN.lower()
ORD = ["000", "100", "010", "001", "110", "011", "111"]            # Revised2's order, that of f1_e4.rep
ref = {}
for l in open(REP):
    if l.startswith("#") or not l.strip(): continue
    f = l.split(); gid = int(f[0]); ref.setdefault(gid, []).append([float(x) for x in f[2:10]])
def mine(d, am):
    out = {}
    for fn in glob.glob(os.path.join(d, "r_%d_*.txt" % am)):
        for l in open(fn):
            f = l.split()
            if not f or f[0] != "REP": continue
            gid = int(f[1]); ns = float(f[10]); a = [float(x) for x in f[12:20]]
            if a[5] != 0: raise SystemExit("atom 101 occupied at gid %d" % gid)
            sh = {"%d%d%d" % (k >> 2, (k >> 1) & 1, k & 1): a[k] / ns for k in range(8)}
            if len(f) > 30 and int(f[26]):                    # wedge W: the exact 111 counter replaces the sample
                sh["111"] = float(f[27]) / float(f[30])
            out.setdefault(gid, []).append([float(f[7])] + [sh[o] for o in ORD])
    return out
def welch(x, y):
    x, y = np.asarray(x), np.asarray(y)
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
    if vx + vy == 0: return 0.0 if x.mean() == y.mean() else float("inf")
    return (x.mean() - y.mean()) / math.sqrt(vx + vy)
M = mine(sys.argv[1], 500)
worst = 0.0
for gid in sorted(M):
    a, b = np.array(M[gid]), np.array(ref[gid])
    print("gid %d   replicates: wfc %d, %s %d" % (gid, len(a), RUN, len(b)))
    print("   %-5s %22s %22s %8s" % ("", "wfc alpha=1/2", "binary " + RUN, "Welch t"))
    for k, name in enumerate(["eff"] + ORD):
        t = welch(a[:, k], b[:, k]); worst = max(worst, abs(t)) if np.isfinite(t) else worst
        print("   %-5s %12.5f +- %7.5f %12.5f +- %7.5f %8.2f" % (name, a[:, k].mean(), a[:, k].std(ddof=1),
                                                            b[:, k].mean(), b[:, k].std(ddof=1), t))
    ma, mb = ORD[int(np.argmax(a[:, 1:].mean(0)))], ORD[int(np.argmax(b[:, 1:].mean(0)))]
    print("   most abundant atom: wfc %s   %s %s   %s" % (ma, RUN, mb, "agree" if ma == mb else "DISAGREE"))
print("run %s: largest |Welch t| over all games and quantities: %.2f" % (RUN, worst))
