#!/usr/bin/env python3
"""check_v4.py [--n N] <dir> <isl>... -- V4 checks on continuum outputs (<isl>, w<isl>, h<isl>) and, if
present, kappa lines k<isl> (kappa.x output).  Exits non-zero on any failure.
--n: the population size (default 1000, run F1; --n 100 for run F2), so maxc = 2N+2 is the table size."""
import sys, os
args = sys.argv[1:]
NPOP = 1000
if args[:1] == ["--n"]:
    NPOP = int(args[1]); args = args[2:]
MAXC = 2 * NPOP + 2
d = args[0]; bad = 0
for isl in args[1:]:
    H = open(os.path.join(d, "h" + isl)).read().split("\n")
    sec, cnt, tot, itend, nper = None, {}, None, None, None
    for l in H:
        if l.startswith("# run:"):
            f = l.split(); nper, itend = int(f[3]), int(f[5])
        elif l.startswith("# atoms"): sec = "a"; cnt[sec] = []
        elif l.startswith("# strata"): sec = "s"; cnt[sec] = []
        elif l.startswith("# corners"): sec = "c"; cnt[sec] = []
        elif l.startswith("# total"): tot = int(l.split()[2]); sec = None
        elif l.strip() and not l.startswith("#"): cnt[sec].append(int(l.split()[1]))
    exp = nper * (itend // 2) * NPOP
    a, s, c = cnt["a"], cnt["s"], cnt["c"]
    corner_strata = [sum(((k >> j) & 1) * 3 ** j for j in range(4)) for k in range(16)]
    ok_h = (len(a), len(s), len(c)) == (8, 81, 16) and sum(a) == sum(s) == tot == exp and a[5] == 0 \
        and all(c[k] == s[corner_strata[k]] for k in range(16))
    R = [l.split() for l in open(os.path.join(d, isl)) if l.strip()]
    dev = max(abs(sum(float(x) for x in r[9:17]) - 1) for r in R)
    s101 = max(float(r[9 + 5]) for r in R)
    W = [l.split() for l in open(os.path.join(d, "w" + isl)) if l.startswith("# ") and l.split()[1].isdigit()]
    nhw = max(int(w[7]) for w in W)
    nnear = sum(int(w[5]) for w in W); ntie = sum(int(w[12]) for w in W)
    ncls = sum(int(w[3]) for w in W); ncall = sum(int(w[4]) for w in W)
    line = ("isl %s: reps %d, h total %d (expected %d), sum atoms %d, sum strata %d, corners=corner strata %s, "
            "max|sum share-1| %.1e, max s101 %g, nhw %d/%d, lineages classified %d, m1cls calls %d, "
            "near ties %d, exact ties %d" % (isl, len(R), tot, exp, sum(a), sum(s),
            all(c[k] == s[corner_strata[k]] for k in range(16)), dev, s101, nhw, MAXC, ncls, ncall, nnear, ntie))
    good = ok_h and len(R) == nper and dev < 1e-9 and s101 == 0 and nhw <= MAXC
    kp = os.path.join(d, "k" + isl)
    if os.path.exists(kp):
        f = open(kp).read().split()
        ka = [float(x) for x in f[7:15]]; kz = [float(x) for x in f[15:96]]
        line += "; kappa M=%s: sum atoms-1 %.1e, sum strata-1 %.1e, k101 %g" % (f[6], sum(ka) - 1, sum(kz) - 1, ka[5])
        good = good and abs(sum(ka) - 1) < 1e-9 and abs(sum(kz) - 1) < 1e-9 and ka[5] == 0
    print(("PASS " if good else "FAIL ") + line)
    bad += not good
sys.exit(bad)
