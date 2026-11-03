#!/usr/bin/env python3
"""
test_cls.py -- proves m1cls.c against two independent references.

    python3 test_cls.py            build, then V1, V2, the Fortran link test, and kappa sanity
    python3 test_cls.py v1|v2|fortran|kappa     one part only
    options:  --n2 N    points per non-corner stratum in V2 (default 48 -> 3120, + 160 corner points)

V1  the 16 corners x 512 games against DiskM1WF/Data/m1_atoms_512.npz (read-only).  Layout of the npz,
    learned by inspection and from DiskM1WF/Opt/m1disk.py: codes[ipt, k] (512 x 16, int8) is the atom
    4*eff + 2*nash + comp of strategy k at game ipt, k in sf1.f's bit convention (bit j of k = C after
    state j, CC CD DC DD), games in games.dat row order (uv[ipt] = (cS, cT - 1)).  Exact rationals there
    (Markov chain tree theorem, m1atoms.py); leading-order GTH + closed forms here.

V2  >= 3000 non-corner points over all 81 strata at 24 games from games.dat (6 per wedge W, S, N, E)
    plus 5 games on the boundary lines, against an EXACT oracle built only from check/m1lib.py:
      * the preimage p of the classifier's input s, in Fractions (s == eps -> 0, s == 1-eps -> 1,
        else (s - eps)/(1 - 2 eps) exactly);
      * N: m1lib.stationary_exact at eps_o = 10^-100 (exact rationals), Nash over the 16 deterministic
        deviators; a non-tie is decided outright, an exact limit tie shows as O(1e-72) or less -- the
        decision threshold is 1e-60*scale (eps_o = 1e-40 is too coarse for the near-face draws, whose
        genuine gaps go down to ~d^2 ~ 1e-28);
      * E, R: m1lib's closed forms (in_eff_cc / in_eff_alt / in_eff_line, in_riv_plus / in_riv_minus,
        all on T = S) on the exact p.
    As a cross-check of the THEORY, E and R are also decided from their definitions at eps_o
    (pi(p,p) = Emax; min over the 16 binary co-players of sign(T-S)(w_DC - w_CD) >= 0) and compared with
    the closed forms.  Interior coordinates are drawn four ways: the production kernel, dyadic p = k/16,
    within 1e-6..1e-14 of a face, and s3 = 1 - s2 exactly (p2 + p3 = 1, the rival boundary).
    A separate block of CONSTRUCTED near-ties (the ZD equaliser sheet and the diagonal x* = u/(u+v)
    of Text/nash.md sec.1, which the double preimage can only approximate) is reported on its own.

READ-ONLY: DiskM1WF/*, check/m1lib.py.  Builds into ./_build_test/.
"""
import os
import random
import subprocess
import sys
import time
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "_build_test")
CHECK = os.path.normpath(os.path.join(HERE, os.pardir, os.pardir, "check"))
sys.dont_write_bytecode = True          # never write a .pyc into the read-only check/ tree
sys.path.insert(0, CHECK)
import m1lib  # noqa: E402

GAMES = "/Users/martin/Documents/CL1/DiskM1WF/Cannon/mw/games.dat"
NPZ = "/Users/martin/Documents/CL1/DiskM1WF/Data/m1_atoms_512.npz"
EPS = 1e-4
EPS_O = F(1, 10 ** 100)
# oracle tie threshold: an exact limit tie shows as O(C eps_o) with C <= ~1/d^2 <= 1e28 for the
# near-face draws (d >= 1e-14), i.e. <= 1e-72; a genuine gap is >= ~d^2 >= 1e-28.  1e-60 separates them
# by twelve orders on each side, and the report prints both extremes.
THR_O = F(1, 10 ** 60)
NAMES = ["000", "001", "010", "011", "100", "101", "110", "111"]


def sh(cmd, **kw):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    if r.returncode:
        raise SystemExit("FAILED: %s\n%s%s" % (cmd, r.stdout, r.stderr))
    return r.stdout


# on this Mac python3 (anaconda) is x86_64 under Rosetta, so a cc it spawns targets x86_64 while
# gfortran is arm64: force the native arch whenever the hardware is Apple silicon
def _cc():
    if sys.platform == "darwin":
        r = subprocess.run(["sysctl", "-n", "hw.optional.arm64"], capture_output=True, text=True)
        if r.stdout.strip() == "1":
            return "cc -arch arm64"
    return "cc"


CC = _cc()


def build():
    os.makedirs(BUILD, exist_ok=True)
    src = HERE
    sh("{cc} -O2 -Wall -Wextra -o {b}/test_cls.x {s}/test_cls.c {s}/m1cls.c -lm".format(cc=CC, b=BUILD, s=src))
    sh("{cc} -O2 -Wall -Wextra -o {b}/kappa.x {s}/kappa.c {s}/m1cls.c -lm".format(cc=CC, b=BUILD, s=src))
    sh("{b}/test_cls.x selftest".format(b=BUILD))


def read_games():
    rows = []
    for line in open(GAMES):
        t = line.split()
        if not t or t[0].startswith("#"):
            continue
        cR, cS, cT = float(t[2]), float(t[3]), float(t[4])
        rows.append((int(t[0]), cS, cT - 1.0, t[7]))
    return rows


# ------------------------------------------------------------------------------------------ V1
def v1():
    z = np.load(NPZ)
    rows = read_games()
    uv = np.array([(r[1], r[2]) for r in rows])
    assert uv.shape == (512, 2) and np.array_equal(uv, z["uv"]), "games.dat and the npz disagree on (u,v)"
    out = sh("{b}/test_cls.x corners {g}".format(b=BUILD, g=GAMES))
    a = np.array([list(map(int, l.split())) for l in out.splitlines()])
    assert a.shape == (8192, 5)
    codes = a[:, 2].reshape(512, 16)
    ref = z["codes"].astype(int)
    mism = int((codes != ref).sum())
    print("V1  corners x games: %d comparisons, %d mismatches   (early exit == full test at all 8192)"
          % (codes.size, mism))
    nt = a[:, 4]
    print("    neartie flags at the corners: 0 -> %d, 1 (rounding-sensitive) -> %d, 2 (exact structural tie) -> %d"
          % ((nt == 0).sum(), (nt == 1).sum(), (nt == 2).sum()))
    print("    atom counts summed over the 512 games (ENR order %s):" % " ".join(NAMES),
          np.bincount(codes.ravel(), minlength=8).tolist(), " npz atoms colsum:",
          z["atoms"].sum(0).tolist())
    if mism:
        bad = np.argwhere(codes != ref)
        for g, k in bad[:20]:
            print("    MISMATCH game %d (%s) strategy %d: C %s  npz %s" % (g, rows[g][3], k, NAMES[codes[g, k]], NAMES[ref[g, k]]))
    return mism


# ------------------------------------------------------------------------------------------ oracle
def preimage(s):
    """exact preimage of the classifier's input, and the stratum code"""
    om = 1.0 - EPS
    p, c = [], []
    e = F(EPS)
    for x in s:
        if x == EPS:
            p.append(F(0)); c.append(0)
        elif x == om:
            p.append(F(1)); c.append(1)
        else:
            p.append((F(x) - e) / (1 - 2 * e)); c.append(2)
    return p, sum(ci * 3 ** k for k, ci in enumerate(c))


def oracle(s, u, v):
    p, st = preimage(s)
    U, V = F(u), F(v)
    scale = max(1.0, abs(u), abs(1 + v))
    thr = THR_O * F(scale)
    wpp = m1lib.stationary_exact(p, p, EPS_O)
    own = m1lib.payoff(wpp, U, V)
    gaps, rivd = [], []
    for q in m1lib.DET_EXACT:
        w = m1lib.stationary_exact(q, p, EPS_O)        # deviator q focal, p the co-player
        gaps.append(m1lib.payoff(w, U, V) - own)
        rivd.append(w[1] - w[2])                        # p's side: w_DC - w_CD = (q side) w_CD - w_DC
    mx = max(gaps)
    N = mx <= thr
    sw = U + V - 1
    if sw < 0:
        E = m1lib.in_eff_cc(p)
    elif sw > 0:
        E = m1lib.in_eff_alt(p)
    else:
        E = m1lib.in_eff_line(p)
    ts = U - V - 1                                      # < 0: T > S
    R = True if ts == 0 else (m1lib.in_riv_plus(p) if ts < 0 else m1lib.in_riv_minus(p))
    emax = max(F(1), (1 + U + V) / 2)
    E_def = own - emax >= -thr
    if ts == 0:
        R_def = True
    else:
        sg = 1 if ts < 0 else -1
        R_def = min(sg * d for d in rivd) >= -thr
    return dict(E=int(E), N=int(N), R=int(R), st=st, mx=mx, E_def=int(E_def), R_def=int(R_def),
                atom=4 * int(E) + 2 * int(N) + int(R))


def run_c_points(pts):
    lines = "".join("%s %s %s %s %s %s %s\n" % (u.hex(), v.hex(), EPS.hex(), *(x.hex() for x in s))
                    for (u, v, s) in pts)
    r = subprocess.run([os.path.join(BUILD, "test_cls.x"), "points"], input=lines, capture_output=True,
                       text=True)
    if r.returncode:
        raise SystemExit("test_cls.x points failed (early exit != full?): " + r.stderr)
    out = []
    for l in r.stdout.splitlines():
        t = l.split()
        out.append(dict(atom=int(t[0]), st=int(t[1]), nt=int(t[2]), E=int(t[3]), N=int(t[4]),
                        R=int(t[5]), mx=float.fromhex(t[6]), atom_e=int(t[7]), st_e=int(t[8]),
                        nt_e=int(t[9]), nd_e=int(t[10])))
    assert len(out) == len(pts)
    return out


def wedge(u, v):
    a, b = u + v - 1, u - v - 1
    if a == 0 or b == 0:
        return "line"
    return {(True, True): "W", (True, False): "S", (False, True): "N", (False, False): "E"}[(a < 0, b < 0)]


def v2_games(rng):
    rows = read_games()
    by = {"W": [], "S": [], "N": [], "E": []}
    for r in rows:
        w = wedge(r[1], r[2])
        if w in by:
            by[w].append(r)
    games = []
    for w in "WSNE":
        pick = rng.sample(by[w], 6)
        games += [(r[1], r[2], w + ":" + r[3]) for r in pick]
    # the boundary loci of Text/friendly.md sec.2, exactly representable
    games += [(0.25, 0.75, "line:switch,T>S"), (1.5, -0.5, "line:switch,T<S"),
              (0.25, -0.75, "line:T=S,below"), (3.0, 2.0, "line:T=S,above"), (1.0, 0.0, "line:corner(1,0)")]
    return games


def draw_interior(rng, pmin, pmax, ua, uw):
    k = rng.random()
    if k < 0.5:                                   # the production kernel's interior
        while True:
            x = ua + rng.random() * uw
            if pmin < x < pmax:
                return x, "kernel"
    if k < 0.7:                                   # dyadic p = j/16
        p = rng.randint(1, 15) / 16.0
        return EPS + (1 - 2 * EPS) * p, "dyadic"
    d = 10.0 ** -rng.uniform(6, 14)               # within 1e-6 .. 1e-14 of a face
    x = (pmin + d) if rng.random() < 0.5 else (pmax - d)
    if x <= pmin or x >= pmax:
        x = (pmin + pmax) / 2
    return x, "nearface"


def v2(n_per=48, seed=20260924, n_corner=10):
    rng = random.Random(seed)
    pmin, pmax = EPS, 1.0 - EPS
    a = 1.0 / 3.0
    uw = (1.0 - 2.0 * pmin) / (1.0 - 2.0 * a)
    ua = pmin - a * uw
    games = v2_games(rng)
    pts, meta = [], []
    gi = 0
    for st in range(81):
        code = [(st // 3 ** k) % 3 for k in range(4)]
        # the 16 corner strata are V1's in full; here they get a few points at the V2 games (incl.
        # the boundary loci V1 does not visit) so that all 81 strata go through the same oracle
        for j in range(n_per if 2 in code else n_corner):
            s, kinds = [], set()
            for c in code:
                if c < 2:
                    s.append(pmin if c == 0 else pmax)
                else:
                    x, kd = draw_interior(rng, pmin, pmax, ua, uw)
                    s.append(x); kinds.add(kd)
            if code[1] == 2 and code[2] == 2 and rng.random() < 0.3:
                s[1] = rng.uniform(0.5, pmax)     # s3 = 1 - s2 exactly: p2 + p3 = 1 exactly
                s[2] = 1.0 - s[1]
                assert pmin < s[2] < pmax and s[1] + s[2] == 1.0
            u, v, tag = games[gi % len(games)]
            gi += 1
            pts.append((u, v, s))
            meta.append((st, tag, "tie23" if (code[1] == 2 and code[2] == 2 and s[1] + s[2] == 1.0) else "",
                         "nearface" in kinds))
    t0 = time.time()
    C = run_c_points(pts)
    tc = time.time() - t0
    t0 = time.time()
    O = [oracle(s, u, v) for (u, v, s) in pts]
    to = time.time() - t0
    nm = dict(atom=0, E=0, N=0, R=0, st=0, early=0, Edef=0, Rdef=0, expl=0)
    bad = []
    for (u, v, s), c, o, m in zip(pts, C, O, meta):
        nm["st"] += c["st"] != o["st"]
        nm["E"] += c["E"] != o["E"]
        nm["N"] += c["N"] != o["N"]
        nm["R"] += c["R"] != o["R"]
        nm["atom"] += c["atom"] != o["atom"]
        nm["early"] += (c["atom"], c["st"], c["nt"]) != (c["atom_e"], c["st_e"], c["nt_e"])
        nm["Edef"] += o["E"] != o["E_def"]
        nm["Rdef"] += o["R"] != o["R_def"]
        scale = max(1.0, abs(u), abs(1 + v))
        # EXPLAINED: only N differs, C flagged a near-tie, and the exact limit margin is itself inside
        # the band (|max gap| <= 1e-9*scale): a genuine but sub-tolerance margin
        expl = (c["atom"] != o["atom"] and c["E"] == o["E"] and c["R"] == o["R"] and c["nt"] != 0
                and abs(float(o["mx"])) <= 1e-9 * scale)
        nm["expl"] += expl
        if (c["atom"] != o["atom"] and not expl) or o["E"] != o["E_def"] or o["R"] != o["R_def"]:
            bad.append((u, v, s, c, o, m))
    strata = sorted({m[0] for m in meta})
    ncorner = sum(1 for m in meta if all((m[0] // 3 ** k) % 3 < 2 for k in range(4)))
    wedges = {}
    for (u, v, s), m in zip(pts, meta):
        wedges[m[1].split(":")[0]] = wedges.get(m[1].split(":")[0], 0) + 1
    print("V2  %d points (%d non-corner + %d corner), %d strata, %d games %s; C %.1fs, exact oracle %.1fs"
          % (len(pts), len(pts) - ncorner, ncorner, len(strata), len(games), wedges, tc, to))
    print("    mismatches C vs oracle: atom %d  (E %d, N %d, R %d, stratum %d);  early-exit != full: %d"
          % (nm["atom"], nm["E"], nm["N"], nm["R"], nm["st"], nm["early"]))
    print("    of the atom mismatches, EXPLAINED (N only, C neartie-flagged, exact |margin| <= 1e-9*scale): %d;"
          " UNEXPLAINED: %d" % (nm["expl"], nm["atom"] - nm["expl"]))
    ex = [(m[0], c["mx"], float(o["mx"]), min(min(abs(x - EPS), abs(1 - EPS - x)) for x in s if x not in (EPS, 1 - EPS)))
          for (u, v, s), c, o, m in zip(pts, C, O, meta) if c["atom"] != o["atom"]]
    for e in ex[:12]:
        print("      stratum %2d  C max-gap %+.3e  exact %+.3e  nearest interior coordinate to a face: %.1e (in s)" % e)
    print("    theory cross-check (m1lib closed form vs definition at eps=1e-100): E %d, R %d mismatches"
          % (nm["Edef"], nm["Rdef"]))
    ntie = sum(1 for m in meta if m[2])
    nnf = sum(1 for m in meta if m[3])
    mm_nf = sum(1 for c, o, m in zip(C, O, meta) if c["atom"] != o["atom"] and m[3])
    print("    points with a coordinate within 1e-6..1e-14 of a face: %d, carrying %d of the %d atom mismatches;"
          " all other points (kernel / dyadic / p2+p3=1 draws): %d mismatches"
          % (nnf, mm_nf, nm["atom"], nm["atom"] - mm_nf))
    at = np.bincount([o["atom"] for o in O], minlength=8)
    print("    exact p2+p3=1 points: %d;  oracle atom counts %s" % (ntie, dict(zip(NAMES, at.tolist()))))
    nt = np.array([c["nt"] for c in C])
    # separation: the smallest |oracle max-gap| among decided non-ties and the largest among ties
    mxs = [abs(float(o["mx"])) for o in O]
    small = sorted(x for x in mxs if x > 1e-60)[:1]
    tiny = max((x for x in mxs if x <= 1e-60), default=0.0)
    print("    neartie flags: 0 -> %d, 1 -> %d, 2 -> %d;  oracle |max gap|: smallest non-tie %.3g, largest tie %.3g"
          % ((nt == 0).sum(), (nt == 1).sum(), (nt == 2).sum(), small[0] if small else float("nan"), tiny))
    # the neartie flag against the oracle: agreement of C's float gap with the exact limit gap
    dev = max(abs(c["mx"] - float(o["mx"])) for c, o in zip(C, O))
    print("    max |C max-gap - exact max-gap| over all points: %.3g" % dev)
    for (u, v, s, c, o, m) in bad[:30]:
        print("    DISAGREE game (%s) stratum %d %s s=%s C=%s O=%s E_def=%d R_def=%d Cmx=%.3g Omx=%.3g nt=%d"
              % (m[1], m[0], m[2], [float.hex(x) for x in s], NAMES[c["atom"]], NAMES[o["atom"]],
                 o["E_def"], o["R_def"], c["mx"], float(o["mx"]), c["nt"]))
    ok = nm["atom"] == nm["expl"] and nm["st"] == 0 and nm["early"] == 0 and nm["Edef"] == 0 and nm["Rdef"] == 0
    v2_constructed(rng, games)
    return 0 if ok else 1


def v2_constructed(rng, games):
    """ZD equaliser sheet and the diagonal x*: genuine limit ties the double preimage can only
    approximate to ~1e-16.  Reported separately: every C/oracle disagreement must carry a neartie flag."""
    pts = []
    for (u, v, tag) in games:
        if not tag.startswith(("W", "S", "N", "E")):
            continue
        if u > 0 and v > 0:                           # x* = u/(u+v) on the diagonal
            x = u / (u + v)
            sx = EPS + (1 - 2 * EPS) * x
            pts.append((u, v, [sx] * 4, "x*"))
        # ZD sheet: A(p1-p2) = v, A(p4-p3) = u, A(1+p3-p1) = 1-u ; pick A and p1, solve
        for _ in range(4):
            A = rng.uniform(max(abs(u), abs(v), abs(1 - u)) * 1.05, max(abs(u), abs(v), abs(1 - u)) * 4 + 1)
            p1 = rng.uniform(0.3, 0.99)
            p3 = p1 - 1 + (1 - u) / A
            p2 = p1 - v / A
            p4 = p3 + u / A
            p = [p1, p2, p3, p4]
            if all(0 < x < 1 for x in p):
                pts.append((u, v, [EPS + (1 - 2 * EPS) * x for x in p], "ZD"))
    if not pts:
        return
    C = run_c_points([(u, v, s) for (u, v, s, _) in pts])
    dis, flagged, nflag, worst = 0, 0, 0, 0.0
    for (u, v, s, k), c in zip(pts, C):
        o = oracle(s, u, v)
        nflag += c["nt"] != 0
        if c["atom"] != o["atom"]:
            dis += 1
            flagged += c["nt"] != 0
            worst = max(worst, abs(float(o["mx"])))
    print("    constructed near-ties (ZD sheet, x*): %d points, neartie-flagged %d; C/oracle disagreements %d, "
          "of which neartie-flagged %d; largest exact |margin| among them %.2e" % (len(pts), nflag, dis, flagged, worst))


# ------------------------------------------------------------------------------------------ Fortran
def fortran():
    # m1cls.c compiled twice -- by clang, and by gfortran's own GCC C front end -- and each object
    # linked with both Fortran call styles
    sh("{cc} -O2 -c -o {b}/m1cls_clang.o {s}/m1cls.c".format(cc=CC, b=BUILD, s=HERE))
    sh("gfortran -O2 -c -o {b}/m1cls_gcc.o {s}/m1cls.c".format(b=BUILD, s=HERE))
    res = {}
    for obj in ("clang", "gcc"):
        for t, src in (("fcall", "test_fcall.f"), ("fbind", "test_fbind.f90")):
            x = "{b}/test_{t}_{o}.x".format(b=BUILD, t=t, o=obj)
            sh("gfortran -O2 -o {x} {s}/{src} {b}/m1cls_{o}.o".format(x=x, s=HERE, src=src, b=BUILD, o=obj))
            res[(t, obj)] = sh(x).split()
    a, b = res[("fcall", "gcc")], res[("fbind", "gcc")]
    syms = sh("nm {b}/m1cls_gcc.o | grep ' T '".format(b=BUILD))
    pts = [(-0.5, 0.5, [EPS] * 4), (-0.1, 0.1, [1 - EPS, EPS, EPS, 1 - EPS]),
           (1.5, -0.25, [1 - EPS, 0.7, 0.6, 0.2]), (-0.5, 0.5, [1 - EPS, 0.75, 0.25, EPS])]
    C = run_c_points(pts)
    want = [str(x) for c in C for x in (c["atom_e"], c["st_e"], c["nt_e"])]
    ok = all(r == want for r in res.values())
    print("Fortran link (gfortran x {clang, gcc} objects): implicit external (m1cls_) %s, bind(C) (m1cls) %s"
          " -> %s;  exported: %s"
          % (" ".join(a), " ".join(b), "OK, all 4 builds agree with C" if ok else "MISMATCH vs C " + " ".join(want),
             " ".join(l.split()[-1] for l in syms.splitlines())))
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------ kappa
def kappa(M=1000000):
    rows = read_games()
    # four representative games, one per wedge: the nearest simulated game to a reference point
    ref = {"W (PD)": (-0.5, 0.5), "S": (1.5, -0.5), "N": (0.5, 1.5), "E": (3.0, 0.5)}
    bad = 0
    print("kappa (kernel measure of the atoms, eps = 1e-4, M = %d draws; atoms in ENR order):" % M)
    for name, (u0, v0) in ref.items():
        ig = min(range(len(rows)), key=lambda i: (rows[i][1] - u0) ** 2 + (rows[i][2] - v0) ** 2)
        isl = 3 * ig + 3
        out = subprocess.run([os.path.join(BUILD, "kappa.x"), str(isl), str(M)], cwd=os.path.dirname(GAMES),
                             capture_output=True, text=True, check=True).stdout.split()
        assert len(out) == 98, len(out)
        fa = np.array(list(map(float, out[7:15])))
        fs = np.array(list(map(float, out[15:96])))
        ok = abs(fa.sum() - 1) < 1e-9 and abs(fs.sum() - 1) < 1e-9 and fa[5] == 0
        bad += not ok
        print("  %-7s isl %4d (u,v)=(%.4f,%.4f) %s  sum=%.12g  strata 1/81 dev %.4f  near1 %s tie2 %s"
              % (name, isl, float(out[4]), float(out[5]),
                 " ".join("%s=%.4f" % (n, x) for n, x in zip(NAMES, fa)), fa.sum(),
                 np.abs(fs - 1 / 81).max() * 81, out[96], out[97]))
    return bad


def main():
    args = sys.argv[1:]
    n2 = 48
    if "--n2" in args:
        i = args.index("--n2"); n2 = int(args[i + 1]); del args[i:i + 2]
    part = args[0] if args else "all"
    build()
    rc = 0
    if part in ("all", "v1"):
        rc |= v1() != 0
    if part in ("all", "v2"):
        rc |= v2(n2)
    if part in ("all", "fortran"):
        rc |= fortran()
    if part in ("all", "kappa"):
        rc |= kappa() != 0
    print("RESULT", "PASS" if rc == 0 else "FAIL")
    sys.exit(rc)


if __name__ == "__main__":
    main()
