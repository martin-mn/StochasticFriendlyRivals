#!/usr/bin/env python3
"""
sitable1.py -- SI Table 1 of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity":
the friendly-rival strata of the memory-two cube, every cell of panels a and b computed from the deposited census output.

    python3 sitable1.py   ->  sitable1.out   the two panels, each cell labelled as in the table; the census's own counts,
                                             the strata the probe added, the dimension histograms, and the checks

A stratum is a support pattern of the 16 coordinates, in the census's code g = sum_j d_j 3^j (d_j = 0: p_j = 0, 1: p_j
free in (0,1), 2: p_j = 1); its dimension is the number of free coordinates.  Inputs (CL1 is the directory that holds
PRFull/; relocate.py points it at the release):
  PartnersRivalsM2f/check/st_eff.txt        the census (Cannon array 45751662, census.c, K = 64 draws per stratum): a line
                                            "g class nP nM K" for every efficient stratum with a confirmed friendly rival,
                                            class C (w_CC = 1), A (w_CD + w_DC = 1) or Z (w_DD = 0 only), nP and nM the
                                            draws confirmed for T > S and for T < S
  PartnersRivalsM2f/cannon/ca1_out/         the dimension probe (array 45752239, dimprobe.c via dim.sbatch: every stratum
    dim_<class>_<sign>_<stripe>.txt         of dimension >= 5 of the class, up to 512 draws; stripe = g mod 64): a line
    dim_45752239_<task>.err                 "g dimension class" for each stratum with a friendly rival, and each task's
                                            count.  On GitHub and in the Zenodo code archive; in the author's tree
                                            PRFull/release/harvest/ca1_out, where assemble_release.py collected it
  PRFull/release/census_class/census_support_class.u8.gz
                                            every stratum's efficiency class (all 3^16 of them), read off the support of
                                            its self-play limit by check/ubound/ub.c (mode classes): the rows T = S, (1,0)
  PartnersRivals/Figures/m2_masks.npz       the binary friendly rivals (effCC, effALT, rivP, rivM)
  PartnersRivalsM2f/check/binary_eff.npz    effLINE: the binary strategies efficient on the switch line

The loci are those of PartnersRivalsM2f/check/census.py: wedge W is class C with T > S, S is C with T < S, N is A with
T > S, E is A with T < S; the ray of the switch line with T > S (T < S) takes the classes C, A and Z with T > S (T < S).
Panel a counts the strata of the census together with those of the probe.  A stratum is binary-invisible if its closed
cell holds no binary friendly rival of its wedge (census.py's visible_flags); in S this is counted over the census's
strata only, as in the table.  On the line T = S every efficient stratum holds friendly rivals: class C below the switch
line, class A above it, and every class at (1,0).  Panel b sets the census's strata of dimension >= 5 against the
probe's, block by block; the deficit counts the probe's strata that are absent from the census, and the rows of the
rays take class Z only, whose strata are efficient on the switch line alone.

Read-only; writes only sitable1.out beside itself, and reads no manuscript.  Every check (listed at the end of
sitable1.out) must pass, or the script says FAILED and exits with status 1.
"""
import collections
import gzip
import math
import os
import re
import sys

sys.dont_write_bytecode = True
import numpy as np                                                 # noqa: E402

CL1 = "/Users/martin/Documents/CL1"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sitable1.out")
M2F = os.path.join(CL1, "PartnersRivalsM2f")
STEFF = os.path.join(M2F, "check", "st_eff.txt")
CA1 = os.path.join(M2F, "cannon", "ca1_out")                       # the release
if not os.path.isdir(CA1):
    CA1 = os.path.join(CL1, "PRFull", "release", "harvest", "ca1_out")   # the author's tree
CA1_NAME = "PartnersRivalsM2f/cannon/ca1_out"                      # as printed, in either tree
CLASSES = os.path.join(CL1, "PRFull", "release", "census_class", "census_support_class.u8.gz")
MASKS = os.path.join(CL1, "PartnersRivals", "Figures", "m2_masks.npz")
BINEFF = os.path.join(M2F, "check", "binary_eff.npz")

N3 = 3 ** 16
DIMJOB = 45752239
NSTRIPE = 64
BLOCKS = [("C", 1), ("C", -1), ("A", 1), ("A", -1), ("Z", 1), ("Z", -1)]   # dim.sbatch: task = 64 * block + stripe
CODE = {"C": 1, "A": 2, "Z": 3}                                    # census_support_class: 0 none, 1 C, 2 A, 3 Z
LOCI = [("W", [("C", 1)]), ("S", [("C", -1)]), ("N", [("A", 1)]), ("E", [("A", -1)]),
        ("ray T>S", [("C", 1), ("A", 1), ("Z", 1)]), ("ray T<S", [("C", -1), ("A", -1), ("Z", -1)])]
BROWS = [("wedge W", ("C", 1)), ("wedge S", ("C", -1)), ("wedge N", ("A", 1)), ("wedge E", ("A", -1)),
         ("ray T>S, efficient only there", ("Z", 1)), ("ray T<S, efficient only there", ("Z", -1))]
CHECKS = []


def check(what, ok):
    CHECKS.append((what, bool(ok)))


def num(n):
    """a count as the table prints it: digits grouped by spaces from 10 000 on"""
    return "{:,}".format(n).replace(",", " ") if n >= 10000 else str(n)


def pct2(x):
    """a share in per cent to two significant digits (panel a, and the probe shares of the legend)"""
    if x == 0:
        return "0"
    return "%.*f%%" % (max(0, 1 - math.floor(math.log10(100 * x))), 100 * x)


def pct_b(x):
    """a share of panel b: per cent to two decimals, 0 as 0"""
    return "0" if x == 0 else "%.2f%%" % (100 * x)


def digits(g):
    """(16, n): the census digits of the strata g"""
    g = np.array(g, dtype=np.int64, ndmin=1)
    d = np.empty((16, g.size), np.int8)
    for j in range(16):
        d[j] = g % 3
        g = g // 3
    return d


def nfree(g):
    return (digits(g) == 1).sum(0)


def pinned(g):
    """M, V: the closed cell of stratum g holds the binary genome b (bit j set: p_j = 1) iff b & M == V"""
    d = digits(g)
    M = np.zeros(d.shape[1], np.int64)
    V = np.zeros(d.shape[1], np.int64)
    for j in range(16):
        M |= (d[j] != 1).astype(np.int64) << j
        V |= (d[j] == 2).astype(np.int64) << j
    return M, V


def mirror(g):
    """the C/D mirror of the strata g, sorted: state s -> 15 - s and p -> 1 - p (census.py's mirror)"""
    d = digits(g)
    m = np.zeros(d.shape[1], np.int64)
    for j in range(16):
        m += (2 - d[15 - j]).astype(np.int64) * 3 ** j
    return np.sort(m)


def code_of(b):
    """the census code of the binary genomes b, the strata of dimension 0"""
    b = np.asarray(b, np.int64)
    return sum(2 * ((b >> j) & 1) * 3 ** j for j in range(16))


def pat(g):
    return "".join("0*1"[x] for x in digits(g)[:, 0])


def invisible(S, ref):
    """the number of strata of S whose closed cell holds no genome of ref, grouped by the pin mask as in census.py"""
    M, V = pinned(S)
    inv = 0
    order = np.argsort(M, kind="stable")
    for idx in np.split(order, np.flatnonzero(np.diff(M[order])) + 1):
        if idx.size:
            inv += int((~np.isin(V[idx], np.unique(ref & M[idx[0]]))).sum())
    return inv


def hist(dims):
    h = np.bincount(np.asarray(dims), minlength=17)
    nz = np.flatnonzero(h)
    return h, int(nz[0]), int(nz[-1])


def hist_txt(h):
    return " ".join("%d:%d" % (k, h[k]) for k in np.flatnonzero(h))


def class_dims(cls):
    """H[class, dimension]: the strata of each class (0 none, 1 C, 2 A, 3 Z) by dimension, over all 3^16 strata in code
    order; and the strata of dimension 0 of each class.  In blocks of 3^12 strata, one per value of the top 4 digits."""
    low = np.zeros(1, np.int8)                                     # the dimensions of the 3^12 values of the low digits
    for _ in range(12):
        low = np.concatenate([low, low + 1, low])
    top = nfree(np.arange(81))
    H, Z0 = np.zeros((4, 17), np.int64), {1: [], 2: [], 3: []}
    for t in range(81):
        cl = cls[t * low.size:(t + 1) * low.size]
        H += np.bincount(cl.astype(np.int64) * 17 + low + top[t], minlength=68).reshape(4, 17)
        if top[t] == 0:
            for k in Z0:
                Z0[k].append(t * low.size + np.flatnonzero((low == 0) & (cl == k)))
    return H, {k: np.concatenate(v) for k, v in Z0.items()}


# --------------------------------------------------------------------------------------------- the census
def load_census():
    g, cls, nP, nM, K, bad = [], [], [], [], [], 0
    with open(STEFF) as fh:
        for line in fh:
            a = line.split()
            if len(a) != 5:
                bad += 1
                continue
            g.append(int(a[0]))
            cls.append(a[1])
            nP.append(int(a[2]))
            nM.append(int(a[3]))
            K.append(int(a[4]))
    check("st_eff.txt: every line has the five fields g class nP nM K", bad == 0)
    g, cls, nP, nM = np.array(g, np.int64), np.array(cls), np.array(nP), np.array(nM)
    check("st_eff.txt: no stratum twice", np.unique(g).size == g.size)
    check("st_eff.txt: every line has a confirmed friendly rival for T>S or T<S", ((nP > 0) | (nM > 0)).all())
    cen = {(c, s): np.sort(g[(cls == c) & ((nP > 0) if s > 0 else (nM > 0))]) for c, s in BLOCKS}
    return dict(g=g, cls=cls, K=sorted(set(K)), cen=cen, n=g.size, bycls=collections.Counter(cls.tolist()))


# --------------------------------------------------------------------------------------------- the probe
ERR_RX = re.compile(r"^stripe (\d+): DMIN=(\d+) want=([CAZ]) sign=(-?1) -> (\d+) strata, max free coords (-?\d+)$")


def load_probe():
    """the strata of each (class, sign) block, with each stripe checked against its task log"""
    pr, nlines, nfiles, dmin, named, counted = {}, 0, 0, set(), [], []
    for b, (c, s) in enumerate(BLOCKS):
        parts = []
        for k in range(NSTRIPE):
            with open(os.path.join(CA1, "dim_%s_%d_%d.txt" % (c, s, k))) as fh:
                rows = [l.split() for l in fh if l.strip()]
            with open(os.path.join(CA1, "dim_%d_%d.err" % (DIMJOB, NSTRIPE * b + k))) as fh:
                err = [l.strip() for l in fh if l.strip()]
            nfiles += 2
            m = ERR_RX.match(err[-1]) if err else None
            named.append(m is not None and (int(m.group(1)), m.group(3), int(m.group(4))) == (k, c, s))
            dmin.add(int(m.group(2)) if m else -1)
            g = np.array([int(r[0]) for r in rows], np.int64)
            nf = np.array([int(r[1]) for r in rows], np.int64)
            counted.append(all(len(r) == 3 and r[2] == c for r in rows) and (nf == nfree(g)).all()
                           and (g % NSTRIPE == k).all() and m is not None and int(m.group(5)) == g.size
                           and int(m.group(6)) == (nf.max() if g.size else -1))
            nlines += g.size
            parts.append(g)
        pr[c, s] = np.sort(np.concatenate(parts))
    check("probe: the task log of each of the %d stripes names its block and stripe" % len(named), all(named))
    check("probe: each stripe holds lines 'g dimension class' of its block and stripe, as many as its task log says, "
          "with the largest dimension it says", all(counted))
    check("probe: every task ran with DMIN = 5", dmin == {5})
    for c, s in BLOCKS:
        check("probe %s_%d: no stratum twice, and every one of dimension >= 5" % (c, s),
              np.unique(pr[c, s]).size == pr[c, s].size and (nfree(pr[c, s]) >= 5).all())
    return pr, nlines, nfiles


# --------------------------------------------------------------------------------------------- the table
def main():
    C = load_census()
    cen = C["cen"]
    P, nprobe, nfiles = load_probe()
    cls_all = np.frombuffer(gzip.decompress(open(CLASSES, "rb").read()), np.uint8)
    check("census_support_class: one byte for each of the 3^16 strata", cls_all.size == N3)
    code = np.array([CODE[x] for x in C["cls"]], np.uint8)
    check("the class the census read from the weights is the class by the support, at all %s census lines"
          % num(C["n"]), (cls_all[C["g"]] == code).all())
    for c, s in BLOCKS:
        check("probe %s_%d: the class the probe read from the weights is the class by the support, at all its strata"
              % (c, s), (cls_all[P[c, s]] == CODE[c]).all())
    M = np.load(MASKS)
    effCC, effALT, rivP, rivM = (M[k] for k in ("effCC", "effALT", "rivP", "rivM"))
    effLINE = np.load(BINEFF)["effLINE"]
    BIN = {"W": np.flatnonzero(effCC & rivP), "S": np.flatnonzero(effCC & rivM),
           "N": np.flatnonzero(effALT & rivP), "E": np.flatnonzero(effALT & rivM)}

    # panel a: the census and the probe together; dimensions, dimension 0, binary-invisible
    A, D, U, INV = [], [], {}, {}
    for name, keys in LOCI:
        ce = np.unique(np.concatenate([cen[k] for k in keys]))
        pr = np.unique(np.concatenate([P[k] for k in keys]))
        u = U[name] = np.union1d(ce, pr)
        nf = nfree(u)
        h, lo, hi = hist(nf)
        row = dict(name=name, n=u.size, lo=lo, hi=hi, d0=int(h[0]), dag=name in ("S", "ray T<S"), inv=None)
        if name in BIN:
            INV[name] = {"census": (invisible(ce, BIN[name]), ce.size),
                         "census + probe": (invisible(u, BIN[name]), u.size)}
            row["inv"], row["base"] = INV[name]["census" if name == "S" else "census + probe"]   # S: the census only
            check("%s: the strata of dimension 0 are the %d binary friendly rivals (m2_masks.npz), as sets"
                  % (name, BIN[name].size), np.array_equal(np.sort(pinned(u[nf == 0])[1]), BIN[name]))
        A.append(row)
        D.append((name, ce, pr, u, h))
    H, Z0 = class_dims(cls_all)
    for name, c, beff in [("T=S, below", [1], effCC), ("T=S, above", [2], effALT), ("(1,0)", [1, 2, 3], effLINE)]:
        h = H[c].sum(0)
        nz = np.flatnonzero(h)
        A.append(dict(name=name, n=int(h.sum()), lo=int(nz[0]), hi=int(nz[-1]), d0=int(h[0]), dag=False, inv=None))
        D.append((name, None, None, None, h))
        check("%s: the strata of dimension 0 are the %d binary strategies efficient there, as sets" % (name, beff.sum()),
              np.array_equal(np.sort(pinned(np.concatenate([Z0[k] for k in c]))[1]), np.flatnonzero(beff)))
    check("the C/D mirror maps N (census and probe) onto E (census and probe), as sets",
          np.array_equal(mirror(U["N"]), U["E"]))
    check("the C/D mirror maps the binary friendly rivals of N onto those of E",
          np.array_equal(np.sort(pinned(mirror(code_of(BIN["N"])))[1]), BIN["E"]))

    # panel b: the census against the probe, strata of dimension >= 5, block by block
    B = []
    for name, k in BROWS:
        ce = cen[k]
        c5 = ce[nfree(ce) >= 5]
        new = np.setdiff1d(P[k], c5)
        B.append(dict(name=name, k=k, c5=c5.size, p5=P[k].size, deficit=new.size, missed=np.setdiff1d(c5, P[k]).size,
                      total=ce.size, share=new.size / ce.size, pshare=new.size / P[k].size))

    # ------------------------------------------------------------------------------------- the output
    fa = "  %-12s %-20s %-12s %-16s %-18s %s"
    fb = "  %-32s %-18s %-17s %-9s %-14s %s"
    L = ["sitable1.py -- SI Table 1 of PRFull (the friendly-rival strata of the memory-two cube), computed from the "
         "census output.", "Written: Figures/sitable1.out.", "",
         "a  by locus: the census completed by the dimension probe",
         fa % ("locus", "strata", "dimensions", "of dimension 0", "binary-invisible", "share invisible")]
    for r in A:
        inv = "--" if r["inv"] is None else num(r["inv"]) + (" dagger" if r["name"] == "S" else "")
        L.append(fa % (r["name"], num(r["n"]) + (" dagger" if r["dag"] else ""), "%d to %d" % (r["lo"], r["hi"]),
                       num(r["d0"]), inv, "--" if r["inv"] is None else pct2(r["inv"] / r["base"])))
    s = A[1]
    L += ["  dagger: a numerical census count.  The binary-invisible count of S is taken over the census's %s strata "
          "alone: %s of them." % (num(s["base"]), pct2(s["inv"] / s["base"])), "",
          "b  the census against the dimension probe: strata with a friendly rival, of dimension at least 5",
          fb % ("locus", "census, dim >= 5", "probe, dim >= 5", "deficit", "census total", "share")]
    for r in B:
        L.append(fb % (r["name"], num(r["c5"]), num(r["p5"]), num(r["deficit"]), num(r["total"]), pct_b(r["share"])))
    L += ["  the deficit as a share of the probe's count: " + ", ".join(
        "%s %s" % (r["name"].split()[1], pct2(r["pshare"])) for r in B[:4]), ""]

    L += ["#" * 110, "HOW EACH CELL IS COMPUTED",
          "  a, strata             the strata of the locus in the census (st_eff.txt: the class, and nP > 0 for T>S or "
          "nM > 0 for T<S),",
          "                        together with those of the probe's blocks of the locus; the rays take the classes C, "
          "A and Z",
          "  a, dimensions         the smallest and the largest number of free coordinates among these strata",
          "  a, of dimension 0     the strata with no free coordinate, each a binary strategy",
          "  a, binary-invisible   the strata whose closed cell holds no binary friendly rival of the wedge "
          "(m2_masks.npz;",
          "                        census.py's visible_flags), in S over the census's strata only; share = "
          "binary-invisible / those strata",
          "  a, T=S and (1,0)      the efficient strata by the support (census_support_class.u8.gz): class C below the "
          "switch line,",
          "                        class A above it, and C, A and Z at (1,0)",
          "  b                     per (class, sign) block of the probe: the census's strata of dimension >= 5, the "
          "probe's strata,",
          "                        the probe's strata not in the census (the deficit), and all the census's strata of "
          "the block;",
          "                        share = deficit / census total", ""]

    L += ["#" * 110, "INPUTS",
          "  PartnersRivalsM2f/check/st_eff.txt: %s lines (%s), K = %s draws at every stratum" % (
              num(C["n"]), ", ".join("%s %s" % (k, num(C["bycls"][k])) for k in "CAZ"),
              "/".join(str(x) for x in C["K"])),
          "  %s: the dimension probe (array %d), %d files dim_<class>_<sign>_<stripe>.txt and dim_%d_<task>.err, "
          "%s lines" % (CA1_NAME, DIMJOB, nfiles, DIMJOB, num(nprobe)),
          "  PRFull/release/census_class/census_support_class.u8.gz: %s strata; by the support none %s, C %s, A %s, "
          "Z %s" % (num(cls_all.size), *(num(int((cls_all == v).sum())) for v in range(4))),
          "  PartnersRivals/Figures/m2_masks.npz: binary friendly rivals W %d, S %d, N %d, E %d; effCC %d, effALT %d"
          % (*(BIN[k].size for k in "WSNE"), effCC.sum(), effALT.sum()),
          "  PartnersRivalsM2f/check/binary_eff.npz: effLINE %d" % effLINE.sum(), ""]

    L += ["#" * 110, "PER LOCUS: the census alone, the probe, what the probe added, and the dimension histogram of panel a"]
    for name, ce, pr, u, h in D:
        if ce is None:
            L.append("  %-12s %s strata; dimensions %s" % (name, num(int(h.sum())), hist_txt(h)))
            continue
        add = np.setdiff1d(pr, ce)
        L.append("  %-12s census %s, probe %s, added by the probe %s, together %s; dimensions %s" % (
            name, num(ce.size), num(pr.size), num(add.size), num(u.size), hist_txt(h)))
        if 0 < add.size <= 8:
            L.append("  %-12s   added: %s" % ("", "  ".join(pat(g) for g in add)))
    L += ["", "  binary-invisible strata, both counts (the table takes census and probe, in S the census alone):"]
    for name in "WSNE":
        for lab, (n, m) in INV[name].items():
            L.append("    %-2s %-15s %10s of %10s  %s" % (name, lab, num(n), num(m), pct2(n / m)))
    L.append("")

    L += ["#" * 110, "PANEL b IN FULL"]
    for r in B:
        L.append("  %-32s block %s_%d: census %s strata, %s of dimension >= 5; probe %s; probe not in census %s (the "
                 "deficit); census not in probe %s; deficit / census total %.4f%%, / probe %.4f%%" % (
                     r["name"], r["k"][0], r["k"][1], num(r["total"]), num(r["c5"]), num(r["p5"]), num(r["deficit"]),
                     num(r["missed"]), 100 * r["share"], 100 * r["pshare"]))
    for c, s in BLOCKS:
        h, lo, hi = hist(nfree(P[c, s]))
        L.append("  probe %s_%d: %s strata of dimension >= 5, largest %d; %s" % (c, s, num(P[c, s].size), hi, hist_txt(h)))
    L.append("")

    bad = [w for w, ok in CHECKS if not ok]
    L += ["#" * 110, "CHECKS"] + ["  %-6s %s" % ("ok" if ok else "FAILED", w) for w, ok in CHECKS]
    L.append("checks failed: %d of %d" % (len(bad), len(CHECKS)))
    with open(OUT, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L[:L.index("#" * 110)]))
    for w in bad:
        print("FAILED: " + w)
    print("checks failed: %d of %d" % (len(bad), len(CHECKS)))
    print("wrote " + os.path.relpath(OUT, os.path.dirname(HERE)).replace(os.sep, "/"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
