#!/usr/bin/env python3
"""
fig1data.py -- the numbers of Figure 1 (friendly rivals over the plane of games, stochastic memory one and memory two),
each read from its source and, where possible, recomputed here and asserted equal.

    python3 fig1data.py        prints the table with its provenance and the checks; fig1.py imports LOCI, M1, M2

The nine loci of the plane of games (R,S,T,P) = (1,u,1+v,0), cut by the switch line u+v = 1 and the line T = S (u-v = 1):

    W     wedge u+v<1, T>S          rayP  switch line, T>S (u<1)       tsLO  T=S below the switch line (u<1)
    S     wedge u+v<1, T<S          rayM  switch line, T<S (u>1)       tsHI  T=S above the switch line (u>1)
    N     wedge u+v>1, T>S                                             pt    the corner (1,0), R=S=T
    E     wedge u+v>1, T<S

Per locus: what efficiency requires; memory one: binary friendly rivals (of 16) and the dimension of the friendly-rival
set in [0,1]^4; memory two: binary friendly rivals (of 65536) and the dimension in [0,1]^16.

SOURCES (all read-only)
  M2 dim, M2 binary, M1 dim   Table 1 of PartnersRivalsM2f/Text/tables.tex, parsed below (its last column is the
                              memory-one continuum answer).  The memory-two dimensions there came from the census;
                              since 2026-09-27 they are proved exactly, with the same values: lower bounds by the
                              certified boxes of PRFull/check/cert/results (W 7, N and E 8, S 15; the N box lies on
                              the ray T>S, the S boxes on the ray T<S), upper bounds by the exclusion certificates of
                              PRFull/check/ubound; on T=S and at (1,0) by the closed forms of efficiency.
  M2 binary, recomputed       PartnersRivals1/Figures/m2_masks.npz (effCC, effALT, rivP, rivM) and
                              PartnersRivalsM2f/check/binary_eff.npz (effLINE = w_DD = 0): all nine counts.
  M1 binary                   PartnersRivals1/Revised2/ms.tex, SI "Memory one, for contrast" (0, 2, 0, 0 on W, S, E, N)
                              and SI "Where the classes live in the plane of games", "The two lines" (none on the ray T>S,
                              CCCC, CCCD, CDCC on the ray T<S; on T=S the efficient ones, three below and none above).
                              The corner (1,0) is not stated there as a count: it is computed here (the five strategies
                              efficient on the switch line, all rivals since T = S), and it equals the 5 of
                              PartnersRivalsM1f/Text/friendly.md, section 2, boundary-locus table.
  M1 binary, recomputed       exactly, below: the 16 vertices, the eps -> 0 limit of the pair chain by the Markov-chain
                              tree theorem with polynomial-in-eps coefficients in exact Fractions, co-players the 16
                              vertices (enough for a memory-one focal strategy: reference SI "Memory one, for contrast").
  M1 dim, cross-checked       PartnersRivalsM1f/Text/friendly.md section 2 (wedge table and boundary-locus table) and
                              PartnersRivalsM1f/Text/ms.tex Prop. 3 (empty for T>S and on open E; dim 3 on S); spot
                              checks below at exact rational points of the claimed sets.
"""
import json
import os
import re
import sys
from fractions import Fraction as F
from itertools import product

sys.dont_write_bytecode = True
import numpy as np

CL1 = "/Users/martin/Documents/CL1"
TABLES = os.path.join(CL1, "PartnersRivalsM2f/Text/tables.tex")
REFMS = os.path.join(CL1, "PartnersRivals1/Revised2/ms.tex")
M1FR = os.path.join(CL1, "PartnersRivalsM1f/Text/friendly.md")
M1MS = os.path.join(CL1, "PartnersRivalsM1f/Text/ms.tex")
MASKS = os.path.join(CL1, "PartnersRivals1/Figures/m2_masks.npz")
# The values this script takes from the manuscripts above (Table 1 of tables.tex) are deposited in
# fig1_reference_values.json, which a normal run reads; the manuscripts are read only by --verify-sources (on the
# author's machine), which checks that file and the sentences below against them.  The code and data release
# contains no manuscript.
REFVALS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig1_reference_values.json")
BEFF = os.path.join(CL1, "PartnersRivalsM2f/check/binary_eff.npz")

LOCI = ["W", "S", "N", "E", "rayP", "rayM", "tsLO", "tsHI", "pt"]
EFF = {"W": "CC", "S": "CC", "N": "ALT", "E": "ALT", "rayP": "Z", "rayM": "Z", "tsLO": "CC", "tsHI": "ALT", "pt": "Z"}
RIV = {"W": "P", "S": "M", "N": "P", "E": "M", "rayP": "P", "rayM": "M", "tsLO": None, "tsHI": None, "pt": None}   # None: T = S


# ============================================================== memory two: Table 1 of tables.tex
def table1_from_source():
    """the nine rows of Table 1 (tab:loci) of tables.tex: {locus: (efficiency, dim, binary, memory-one answer)}"""
    txt = open(TABLES).read()
    t1 = txt[txt.index(r"\label{tab:loci}"):]
    body = t1[t1.index(r"\begin{tabular}"):t1.index(r"\end{tabular}")]
    rows = [l for l in body.splitlines() if "&" in l and "locus" not in l]
    assert len(rows) == 9, rows
    pat = [r"wedge \$W", r"wedge \$S", r"wedge \$N", r"wedge \$E", r"ray \$u\+v=1\$, \$u<1\$ \(\$T>S\$\)",
           r"ray \$u\+v=1\$, \$u>1\$ \(\$T<S\$\)", r"line \$T=S\$, \$u\+v<1\$", r"line \$T=S\$, \$u\+v>1\$", r"point \$\(u,v\)=\(1,0\)\$"]
    out = {}
    for key, p, row in zip(LOCI, pat, rows):
        assert re.match(r"\s*" + p, row), (key, row)
        c = [x.strip() for x in row.rstrip().rstrip("\\").split("&")]
        dim = int(re.fullmatch(r"\$\\mathbf\{(\d+)\}\$", c[2]).group(1))
        binary = int(re.fullmatch(r"\$(\d+)\$", c[4].replace(r"\,", "")).group(1))
        m1 = c[5]
        m1dim = None if m1 == "empty" else int(re.match(r"dim \$(\d+)\$", m1).group(1))
        out[key] = (c[1], dim, binary, m1dim, m1)
    return out


def table1():
    """the nine rows of Table 1, from fig1_reference_values.json (written from tables.tex by --write-reference-values)"""
    d = json.load(open(REFVALS))["table1"]
    return {k: tuple(d[k]) for k in LOCI}


T1 = table1()
M2 = {k: {"dim": T1[k][1], "bin": T1[k][2]} for k in LOCI}


def m2_recount():
    """the nine binary memory-two counts from the masks"""
    M = np.load(MASKS, allow_pickle=False)
    B = np.load(BEFF, allow_pickle=False)
    eff = {"CC": M["effCC"], "ALT": M["effALT"], "Z": B["effLINE"]}
    W = B["W"]                                   # self-play limits (w_CC, w_CD, w_DC, w_DD), same genome indexing?
    assert ((W[:, 0] > 1 - 1e-12) == eff["CC"]).all(), "binary_eff.npz and m2_masks.npz disagree on the indexing"
    riv = {"P": M["rivP"], "M": M["rivM"], None: np.ones(65536, bool)}
    return {k: int((eff[EFF[k]] & riv[RIV[k]]).sum()) for k in LOCI}


# ============================================================== memory one: exact eps -> 0 limits at any rational p, q
# States 0 CC, 1 CD, 2 DC, 3 DD, focal action first; p = (p1, p2, p3, p4) the probability of C after each.
# Error p -> eps + (1-2eps) p on both players.  A polynomial in eps is a tuple of Fractions, lowest order first.
SWAP = (0, 2, 1, 3)


def pmul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def coop(x):
    """(C, D) probabilities of a clamped move, as polynomials in eps"""
    x = F(x)
    return {1: [x, 1 - 2 * x], 0: [1 - x, 2 * x - 1]}


# the 16 in-trees on four states for each root: parent maps of the other three states that reach the root
TREES = {}
for r in range(4):
    others = [s for s in range(4) if s != r]
    TREES[r] = []
    for par in product(range(4), repeat=3):
        f = dict(zip(others, par))
        if any(f[s] == s for s in others):
            continue
        ok = True
        for s in others:
            seen, x = set(), s
            while x != r:
                if x in seen:
                    ok = False
                    break
                seen.add(x)
                x = f[x]
            if not ok:
                break
        if ok:
            TREES[r].append(f)
assert all(len(TREES[r]) == 16 for r in range(4))


def limit(p, q):
    """w = lim_{eps->0} of the stationary distribution of the pair (p focal, q co-player), exact Fractions"""
    P = {}
    for s in range(4):
        a, b = coop(p[s]), coop(q[SWAP[s]])
        for t in range(4):
            P[s, t] = pmul(a[1 if t in (0, 1) else 0], b[1 if t in (0, 2) else 0])
    N = []
    for r in range(4):
        acc = [F(0)]
        for f in TREES[r]:
            term = [F(1)]
            for s, t in f.items():
                term = pmul(term, P[s, t])
            acc = padd(acc, term)
        N.append(acc)
    tot = [F(0)]
    for n in N:
        tot = padd(tot, n)
    m = next(i for i, c in enumerate(tot) if c != 0)
    assert tot[m] > 0
    return tuple((N[r][m] if m < len(N[r]) else F(0)) / tot[m] for r in range(4))


VERT = [tuple(int(c) for c in "{:04b}".format(g)) for g in range(16)]        # "1011" = C at CC, D at CD, C at DC, C at DD


def name(p):
    return "".join("C" if x == 1 else "D" for x in p)


def classify(p):
    """the three efficiency classes and the two rival classes of a memory-one strategy p (rational coordinates)"""
    w = limit(p, p)
    assert w[1] == w[2], ("self-play role symmetry fails", p, w)
    eff = {"CC": w[0] == 1, "ALT": w[1] == F(1, 2) and w[2] == F(1, 2), "Z": w[3] == 0}
    marg = [ (lambda x: x[2] - x[1])(limit(p, q)) for q in VERT ]        # w_DC - w_CD against each vertex co-player
    riv = {"P": min(marg) >= 0, "M": max(marg) <= 0, None: True}
    return eff, riv, w


def m1_binary():
    """the nine binary memory-one counts, and the members"""
    C = {v: classify(v) for v in VERT}
    mem = {k: [name(v) for v in VERT if C[v][0][EFF[k]] and C[v][1][RIV[k]]] for k in LOCI}
    return {k: len(mem[k]) for k in LOCI}, mem, C


# the memory-one closed forms of PartnersRivalsM1f (ms.tex Props. 1, 2; friendly.md section 1), for the spot checks
def effCC(p):
    p1, p2, p3, p4 = p
    return p1 == 1 and (p2, p3) not in ((0, 1), (1, 0)) and (p4 > 0 or max(p2, p3) == 1)


def effALT(p):
    p1, p2, p3, p4 = p
    return (p2, p3) in ((0, 1), (1, 0)) and p1 < 1 and p4 > 0 and (p1, p4) != (0, 1)


def rivP(p):
    return p[3] == 0 and p[1] + p[2] <= 1 and p[1] < 1


def rivM(p):
    return p[0] == 1 and p[1] + p[2] >= 1 and p[2] > 0


def m1_spot(seed=1, per=2):
    """exact spot checks of the memory-one continuum.  (i) every one of the 81 support strata of [0,1]^4 (each
    coordinate 0, 1 or free), `per` random rational interior draws each: the computed efficiency and rivalry against the
    closed forms, and friendly rivals in W, N, E and on the ray T>S (claimed empty); (ii) open 3-dimensional and
    2-dimensional pieces: random points of the facet p1 = 1 near (1, 7/10, 7/10, 1/2), claimed friendly rivals on S, on
    the ray T<S, on T=S below the switch line and at (1,0); and of the sheet (p1, 0, 1, p4), claimed friendly rivals on
    T=S above the switch line, with the off-corner control (p1, 1/10, 1, p4) not efficient there."""
    rng = np.random.default_rng(seed)
    rnd = lambda lo=1, hi=99: F(int(rng.integers(lo, hi + 1)), 100)
    bad, hits, n = 0, {k: 0 for k in ("W", "N", "E", "rayP", "S")}, 0
    for st in product((0, 1, 2), repeat=4):
        for _ in range(per):
            p = tuple(F(x) if x < 2 else rnd() for x in st)
            eff, riv, w = classify(p)
            n += 1
            bad += (eff["CC"] != effCC(p)) + (eff["ALT"] != effALT(p)) + (riv["P"] != rivP(p)) + (riv["M"] != rivM(p))
            for k in hits:
                hits[k] += bool(eff[EFF[k]] and riv[RIV[k]])
    box3 = [(F(1), rnd(55, 85), rnd(55, 85), rnd(20, 80)) for _ in range(12)]
    ok3 = all((lambda e, r: e["CC"] and e["Z"] and r["M"] and not r["P"])(*classify(p)[:2]) for p in box3)
    box2 = [(rnd(20, 80), F(0), F(1), rnd(20, 80)) for _ in range(8)] + [(rnd(20, 80), F(1), F(0), rnd(20, 80)) for _ in range(4)]
    ok2 = all(classify(p)[0]["ALT"] for p in box2)
    ctrl = [(rnd(20, 80), F(1, 10), F(1), rnd(20, 80)) for _ in range(6)]
    okc = not any(classify(p)[0]["ALT"] for p in ctrl)
    return dict(n=n, mismatches=bad, hits=hits, box3=ok3, box2=ok2, ctrl=okc)


# ============================================================== memory one: the values, with their sources
# binary counts, as stated in the reference paper (PartnersRivals1/Revised2/ms.tex); the corner computed (see above)
M1BIN_REF = {"W": 0, "S": 2, "N": 0, "E": 0, "rayP": 0, "rayM": 3, "tsLO": 3, "tsHI": 0}
# REF_SENTENCES: sentences checked verbatim in the reference paper; kept in fig1_source_quotes.json (author only)
# continuum dimensions in [0,1]^4 (None = empty): PartnersRivalsM1f/Text/friendly.md section 2
M1DIM_FR = {"W": None, "S": 3, "N": None, "E": None, "rayP": None, "rayM": 3, "tsLO": 3, "tsHI": 2, "pt": 3}
# M1_SENTENCES: the same for the memory-one sources (PartnersRivalsM1f/Text/ms.tex and friendly.md), in the same file

M1 = {k: {"dim": M1DIM_FR[k]} for k in LOCI}      # "bin" filled in by check() below, from the exact computation


def verify_sources(say=print):
    """author's machine only: the deposited reference values and the quoted sentences against the manuscripts"""
    assert table1_from_source() == T1, "fig1_reference_values.json differs from Table 1 of tables.tex"
    Q = json.load(open(os.path.join(os.path.dirname(REFVALS), "fig1_source_quotes.json")))
    ref = open(REFMS).read()
    for s in Q["reference"]:
        assert s in ref, "reference no longer says: " + s
    for f, s in Q["memory_one"]:
        assert s in open(os.path.join(CL1, f)).read(), "%s no longer says: %s" % (os.path.basename(f), s)
    say("sources verified: fig1_reference_values.json equals Table 1 of tables.tex, and all %d quoted sentences are "
        "present" % (len(Q["reference"]) + len(Q["memory_one"])))


def check(verbose=True):
    say = print if verbose else (lambda *a, **k: None)
    # memory one, binary: exact recomputation against the reference
    cnt, mem, C = m1_binary()
    for k, v in M1BIN_REF.items():
        assert cnt[k] == v, (k, cnt[k], v)
    for k in LOCI:
        M1[k]["bin"] = cnt[k]
    assert sorted(mem["rayM"]) == sorted(["CCCC", "CCCD", "CDCC"]) and sorted(mem["tsLO"]) == sorted(["CCCC", "CCCD", "CDDC"])
    assert sorted(mem["pt"]) == sorted(["CCCC", "CCCD", "CDDC", "CCDC", "CDCC"])
    # memory one, continuum: friendly.md against the memory-one column of tables.tex
    for k in LOCI:
        assert T1[k][3] == M1DIM_FR[k], (k, T1[k], M1DIM_FR[k])
    assert T1["tsHI"][4] == "dim $2$, no vertex"
    sp = m1_spot()
    assert sp["mismatches"] == 0 and not any(sp["hits"][k] for k in ("W", "N", "E", "rayP")) and sp["box3"] and sp["box2"] and sp["ctrl"], sp
    # memory two, binary: masks against Table 1
    rc = m2_recount()
    for k in LOCI:
        assert rc[k] == M2[k]["bin"], (k, rc[k], M2[k]["bin"])
    say("%-5s %-12s | %-18s | %-18s" % ("locus", "efficiency", "M1: binary / dim", "M2: binary / dim"))
    for k in LOCI:
        say("%-5s %-12s | %5d  %-11s | %6d  %-10s" % (k, T1[k][0].replace("$", ""), M1[k]["bin"],
                                                    "empty" if M1[k]["dim"] is None else "dim %d" % M1[k]["dim"],
                                                    M2[k]["bin"], "dim %d" % M2[k]["dim"]))
    say("memory-one binary friendly rivals: " + "; ".join("%s %s" % (k, ",".join(mem[k]) or "-") for k in LOCI))
    say("memory-one spot checks: %d strategies over the 81 strata, %d closed-form mismatches, friendly rivals W %d N %d E %d "
        "rayP %d (S %d); facet box dim-3 test %s; alternating sheet dim-2 test %s; off-corner control %s"
        % (sp["n"], sp["mismatches"], sp["hits"]["W"], sp["hits"]["N"], sp["hits"]["E"], sp["hits"]["rayP"], sp["hits"]["S"],
           sp["box3"], sp["box2"], sp["ctrl"]))
    say("memory-two binary counts recomputed from the masks: " + ", ".join("%s %d" % (k, rc[k]) for k in LOCI) + "  (= Table 1)")
    return mem, sp


if __name__ == "__main__":
    if "--write-reference-values" in sys.argv[1:]:       # author's machine: from tables.tex
        t = table1_from_source()
        json.dump({"note": "Table 1 (tab:loci) of PartnersRivalsM2f/Text/tables.tex, parsed by fig1data.py: per locus "
                           "(efficiency, dimension, binary count, memory-one dimension or null, memory-one entry)",
                   "table1": {k: list(t[k]) for k in LOCI}}, open(REFVALS, "w"), indent=1)
        print("wrote " + REFVALS)
    elif "--verify-sources" in sys.argv[1:]:
        verify_sources()
    else:
        check()
