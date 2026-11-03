"""
The exact eps -> 0 layer over the whole plane, vectorised.

For any game (u,v) this returns the four boolean masks over the 65536
strategies that the vocabulary of section 2 needs:

    eff   attains Emax = max(1, (1+u+v)/2) in self play
    riv   Riv+ if T > S, Riv- if T < S, everything on T = S
    ne    an eps -> 0 Nash equilibrium at that game
    fr    eff & riv, the friendly rivals    (= partner & rival, see PROJECT.md)

and partner = eff & ne.  Everything comes from CL/BinM2NE1's exact census and
Donation's ca17-certified rival mask; nothing here is simulated.

Facets are integer triples (a,b,c) meaning a + b u + c v <= 0.  The tie strata
(dim 1 and dim 0) have measure zero in the plane, so a raster point never lands
on one; they are evaluated with a tolerance and counted separately so that the
approximation is visible rather than silent.
"""
import csv, os
import numpy as np

# The eps -> 0 LIMIT Nash oracle, and it must be CL/BinM2NE1's -- that project
# owns the equilibrium layer for both error regimes.  This file is the exact
# census (Markov chain tree theorem in integer arithmetic, no floating point,
# no extrapolation); verified against BinM2NE1's own MANIFEST.sha256 on
# 2026-09-04, ec98f6d2f91beb49a03632e0d18b5de81e6bf6f6cb0515b6d279a6e23e9dce01.
# The FIXED-eps layer is a different instrument and lives in
# Figures/nashe4.py, on that project's exactcoef.py; the two are not
# interchangeable (the 16-flip reduction is valid at eps > 0 and invalid in
# the limit).
CENSUS = os.environ.get(
    "BINM2NE1_CENSUS",
    "/Users/martin/Documents/CL/BinM2NE1/Final/Github/exact/ca1_regions.csv")
RIVAL = "/Users/martin/Documents/CL1/Donation/Data/rival_mask_ca17.npy"
NC = 65536
EIGHT = [19079, 19087, 20103, 20111, 23175, 23183, 24199, 24207]
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "Data", "exactlayer.npz")


def bar_table():
    b = np.empty(NC, dtype=np.int32)
    for s in range(NC):
        r = 0
        for k in range(16):
            if s >> k & 1:
                r |= 1 << (15 - k)
        b[s] = 65535 - r
    return b


def build(path=CENSUS):
    dim = np.zeros(NC, dtype=np.int8)
    W = np.zeros((NC, 4), dtype=np.int64)
    den = np.zeros(NC, dtype=np.int64)
    fa, fb, fc, cnt = [], [], [], np.zeros(NC, dtype=np.int32)
    ln = np.zeros((NC, 3), dtype=np.int64)          # dim 1: A u + C v = E
    seg = np.zeros((NC, 2), dtype=float)            # its endpoints, nan = open
    seg[:] = np.nan
    pnt = np.zeros((NC, 2), dtype=float)            # dim 0: the single game
    pnt[:] = np.nan
    for r in csv.DictReader(open(path)):
        g = int(r["genotype"]); dim[g] = int(r["dim"])
        W[g] = [int(r["w%d" % i]) for i in range(4)]
        den[g] = int(r["wden"])
        if r["facets"]:
            for t in r["facets"].split(";"):
                a, b, c = (int(x) for x in t.split(","))
                fa.append(a); fb.append(b); fc.append(c)
            cnt[g] = len(r["facets"].split(";"))
        if r["line"]:
            ln[g] = [int(x) for x in r["line"].split(",")]
            from fractions import Fraction as F
            seg[g, 0] = float(F(r["x0"])) if r["x0"] else np.nan
            seg[g, 1] = float(F(r["x1"])) if r["x1"] else np.nan
        if r["pu"]:
            from fractions import Fraction as F2
            pnt[g] = [float(F2(r["pu"])), float(F2(r["pv"]))]
    off = np.zeros(NC + 1, dtype=np.int64)
    off[1:] = np.cumsum(cnt)
    # facets are stored genotype by genotype in file order, which is genotype
    # order, so the offsets line up with the concatenation
    return dict(dim=dim, W=W, den=den, ln=ln, seg=seg, pnt=pnt, cnt=cnt, off=off,
                fa=np.array(fa, dtype=np.int64),
                fb=np.array(fb, dtype=np.int64),
                fc=np.array(fc, dtype=np.int64))


def load():
    p = os.path.abspath(CACHE)
    if os.path.exists(p):
        z = np.load(p)
        d = {k: z[k] for k in z.files}
    else:
        d = build()
        os.makedirs(os.path.dirname(p), exist_ok=True)
        np.savez_compressed(p, **d)
    d["rivp"] = np.load(RIVAL)
    d["bar"] = bar_table()
    rivm = np.zeros(NC, dtype=bool)
    rivm[d["bar"][np.nonzero(d["rivp"])[0]]] = True
    d["rivm"] = rivm
    w0, w1, w3 = d["W"][:, 0], d["W"][:, 1], d["W"][:, 3]
    d["eff_lo"] = (w1 == 0) & (w3 == 0)      # all-CC self play,   7639
    d["eff_hi"] = (w0 == 0) & (w3 == 0)      # perfect alternation, 3072
    d["eff_on"] = (w3 == 0)                  # on the switch line, 14757
    return d


def masks(d, u, v, tol=0.0):
    """eff, riv, ne at the game (u,v)"""
    if u + v < -tol + 1:
        eff = d["eff_lo"]
    elif u + v > tol + 1:
        eff = d["eff_hi"]
    else:
        eff = d["eff_on"]

    if u - v < -tol + 1:
        riv = d["rivp"]
    elif u - v > tol + 1:
        riv = d["rivm"]
    else:
        riv = np.ones(NC, dtype=bool)

    val = d["fa"] + d["fb"] * u + d["fc"] * v
    ok = (val <= tol).astype(np.int32)
    # a dim-2 resident is an equilibrium iff every one of its facets holds
    tot = np.add.reduceat(ok, d["off"][:-1]) if len(ok) else np.zeros(NC)
    tot = np.where(d["cnt"] > 0, tot, 0)
    ne = (d["dim"] == 2) & (tot == d["cnt"])
    # the dim-1 and dim-0 strata are measure zero; include them if the point
    # happens to lie on them
    m1 = np.nonzero(d["dim"] == 1)[0]
    if len(m1):
        A, C, E = d["ln"][m1, 0], d["ln"][m1, 1], d["ln"][m1, 2]
        x = np.where(C != 0, u, v)      # geom_ca1's segment parameter
        x0, x1 = d["seg"][m1, 0], d["seg"][m1, 1]
        hit = (np.abs(A * u + C * v - E) <= tol) \
              & (np.isnan(x0) | (x >= x0 - tol)) \
              & (np.isnan(x1) | (x <= x1 + tol))
        ne[m1[hit]] = True
    m0 = np.nonzero(d["dim"] == 0)[0]
    if len(m0):
        hit0 = (np.abs(d["pnt"][m0, 0] - u) <= tol) \
               & (np.abs(d["pnt"][m0, 1] - v) <= tol)
        ne[m0[hit0]] = True             # 52 games only; never on a raster
    return eff, riv, ne


def counts(d, u, v):
    eff, riv, ne = masks(d, u, v)
    fr = eff & riv
    return dict(eff=int(eff.sum()), riv=int(riv.sum()), ne=int(ne.sum()),
                partner=int((eff & ne).sum()), fr=int(fr.sum()),
                fr_in_ne=bool((fr & ~ne).sum() == 0))


if __name__ == "__main__":
    d = load()
    print("census loaded: %d facets over %d full-dimensional residents"
          % (len(d["fa"]), int((d["dim"] == 2).sum())))
    for (u, v) in [(-2, 2), (-0.5, 0.5), (0, 0), (-20, -20), (3, 0), (1, 2),
                   (2, -3), (-8, 1)]:
        c = counts(d, u, v)
        print("  (%7.2f,%7.2f)  eff %6d  riv %6d  NE %6d  partner %6d"
              "  FR %5d  FR<=NE %s" %
              (u, v, c["eff"], c["riv"], c["ne"], c["partner"], c["fr"],
               c["fr_in_ne"]))
