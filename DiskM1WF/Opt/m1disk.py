#!/usr/bin/env python3
"""
The seven atoms of (efficient, nash, competitive) for the 16 binary memory-one strategies at each of
the 512 sunflower games, exactly, in the limit eps -> 0 -- the classification layer of this project's
figures, and the thing that turns a WF abundance histogram into an atom share.

    python3 m1disk.py            build and check the cache, print the census   -> ../Data/m1_atoms_512.npz
    python3 m1disk.py rebuild    ignore an existing cache and rebuild it

WHERE THE CLASSIFICATION COMES FROM.  PartnersRivals/FinalFigures/m1atoms.py does the exact work and
is used read-only: for every ordered pair of the 16 it takes the eps -> 0 limit of the four-state
chain in Fractions (Markov chain tree theorem), so the limit payoff is linear in (u, v) with rational
coefficients, and `codes_at(u, v)` then decides the three properties at a rational game with no
threshold anywhere:

    efficient    pi(i, i) = Emax = max(R, (S + T)/2), the largest average any pair of the 16 can earn
    nash         pi(i, i) >= pi(j, i) for all 16 j          (ties allowed)
    competitive  pi(i, j) >= pi(j, i) for all 16 j

The atom code is 4*efficient + 2*nash + competitive, so 110 is a partner and 111 a friendly rival.
The same module checked (its `m2` mode) that each of the 16 has these same three properties against
all 65536 memory-two co-players, on all 45 faces of the plane, so the atom of a memory-one strategy
does not depend on cutting the space down to 16.

THE GAMES ARE THE SIMULATED ONES.  They are read from the kit's own ../Cannon/mw/games.dat -- the
512 rows the Fortran read -- not from a points file, so a classification can never be attached to a
game that was not run.  Every u and v is a multiple of 1/1024 (ca6 snapped them so), which this
module asserts and then uses: the Fraction handed to `codes_at` is exactly the double the simulator
used, so "the atom at this game" means the atom at the game that was simulated.

THE CHECKS (all of them run on every build; the cache is not written if one fails):
  * every u, v is an exact multiple of 1/1024, and u = cS, 1 + v = cT of the same row;
  * the eight atom counts sum to 16 at every game;
  * 101 -- efficient and competitive but not nash -- is empty at all 512 games (an efficient
    competitive strategy is nash), the premise the seven-panel figures rest on;
  * INDEPENDENT ROUTE: each game is placed in the 45-face arrangement m1atoms.py built for the whole
    plane, by the sign vector of the 11 active lines, and its 16 codes must equal that face's codes.
    Two things have to be right for that to pass -- this module's own exact evaluation, and the
    arrangement -- and they were built by different code paths.
  * the patterns of non-empty atoms that occur among the 512 are a subset of the twelve that occur
    on the plane (SI Figure 4's panels).
"""
import os
import sys
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FF = "/Users/martin/Documents/CL1/PartnersRivals/FinalFigures"
FIGDIR = "/Users/martin/Documents/CL1/PartnersRivals/Figures"
GAMES = os.path.join(HERE, os.pardir, "Cannon", "mw", "games.dat")
DUMP = os.path.join(HERE, os.pardir, "Data", "m1_atoms_512.npz")
ARR = os.path.join(FIGDIR, "m1_atoms_k4.npz")            # m1atoms.py's 45-face arrangement of the plane
GRID = 1024                                              # ca6 snapped every game to a multiple of 1/1024
NSTRAT = 16
ORDER = ["000", "100", "010", "001", "110", "011", "111"]
NAME = {"000": "none of the three", "100": "efficient only", "010": "nash only",
        "001": "competitive only", "110": "efficient and nash only",
        "011": "nash and competitive only", "111": "efficient, nash and competitive"}
NICK = {15: "ALLC", 0: "ALLD", 5: "TFT", 9: "WSLS", 1: "Grim"}


def genome(s):
    """the strategy's four answers, own action first: CC, CD, DC, DD"""
    return "".join("C" if (s >> k) & 1 else "D" for k in range(4))


def read_games(path=GAMES):
    """the 512 rows the simulator read: (ipt, gid, tag, u, v) with u, v exact Fractions"""
    out = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        gid, cR, cS, cT = int(f[0]), float(f[2]), float(f[3]), float(f[4])
        u, v, tag = float(f[5]), float(f[6]), f[7]
        assert cR == 1.0, "%s: cR is %r, not 1" % (tag, cR)
        assert cS == u and cT == 1.0 + v, "%s: (cS, cT) = (%r, %r) is not (u, 1+v)" % (tag, cS, cT)
        uf, vf = F(u).limit_denominator(GRID), F(v).limit_denominator(GRID)
        assert float(uf) == u and float(vf) == v and uf.denominator <= GRID \
            and vf.denominator <= GRID, "%s: (u, v) = (%r, %r) is not on the 1/%d grid" % (tag, u, v, GRID)
        out.append((len(out), gid, tag, uf, vf))
    assert len(out) == 512, "%s has %d games, expected 512" % (path, len(out))
    return out


def face_of(u, v, lines, uvc):
    """the face of m1atoms.py's arrangement that holds (u, v): the sign vector of the 11 lines"""
    def sgn(x, y):
        return tuple(1 if a * x + b * y + c > 0 else -1 for a, b, c in lines)
    want = sgn(F(u), F(v))
    hit = [f for f in range(len(uvc)) if sgn(uvc[f][0], uvc[f][1]) == want]
    assert len(hit) == 1, "(%s, %s) lands in %d faces" % (u, v, len(hit))
    return hit[0]


def build():
    sys.path.insert(0, FF)
    import m1atoms                                                    # noqa: E402  (read-only)
    G = read_games()
    A = np.load(ARR)
    lines = [tuple(int(t) for t in r) for r in A["lines"]]
    uvc = [(F(x).limit_denominator(10 ** 12), F(y).limit_denominator(10 ** 12)) for x, y in A["uvc"]]
    codes = np.zeros((512, NSTRAT), np.int8)
    atoms = np.zeros((512, 8), np.int16)
    uv = np.zeros((512, 2))
    face = np.zeros(512, np.int16)
    tags = []
    for ipt, gid, tag, u, v in G:
        c = m1atoms.codes_at(u, v)                                    # exact, at the simulated game
        codes[ipt] = c
        atoms[ipt] = np.bincount(c, minlength=8)
        uv[ipt] = (float(u), float(v))
        tags.append(tag)
        f = face_of(u, v, lines, uvc)
        face[ipt] = f
        assert list(A["codes"][f]) == list(c), \
            "%s: the exact codes differ from face %d of m1atoms.py's arrangement" % (tag, f)
    assert (atoms.sum(1) == NSTRAT).all(), "the atom counts do not sum to 16 at every game"
    assert not atoms[:, 5].any(), "101 is NOT empty at every game -- the seven-panel premise fails"
    pat = {tuple((atoms[i] > 0).astype(int)) for i in range(512)}
    allpat = {tuple((A["atoms"][f] > 0).astype(int)) for f in range(len(A["atoms"]))}
    assert pat <= allpat, "a pattern occurs at a game that m1atoms.py did not find on the plane"
    np.savez_compressed(DUMP, uv=uv, codes=codes, atoms=atoms, face=face,
                        tags=np.array(tags), lines=np.array(lines))
    return dict(uv=uv, codes=codes, atoms=atoms, face=face, tags=np.array(tags),
                lines=np.array(lines), npat=len(pat))


def load(rebuild=False):
    if rebuild or not os.path.exists(DUMP):
        return build()
    d = np.load(DUMP, allow_pickle=False)
    return {k: d[k] for k in d.files}


def census(d):
    atoms, codes = d["atoms"], d["codes"]
    u, v = d["uv"][:, 0], d["uv"][:, 1]
    quad = np.where((u < 0) & (v > 0), "PD", np.where((u > 0) & (v > 0), "SD",
                    np.where((u < 0) & (v < 0), "SH", "HG")))
    print("the seven atoms of the 16 binary memory-one strategies at the 512 sunflower games")
    print("%-5s %-34s %8s %7s %7s   %s" % ("atom", "", "per game", "empty", "games", "strategies"))
    for c in ORDER:
        i = int(c, 2)
        n = atoms[:, i]
        who = sorted({k for g in range(512) for k in range(NSTRAT) if codes[g][k] == i},
                     key=lambda k: (-sum(codes[g][k] == i for g in range(512)), k))
        nm = " ".join(genome(k) + ("(%s)" % NICK[k] if k in NICK else "") for k in who)
        print("%-5s %-34s %4d-%-3d %7d %7d   %d: %s"
              % (c, NAME[c], n.min(), n.max(), int((n == 0).sum()), int((n > 0).sum()), len(who), nm))
    print("101 (efficient and competitive, not nash): empty at all 512 games")
    pat = {}
    for g in range(512):
        pat.setdefault(tuple((atoms[g] > 0).astype(int)), []).append(g)
    print("%d patterns of non-empty atoms occur among the 512 games (twelve occur on the plane):" % len(pat))
    for k in sorted(pat, key=lambda k: -len(pat[k])):
        pres = " ".join(ORDER[j] for j, c in enumerate(ORDER) if k[int(c, 2)])
        q = " ".join("%s %d" % (qq, sum(quad[g] == qq for g in pat[k]))
                     for qq in ("PD", "SD", "SH", "HG") if any(quad[g] == qq for g in pat[k]))
        print("   %4d games   %-32s  %s" % (len(pat[k]), pres, q))
    print("faces of m1atoms.py's 45-face arrangement visited: %d" % len(set(d["face"].tolist())))


if __name__ == "__main__":
    census(load("rebuild" in sys.argv))
    print("cache: %s" % DUMP)
