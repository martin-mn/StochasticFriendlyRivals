#!/usr/bin/env python3
"""
disc_check.py -- DISCOUNTED friendly rivals over the memory-one continuum.

Route 4, numerical.  Everything here is MEASURED; the closed forms printed alongside the
measurements are predictions that the script checks, not assumptions it uses.

Conventions (those of ../m1lib.py, not redefined):
  * base game (R,S,T,P) = (1, u, 1+v, 0);  T-S = 1+v-u;  2R-(T+S) = 1-u-v
  * p = (p1,p2,p3,p4) = P(C | CC, CD, DC, DD), focal action first; states CC=0 CD=1 DC=2 DD=3
  * convention A: implementation error acts as x -> eps + (1-2 eps) x, applied to p1..p4 AND
    to the round-0 unconditional cooperation probability p0.
  * a strategy is the 5-tuple (p0; p1,p2,p3,p4).

Discounted payoff (this is the only thing new):
      v  = (1-delta) v0 (I - delta M)^{-1}          (row vector, sum(v) = 1)
      v0 = (p0 q0, p0(1-q0), (1-p0) q0, (1-p0)(1-q0))
      pi_delta(p,q) = v.(R,S,T,P) = v.(1,u,1+v,0)
      pi_delta(q,p) = v.(R,T,S,P) = v.(1,1+v,u,0)
      pi_delta(p,q) - pi_delta(q,p) = (T-S)(v_DC - v_CD)

Definitions used, unchanged from the undiscounted project:
  efficient  : pi_delta(p,p) = Emax = max{R,(S+T)/2} = max{1,(1+u+v)/2}
  rival      : pi_delta(p,q) >= pi_delta(q,p) for EVERY co-player (q0;q)
  friendly rival : both.

Run:  python3 disc_check.py            (about 3-5 minutes on 12 cores)
"""

import itertools
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import numpy as np

T0 = time.time()
def el():
    return "[%6.1fs]" % (time.time() - T0)

SWAP = (0, 2, 1, 3)                      # state seen from the co-player's side
DET = np.array(list(itertools.product([0, 1], repeat=4)), float)   # LSB-first genome order
NAMES = {"0000": "ALLD", "1111": "ALLC", "1010": "TFT", "1001": "WSLS", "1000": "Grim",
         "0101": "ATFT", "0110": "AWSLS"}

def genome(p):
    return "".join("1" if x > 0.5 else "0" for x in p)

def gname(p):
    g = genome(p)
    return "%s %s" % (g, NAMES.get(g, ""))

# --------------------------------------------------------------------- float route (batched)

def clamp(X, eps):
    return eps + (1.0 - 2.0 * eps) * np.asarray(X, float)

def trans(P, Q, eps):
    """(n,4,4) row-stochastic pair-chain matrices, focal P against co-player Q."""
    P = clamp(P, eps); Q = clamp(Q, eps)
    n = P.shape[0]
    M = np.empty((n, 4, 4))
    for s in range(4):
        a = P[:, s]; b = Q[:, SWAP[s]]
        M[:, s, 0] = a * b
        M[:, s, 1] = a * (1.0 - b)
        M[:, s, 2] = (1.0 - a) * b
        M[:, s, 3] = (1.0 - a) * (1.0 - b)
    return M

def disc_v(P, p0, Q, q0, delta, eps):
    """v = (1-delta) v0 (I - delta M)^{-1}, batched over rows of P (and of Q if 2-D)."""
    P = np.atleast_2d(np.asarray(P, float))
    Q = np.atleast_2d(np.asarray(Q, float))
    n = max(P.shape[0], Q.shape[0])
    if P.shape[0] == 1 and n > 1: P = np.repeat(P, n, 0)
    if Q.shape[0] == 1 and n > 1: Q = np.repeat(Q, n, 0)
    a0 = clamp(np.broadcast_to(np.asarray(p0, float), (n,)), eps)
    b0 = clamp(np.broadcast_to(np.asarray(q0, float), (n,)), eps)
    v0 = np.stack([a0 * b0, a0 * (1 - b0), (1 - a0) * b0, (1 - a0) * (1 - b0)], axis=1)
    M = trans(P, Q, eps)
    A = np.eye(4)[None, :, :] - delta * M           # v (I - delta M) = (1-delta) v0
    x = np.linalg.solve(np.swapaxes(A, 1, 2), ((1.0 - delta) * v0)[:, :, None])[:, :, 0]
    return x

def pay(v, u, vv):
    return v[:, 0] + v[:, 1] * u + v[:, 2] * (1.0 + vv)

def pay_co(v, u, vv):
    return v[:, 0] + v[:, 1] * (1.0 + vv) + v[:, 2] * u

def Emax(u, vv):
    return max(1.0, (1.0 + u + vv) / 2.0)

def eff_deficit(P, p0, u, vv, delta, eps):
    v = disc_v(P, p0, P, p0, delta, eps)
    return pay(v, u, vv) - Emax(u, vv)

def margin_vs_set(P, p0, Qs, q0s, delta, eps, u, vv):
    """min over the given co-player list of pi_delta(p,q) - pi_delta(q,p), per candidate."""
    P = np.atleast_2d(P)
    n = P.shape[0]
    best = np.full(n, np.inf)
    for q, q0 in zip(Qs, q0s):
        v = disc_v(P, p0, np.asarray(q, float)[None, :], q0, delta, eps)
        np.minimum(best, pay(v, u, vv) - pay_co(v, u, vv), out=best)
    return best

# the 32 deterministic co-players: 16 vertices x q0 in {0,1}.  For a DISCOUNTED MDP the
# co-player's optimal reply is a deterministic stationary policy plus a deterministic q0,
# so these 32 are the exact minimisers over the whole box -- which the random stage checks.
DET32_Q = np.repeat(DET, 2, axis=0)
DET32_Q0 = np.tile(np.array([0.0, 1.0]), 16)

# ------------------------------------------------------------------------ exact route

def _exact_solve(A, b):
    n = len(b)
    A = [row[:] for row in A]; b = b[:]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]; b[c], b[piv] = b[piv], b[c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
                b[r] -= f * b[c]
    x = [Fraction(0)] * n
    for r in reversed(range(n)):
        x[r] = (b[r] - sum(A[r][k] * x[k] for k in range(r + 1, n))) / A[r][r]
    return x

def disc_v_exact(p, p0, q, q0, delta, eps):
    d = Fraction(delta); e = Fraction(eps)
    cl = lambda x: e + (1 - 2 * e) * Fraction(x)
    pc = [cl(x) for x in p]; qc = [cl(x) for x in q]
    a0 = cl(p0); b0 = cl(q0)
    v0 = [a0 * b0, a0 * (1 - b0), (1 - a0) * b0, (1 - a0) * (1 - b0)]
    M = [[None] * 4 for _ in range(4)]
    for s in range(4):
        a, b = pc[s], qc[SWAP[s]]
        M[s] = [a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)]
    # v (I - d M) = (1-d) v0   ->   (I - d M)^T v = (1-d) v0
    A = [[(Fraction(1) if i == j else Fraction(0)) - d * M[j][i] for j in range(4)]
         for i in range(4)]
    return _exact_solve(A, [(1 - d) * x for x in v0])

def pay_exact(v, u, vv):
    u = Fraction(u); vv = Fraction(vv)
    return v[0] + v[1] * u + v[2] * (1 + vv)

def pay_co_exact(v, u, vv):
    u = Fraction(u); vv = Fraction(vv)
    return v[0] + v[1] * (1 + vv) + v[2] * u

# ------------------------------------------------------------------------ the games

GAMES = [
    # label            u     v     wedge / line
    ("W  (0.2,0.2)",   0.2,  0.2),     # u+v<1, u-v<1   PD, T>S, CC optimal
    ("S  (0.4,-1.0)",  0.4, -1.0),     # u+v<1, u-v>1   T<S, CC optimal
    ("E  (2.0,0.5)",   2.0,  0.5),     # u+v>1, u-v>1   T<S, alternation optimal
    ("N  (0.5,1.5)",   0.5,  1.5),     # u+v>1, u-v<1   T>S, alternation optimal
    ("L+ (0.3,0.7)",   0.3,  0.7),     # u+v=1  (W|N boundary)
    ("L+ (1.5,-0.5)",  1.5, -0.5),     # u+v=1  (S|E boundary)
    ("L- (0.5,-0.5)",  0.5, -0.5),     # u-v=1, T=S  (W|S boundary)
    ("L- (1.5,0.5)",   1.5,  0.5),     # u-v=1, T=S  (N|E boundary)
]

DELTAS = [0.5, 0.9, 0.99, 0.999, 0.99999]
EPSES = [0.0, 1e-7, 1e-5, 1e-3]        # eps = 0 is legitimate here: I - delta M is
                                       # invertible for every delta < 1, ergodic or not.

# --------------------------------------------------------------- the candidate continuum

def build_candidates(seed=20260912):
    rng = np.random.default_rng(seed)
    P = []; p0 = []; tag = []
    # (a) all 16 binary vertices x p0 in {0, 1/4, 1/2, 3/4, 1}
    for b in DET:
        for z in (0.0, 0.25, 0.5, 0.75, 1.0):
            P.append(b); p0.append(z); tag.append("vertex")
    # (b) regular 5^4 grid x the same five p0
    g = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    g7 = np.linspace(0.0, 1.0, 7)
    for c in itertools.product(g7, repeat=4):
        for z in (0.0, 1.0):
            P.append(c); p0.append(z); tag.append("grid7")
    for c in itertools.product(g, repeat=4):
        for z in g:
            P.append(c); p0.append(z); tag.append("grid")
    # (c) stratified random: all 81 strata (each p_i at 0, interior, or 1), 250 each,
    #     p0 cycling through {0, 1, uniform}
    K = 600
    for st in itertools.product([0, 1, 2], repeat=4):
        for k in range(K):
            c = [0.0 if s == 0 else (1.0 if s == 2 else float(rng.random())) for s in st]
            z = (0.0, 1.0, float(rng.random()))[k % 3]
            P.append(c); p0.append(z); tag.append("strat")
    # (d) extra resolution where the undiscounted theory puts the friendly rivals:
    #     the facet p1 = 1 with p0 = 1, 6000 uniform points
    for _ in range(20000):
        P.append([1.0, float(rng.random()), float(rng.random()), float(rng.random())])
        p0.append(1.0); tag.append("facet")
    return np.array(P, float), np.array(p0, float), np.array(tag)

CAND_P, CAND_P0, CAND_TAG = build_candidates()
NCAND = CAND_P.shape[0]

# ------------------------------------------------------------------------ the sweep cell

NRAND_Q = 3000          # random co-players applied to the stage-2 survivors
MAX_STAGE3 = 400        # cap on how many survivors get the random-co-player treatment

def cell(args):
    gi, di, ei = args
    label, u, vv = GAMES[gi]
    delta = DELTAS[di]; eps = EPSES[ei]
    tol = 10.0 * eps + 1e-12                      # the STATED tolerance
    tight = 1e-9

    dfc = eff_deficit(CAND_P, CAND_P0, u, vv, delta, eps)
    eff = dfc >= -tol
    eff_t = dfc >= -tight

    out = {"game": label, "u": u, "v": vv, "delta": delta, "eps": eps,
           "n_cand": NCAND, "n_eff": int(eff.sum()), "n_eff_tight": int(eff_t.sum()),
           "best_deficit": float(dfc.max())}

    def stage23(mask, tol_r):
        idx = np.nonzero(mask)[0]
        if idx.size == 0:
            return idx, np.array([]), 0
        m = margin_vs_set(CAND_P[idx], CAND_P0[idx], DET32_Q, DET32_Q0, delta, eps, u, vv)
        keep = m >= -tol_r
        idx2 = idx[keep]; m2 = m[keep]
        if idx2.size == 0:
            return idx2, m2, 0
        rng = np.random.default_rng(1000 * gi + 10 * di + ei)
        sub = idx2 if idx2.size <= MAX_STAGE3 else idx2[rng.choice(idx2.size, MAX_STAGE3, replace=False)]
        Qr = rng.random((NRAND_Q, 4)); q0r = rng.random(NRAND_Q)
        mr = margin_vs_set(CAND_P[sub], CAND_P0[sub], Qr, q0r, delta, eps, u, vv)
        broken = int((mr < -tol_r).sum())
        return idx2, m2, broken

    idx2, m2, broken = stage23(eff, tol)
    out["n_fr"] = int(idx2.size)
    out["n_fr_broken_by_random"] = broken
    out["best_margin_eff"] = float(m2.max()) if m2.size else float("nan")

    idx2t, m2t, broken_t = stage23(eff_t, tight)
    out["n_fr_tight"] = int(idx2t.size)
    out["n_fr_tight_broken"] = broken_t

    # best rivalry margin attainable among TIGHTLY-efficient candidates (item 4)
    if eff_t.any():
        mm = margin_vs_set(CAND_P[eff_t], CAND_P0[eff_t], DET32_Q, DET32_Q0, delta, eps, u, vv)
        out["best_riv_margin_among_eff"] = float(mm.max())
    else:
        out["best_riv_margin_among_eff"] = float("nan")

    # the closed form this script predicts for the discounted friendly rivals (checked, not used)
    pr = predict(CAND_P, CAND_P0, u, vv, delta)
    out["n_pred"] = int(pr.sum())
    agree = int((pr == (eff_t & _mask_from(idx2t, NCAND))).sum())
    out["pred_disagree"] = NCAND - agree
    return out

def _mask_from(idx, n):
    m = np.zeros(n, bool); m[idx] = True; return m

def predict(P, p0, u, vv, delta, tol=1e-12):
    """PREDICTED discounted friendly rivals (eps = 0):

       exact efficiency  <=>  u+v <= 1, p0 = 1, p1 = 1          (CC forever from round 0)
       rivalry, T < S    <=>  p2 + p3 >= 2 - 1/delta.  NOTE: the undiscounted exclusion
                              "not (p2=1 and p3=0)" is DROPPED -- for delta < 1 the identity
                              forces v_DC = v_DD = 0 there, so that edge is a rival after all.
       rivalry, T = S    <=>  vacuous
       rivalry, T > S    <=>  impossible once p0 = 1 (ALLD gives -(1-delta)(T-S) < 0)
    """
    p1, p2, p3, p4 = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    base = (p0 >= 1 - tol) & (p1 >= 1 - tol) & (u + vv <= 1 + tol)
    TS = 1.0 + vv - u                      # T - S
    if TS > tol:                           # T > S
        return np.zeros(P.shape[0], bool)
    if abs(TS) <= tol:                     # T = S, rivalry vacuous
        return base
    thr = 2.0 - 1.0 / delta
    return base & (p2 + p3 >= thr - tol)

# ------------------------------------------------------------------------------- main

def main():
    print("=" * 100)
    print("DISCOUNTED friendly rivals over the memory-one continuum -- measured")
    print("candidates: %d  (16 vertices x5 p0; 5^4 grid x5 p0; 7^4 grid x2 p0; 81 strata x600; 20000 on p1=1)"
          % NCAND)
    print("co-players: 32 deterministic (16 vertices x q0 in {0,1}) then %d random per survivor"
          % NRAND_Q)
    print("=" * 100)

    # ---------------------------------------------------------------- 0. exact receipts
    print("\n%s === 0. EXACT (Fraction) SPOT CHECKS =====================================" % el())
    ok = True
    for delta in (Fraction(1, 2), Fraction(9, 10), Fraction(99, 100)):
        for eps in (Fraction(0), Fraction(1, 1000)):
            for p, p0, q, q0 in [((1, 0, 1, 0), 1, (0, 0, 0, 0), 0),
                                 ((1, 1, 1, 1), 1, (1, 0, 0, 1), 0),
                                 ((Fraction(3, 7), Fraction(1, 3), 1, Fraction(2, 5)),
                                  Fraction(4, 9), (0, 1, 1, 0), 1)]:
                v = disc_v_exact(p, p0, q, q0, delta, eps)
                if sum(v) != 1:
                    print("   FAIL normalisation", p, q, delta, eps, sum(v)); ok = False
                # discounted Akin identity:
                #   v_CC(1-d p1)+v_CD(1-d p2) - d p3 v_DC - d p4 v_DD = (1-d) p0   (clamped)
                cl = lambda x: eps + (1 - 2 * eps) * Fraction(x)
                pc = [cl(x) for x in p]
                res = (v[0] * (1 - delta * pc[0]) + v[1] * (1 - delta * pc[1])
                       - delta * pc[2] * v[2] - delta * pc[3] * v[3] - (1 - delta) * cl(p0))
                if res != 0:
                    print("   FAIL Akin_delta", p, q, delta, eps, res); ok = False
    print("   normalisation sum(v) == 1 exactly, and the discounted Akin identity")
    print("     v_CC(1-d p1) + v_CD(1-d p2) - d p3 v_DC - d p4 v_DD = (1-d) p0")
    print("   hold exactly on all spot checks:", "OK" if ok else "FAILED")

    # float normalisation over the whole candidate set
    worst = 0.0
    for delta in DELTAS:
        for eps in EPSES:
            v = disc_v(CAND_P, CAND_P0, CAND_P, CAND_P0, delta, eps)
            worst = max(worst, float(np.abs(v.sum(1) - 1).max()))
    print("   float: max |sum(v) - 1| over all %d candidates x %d delta x %d eps = %.2e"
          % (NCAND, len(DELTAS), len(EPSES), worst))

    # ---------------------------------------------------------------- 1. the big table
    print("\n%s === 1. SWEEP: efficient AND unbeaten, by (game, delta, eps) =============" % el())
    jobs = [(gi, di, ei) for gi in range(len(GAMES)) for di in range(len(DELTAS))
            for ei in range(len(EPSES))]
    with Pool(12) as pool:
        res = pool.map(cell, jobs)
    R = {(r["game"], r["delta"], r["eps"]): r for r in res}

    hdr = "  %-16s %-9s" % ("game (u,v)", "delta") + "".join("   eps=%-8g" % e for e in EPSES)
    print("\n  TABLE 1a.  count of candidates that are efficient AND unbeaten")
    print("             tolerance = 10*eps + 1e-12 on BOTH predicates;  out of %d" % NCAND)
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for label, u, vv in GAMES:
        for delta in DELTAS:
            row = "  %-16s %-9g" % (label if delta == DELTAS[0] else "", delta)
            for e in EPSES:
                r = R[(label, delta, e)]
                row += "   %11d" % r["n_fr"]
            print(row)
        print()

    print("  TABLE 1b.  the same count with a TIGHT absolute tolerance 1e-9 on both predicates")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for label, u, vv in GAMES:
        for delta in DELTAS:
            row = "  %-16s %-9g" % (label if delta == DELTAS[0] else "", delta)
            for e in EPSES:
                row += "   %11d" % R[(label, delta, e)]["n_fr_tight"]
            print(row)
        print()

    print("  TABLE 1c.  how many of the counted ones were then BROKEN by a random co-player")
    print("             (0 everywhere = the 32 deterministic co-players really are the minimisers)")
    tot = sum(r["n_fr_broken_by_random"] + r["n_fr_tight_broken"] for r in res)
    print("             total broken over all %d cells: %d" % (len(res), tot))

    print("\n  TABLE 1d.  candidates that are efficient at tolerance 1e-9 (before rivalry)")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for label, u, vv in GAMES:
        for delta in DELTAS:
            row = "  %-16s %-9g" % (label if delta == DELTAS[0] else "", delta)
            for e in EPSES:
                row += "   %11d" % R[(label, delta, e)]["n_eff_tight"]
            print(row)
        print()

    print("  TABLE 1e.  best (largest) efficiency deficit achievable, tolerance-free, eps=0")
    print("             = max_p [ pi_delta(p,p) - Emax ];  0 means exact efficiency is attainable")
    print("  %-16s" % "game" + "".join("  delta=%-9g" % d for d in DELTAS))
    for label, u, vv in GAMES:
        row = "  %-16s" % label
        for d in DELTAS:
            row += "  %15.3e" % R[(label, d, 0.0)]["best_deficit"]
        print(row)

    print("\n  TABLE 1f.  closed-form prediction vs measurement, eps = 0, tolerance 1e-9")
    print("             prediction: efficient <=> u+v<=1 & p0=1 & p1=1;")
    print("                         unbeaten  <=> T<S and p2+p3 >= 2-1/delta, or T=S (vacuous)")
    print("  %-16s %-9s %10s %10s %10s" % ("game", "delta", "predicted", "measured", "disagree"))
    baddis = 0
    for label, u, vv in GAMES:
        for d in DELTAS:
            r = R[(label, d, 0.0)]
            baddis += r["pred_disagree"]
            print("  %-16s %-9g %10d %10d %10d"
                  % (label, d, r["n_pred"], r["n_fr_tight"], r["pred_disagree"]))
    print("  TOTAL DISAGREEMENTS (eps=0):", baddis)

    # ---------------------------------------------------- 2. the undiscounted wedge-S set
    print("\n%s === 2. DO THE UNDISCOUNTED WEDGE-S FRIENDLY RIVALS SURVIVE delta < 1? ====" % el())
    label, u, vv = GAMES[1]
    print("  game %s:  R=1 S=%g T=%g P=0,  T-S=%g (<0),  Emax=%g"
          % (label, u, 1 + vv, 1 + vv - u, Emax(u, vv)))
    probes = [("ALLC   1111", [1, 1, 1, 1]), ("       1110", [1, 1, 1, 0]),
              ("       1101", [1, 1, 0, 1]), ("       1011", [1, 0, 1, 1]),
              ("       1010 TFT", [1, 0, 1, 0]), ("       1000 Grim", [1, 0, 0, 0]),
              ("       1001 WSLS", [1, 0, 0, 1]),
              ("interior .6/.6/.4", [1, .6, .6, .4]), ("interior .5/.5/0", [1, .5, .5, 0]),
              ("interior .8/.3/.5", [1, .8, .3, .5]), ("interior .3/.3/.5", [1, .3, .3, .5]),
              ("interior .2/.2/.2", [1, .2, .2, .2]), ("p1=.99 .6/.6/.4", [.99, .6, .6, .4])]
    for p0v in (1.0, 0.5):
        print("\n  --- p0 = %g   (columns: efficiency deficit / rivalry margin;"
              "  FR = both >= -1e-9) ---" % p0v)
        print("  %-18s %8s" % ("strategy", "p2+p3") +
              "".join("      d=%-11g" % d for d in DELTAS))
        for nm, p in probes:
            row = "  %-18s %8.2f" % (nm, p[1] + p[2])
            for d in DELTAS:
                P = np.array([p], float)
                df = float(eff_deficit(P, p0v, u, vv, d, 0.0)[0])
                mg = float(margin_vs_set(P, p0v, DET32_Q, DET32_Q0, d, 0.0, u, vv)[0])
                fr = "FR " if (df >= -1e-9 and mg >= -1e-9) else "   "
                row += "  %s%+.1e/%+.1e" % (fr, df, mg)
            print(row)

    print("\n  --- measured rivalry threshold on the segment p=(1,s,1-s+g,0.5), p0=1, eps=0 ---")
    print("      (bisection in x = p2+p3 on the family p2 = p3 = x/2, p4 = 0.5)")
    print("  %-9s %14s %14s %14s" % ("delta", "measured x*", "2 - 1/delta", "difference"))
    for d in DELTAS + [0.6, 0.75, 0.8]:
        lo, hi = 0.0, 2.0
        for _ in range(60):
            mid = (lo + hi) / 2
            P = np.array([[1.0, mid / 2, mid / 2, 0.5]])
            if float(margin_vs_set(P, 1.0, DET32_Q, DET32_Q0, d, 0.0, u, vv)[0]) >= 0:
                hi = mid
            else:
                lo = mid
        pred = max(0.0, 2.0 - 1.0 / d)
        print("  %-9g %14.9f %14.9f %14.2e" % (d, hi, pred, hi - pred))

    # ----------------------------------------------- 3. the classic PD candidates
    print("\n%s === 3. THE PD SIDE: exact losing margin against ALLD =====================" % el())
    label, u, vv = GAMES[0]
    TS = 1 + vv - u
    print("  game %s:  R=1 S=%g T=%g P=0,  T-S=%g,  Emax=%g  (a Prisoner's Dilemma)"
          % (label, u, 1 + vv, TS, Emax(u, vv)))
    print("\n  margin against ALLD (q0=0), eps=0, p0=1;  predicted -(1-delta)(T-S) for any")
    print("  strategy with p1=1 that never cooperates again after being suckered (p2=p4=0)")
    print("  %-16s" % "strategy" + "".join("      delta=%-9g" % d for d in DELTAS))
    allD = np.array([[0., 0., 0., 0.]])
    pdprobes = [("1010 TFT", [1, 0, 1, 0]), ("1000 Grim", [1, 0, 0, 0]),
                ("1001 WSLS", [1, 0, 0, 1]), ("1111 ALLC", [1, 1, 1, 1]),
                ("1100", [1, 1, 0, 0]), ("1011", [1, 0, 1, 1]),
                ("(1,0,.5,0)", [1, 0, .5, 0]), ("(1,.2,.3,0)", [1, .2, .3, 0])]
    for nm, p in pdprobes:
        row = "  %-16s" % nm
        for d in DELTAS:
            v = disc_v(np.array([p], float), 1.0, allD, 0.0, d, 0.0)
            row += "  %17.10e" % float(pay(v, u, vv)[0] - pay_co(v, u, vv)[0])
        print(row)
    row = "  %-16s" % "-(1-d)(T-S)"
    for d in DELTAS:
        row += "  %17.10e" % (-(1 - d) * TS)
    print(row)

    print("\n  exact (Fraction) identity check, TFT/Grim with p0=1 vs ALLD, eps=0:")
    for nm, p in (("TFT", (1, 0, 1, 0)), ("Grim", (1, 0, 0, 0)), ("(1,0,1/2,0)", (1, 0, Fraction(1, 2), 0))):
        for d in (Fraction(1, 2), Fraction(9, 10), Fraction(999, 1000), Fraction(99999, 100000)):
            v = disc_v_exact(p, 1, (0, 0, 0, 0), 0, d, Fraction(0))
            got = pay_exact(v, Fraction(1, 5), Fraction(1, 5)) - pay_co_exact(v, Fraction(1, 5), Fraction(1, 5))
            want = -(1 - d) * (1 + Fraction(1, 5) - Fraction(1, 5))
            print("    %-12s delta=%-10s margin = %-22s  -(1-d)(T-S) = %-22s  %s"
                  % (nm, d, got, want, "EQUAL" if got == want else "DIFFER"))
    print("  (WSLS differs: p4=1 makes it re-cooperate from DD, so it is suckered every other")
    print("   round -- margin -(T-S)/(1+delta), which does NOT vanish as delta -> 1.)")
    print("  WSLS exact:", end=" ")
    for d in (Fraction(9, 10), Fraction(999, 1000)):
        v = disc_v_exact((1, 0, 0, 1), 1, (0, 0, 0, 0), 0, d, Fraction(0))
        got = pay_exact(v, Fraction(1, 5), Fraction(1, 5)) - pay_co_exact(v, Fraction(1, 5), Fraction(1, 5))
        print("d=%s -> %s  (-1/(1+d) = %s)" % (d, got, -Fraction(1) / (1 + d)), end="   ")
    print()

    print("\n  the WHOLE facet p1 = 1 (p0 = 1), eps = 0: best margin against the 32 co-players")
    print("  over 20000 uniform points of the facet, vs the predicted ceiling -(1-delta)(T-S)")
    fac = CAND_TAG == "facet"
    print("  %-9s %18s %18s %14s" % ("delta", "best margin", "-(1-delta)(T-S)", "argmax p"))
    for d in DELTAS:
        m = margin_vs_set(CAND_P[fac], 1.0, DET32_Q, DET32_Q0, d, 0.0, u, vv)
        j = int(np.argmax(m))
        print("  %-9g %18.10e %18.10e   %s" % (d, float(m.max()), -(1 - d) * TS,
              np.round(CAND_P[fac][j], 3)))
    print("  and the two exact PD-optimal points (1,0,p3,0), p3 free:")
    for d in DELTAS:
        P = np.array([[1., 0., t, 0.] for t in (0., .25, .5, .75, 1.)])
        m = margin_vs_set(P, 1.0, DET32_Q, DET32_Q0, d, 0.0, u, vv)
        print("    delta=%-9g margins %s   ceiling %.6e" % (d, np.array2string(m, precision=6), -(1 - d) * TS))

    # -------------------------------------------------- 4. scaling of the best margin
    print("\n%s === 4. SCALING OF THE BEST ACHIEVABLE RIVALRY MARGIN AS delta -> 1 =======" % el())
    print("  best_riv = max over EXACTLY-efficient candidates of min over the 32 co-players")
    print("  of pi_delta(p,q) - pi_delta(q,p).   eps = 0.   nan = no efficient candidate exists.")
    ds = [0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999]
    print("  %-16s" % "game" + "".join("  delta=%-10g" % d for d in ds))
    for label, u, vv in GAMES:
        row = "  %-16s" % label
        for d in ds:
            df = eff_deficit(CAND_P, CAND_P0, u, vv, d, 0.0)
            m = df >= -1e-9
            if not m.any():
                row += "  %16s" % "nan"
            else:
                mm = margin_vs_set(CAND_P[m], CAND_P0[m], DET32_Q, DET32_Q0, d, 0.0, u, vv)
                row += "  %16.7e" % float(mm.max())
        print(row)
    print("\n  ratio best_riv / (1-delta) -- a finite non-zero limit means the failure (or the")
    print("  success) is exactly of order (1-delta):")
    print("  %-16s" % "game" + "".join("  delta=%-10g" % d for d in ds))
    for label, u, vv in GAMES:
        row = "  %-16s" % label
        for d in ds:
            df = eff_deficit(CAND_P, CAND_P0, u, vv, d, 0.0)
            m = df >= -1e-9
            if not m.any():
                row += "  %16s" % "nan"
            else:
                mm = margin_vs_set(CAND_P[m], CAND_P0[m], DET32_Q, DET32_Q0, d, 0.0, u, vv)
                row += "  %16.7f" % (float(mm.max()) / (1 - d))
        print(row)
    print("\n  T-S per game:", ", ".join("%s: %.3g" % (g[0], 1 + g[2] - g[1]) for g in GAMES))


    # ------------------------------------------- 5. VOLUME of the discounted friendly-rival set
    print("\n%s === 5. VOLUME of the friendly-rival set inside the facet p1 = 1 =========" % el())
    print("  wedge S, p0 = 1, eps = 0.  Monte Carlo over (p2,p3,p4) uniform in the unit cube,")
    print("  N = 400000.  Predicted volume = 1 - max(0, 2-1/delta)^2 / 2, which is 1 at")
    print("  delta <= 1/2 and decreases monotonically to the undiscounted 1/2 at delta = 1.")
    label, u, vv = GAMES[1]
    rng = np.random.default_rng(7)
    N = 400000
    Pm = np.column_stack([np.ones(N), rng.random((N, 3))])
    print("  %-9s %14s %14s %12s %14s" % ("delta", "measured vol", "1-(2-1/d)^2/2", "diff", "undisc. 1/2"))
    for d in [0.25, 0.5, 0.6, 0.75, 0.9, 0.99, 0.999, 0.99999, 0.9999999]:
        m = margin_vs_set(Pm, 1.0, DET32_Q, DET32_Q0, d, 0.0, u, vv)
        vol = float((m >= -1e-12).mean())
        thr = max(0.0, 2.0 - 1.0 / d)
        pred = 1.0 - thr * thr / 2.0
        print("  %-9g %14.6f %14.6f %12.2e %14.6f" % (d, vol, pred, vol - pred, 0.5))

    # ---------------------------------------------- 6. the eps <-> (1-delta) order of limits
    print("\n%s === 6. ORDER OF LIMITS: eps vs (1-delta) =================================" % el())
    print("  wedge S, p0 = 1.  Efficiency deficit and rivalry margin of three wedge-S friendly")
    print("  rivals at finite eps.  Both are O(eps) at fixed delta: the eps -> 0 friendly")
    print("  rivals are eps-friendly-rivals at finite eps, exactly as undiscounted.")
    label, u, vv = GAMES[1]
    for nm, p in (("ALLC 1111", [1., 1., 1., 1.]), ("1110", [1., 1., 1., 0.]),
                  ("(1,.6,.6,.4)", [1., .6, .6, .4]), ("TFT 1010", [1., 0., 1., 0.])):
        print("  %s" % nm)
        print("    %-9s" % "delta" + "".join("     eps=%-9g" % e for e in EPSES[1:]))
        for d in DELTAS:
            row = "    %-9g" % d
            for e in EPSES[1:]:
                P = np.array([p])
                df = float(eff_deficit(P, 1.0, u, vv, d, e)[0])
                mg = float(margin_vs_set(P, 1.0, DET32_Q, DET32_Q0, d, e, u, vv)[0])
                row += "  %+.1e/%+.1e" % (df / e, mg / e)
            print(row + "      (values are deficit/eps and margin/eps)")

    print("\n  wedge E: exact efficiency is unattainable for EVERY p0 and every delta < 1,")
    print("  because round 0 has both players drawing independently from the same p0 and so")
    print("  cannot be perfectly anti-correlated.  deficit / (1-delta), eps = 0:")
    for gi in (2, 3, 7):
        label, u, vv = GAMES[gi]
        row = "  %-16s" % label
        for d in [0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999]:
            row += "  %12.6f" % (float(eff_deficit(CAND_P, CAND_P0, u, vv, d, 0.0).max()) / (1 - d))
        print(row)
    print("  %-16s" % "" + "".join("  delta=%-6g" % d for d in [0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999]))
    print("  max over p0 of the round-0 payoff p0^2 + p0(1-p0)(1+u+v), against Emax:")
    for gi in (2, 3, 7):
        label, u, vv = GAMES[gi]
        z = np.linspace(0, 1, 200001)
        f = z * z + z * (1 - z) * (1 + u + vv)
        print("    %-16s max round-0 = %.6f   Emax = %.6f   gap = %.6f"
              % (label, f.max(), Emax(u, vv), f.max() - Emax(u, vv)))

    # p0 clamping sensitivity
    label, u, vv = GAMES[1]                       # back to wedge S
    print("\n  DO THE LIMITS COMMUTE?  eps=0 efficiency for delta<1 is just {p0=1, p1=1}, which")
    print("  is STRICTLY BIGGER than the undiscounted Eff_CC: the alternating corners")
    print("  (p2,p3) in {(0,1),(1,0)} are excluded undiscounted (they trap the eps-excursion for")
    print("  Theta(1/eps) rounds) but are harmless at delta < 1 with eps = 0.  Wedge S deficit:")
    print("  %-16s %14s %14s %14s %14s" % ("strategy", "eps=0,d=.99999", "eps=1e-7,d->1",
                                           "eps=1e-5,d->1", "undisc. eff?"))
    for nm, p, ue in (("TFT 1010", [1., 0., 1., 0.], "NO (alt corner)"),
                      ("0110-like (1,0,1,1)", [1., 0., 1., 1.], "NO (alt corner)"),
                      ("(1,1,0,.5)", [1., 1., 0., .5], "NO (alt corner)"),
                      ("ALLC 1111", [1., 1., 1., 1.], "yes"),
                      ("1110", [1., 1., 1., 0.], "yes"),
                      ("(1,.6,.6,.4)", [1., .6, .6, .4], "yes")):
        P = np.array([p])
        a = float(eff_deficit(P, 1.0, u, vv, 0.99999, 0.0)[0])
        b = float(eff_deficit(P, 1.0, u, vv, 0.9999999, 1e-7)[0])
        c = float(eff_deficit(P, 1.0, u, vv, 0.9999999, 1e-5)[0])
        print("  %-16s %14.3e %14.3e %14.3e   %s" % (nm, a, b, c, ue))

    print("\n  round-0 clamping: the script applies eps to p0 as well.  Effect at eps=1e-3,")
    print("  wedge S, ALLC, p0=1:  deficit with / without clamping of p0:")
    for d in DELTAS:
        P = np.array([[1., 1., 1., 1.]])
        a = float(eff_deficit(P, 1.0, u, vv, d, 1e-3)[0])
        b = float(eff_deficit(P, 1.0 / (1 - 2e-3) - 1e-3 / (1 - 2e-3), u, vv, d, 1e-3)[0])
        print("    delta=%-9g clamped %.8e   unclamped %.8e" % (d, a, b))

    print("\n%s DONE." % el())

if __name__ == "__main__":
    main()
