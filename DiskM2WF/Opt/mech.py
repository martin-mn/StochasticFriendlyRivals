"""
The mechanism instrument: exact per-strategy tables, batched over all 65536
codes at once, and the linear functionals of an abundance vector built from
them.

WHY THIS EXISTS.  At u = 1e-2 the WF process does not select a strategy, it
maintains a cloud (PROJECT.md §6b: the reported leader was the rank-1 of 1 to 5
of 10 replicates, and at the square's (+8,+8) corner of none of them twice).  So
the deliverables measure the cloud, and the way to compare a WF cloud with
DiskM2IN's monomorphic-time vector -- two genuinely different objects, §1 -- is
to compare LINEAR FUNCTIONALS of them.  For any per-strategy table phi(s) that
is a property of the STRATEGY and not of the process,

    <phi> = sum_s pi(s) phi(s)

is the same functional of WF's abundance vector and of IN's stationary vector,
and is honestly comparable cell for cell even though the two `pi` mean
different things.  `Opt/exact.py`'s eps->0 masks are already such a table; this
module adds the ones that come from self play.

WHAT IS COMPUTED, all exactly and in one batched solve.

Every code is a deterministic-up-to-eps memory-two rule, so a player of code s
playing ITSELF is a 16-state Markov chain on k = 4a + b (a the last round, b
the round before, own action first, 0=CC 1=CD 2=DC 3=DD).  Its stationary
distribution PI16[s] gives, exactly:

  W[s]     the last-round outcome frequencies (w_CC, w_CD, w_DC, w_DD),
           W = PI16 summed over b.  Self-play payoff is cR*w_CC + cS*w_CD +
           cT*w_DC, which is `payf3`'s quantity and the thing DiskM2IN's `eff`
           column is built from.

  ALT[s]   PI16 at (last CD, before DC) + (last DC, before CD) -- the pair is
           TURN-TAKING, alternating who concedes.
  FRZ[s]   PI16 at (last CD, before CD) + (last DC, before DC) -- the pair is
           FROZEN in a fixed exploiter/victim role.

ALT and FRZ are the point.  Above the switch line u+v = 1 the attainable
optimum is Emax = (S+T)/2, reached by one-sided play, and a strategy reaches it
whether it turn-takes or freezes: both have the same LAST-ROUND frequencies, so
`w`, the payoff, `eff` and the outcome mixture are all identical and no
existing panel can tell them apart.  The two-round occupancy can.  That
distinction is exactly where the two processes part company (§6c).

    tau = <ALT> / (<ALT> + <FRZ>)      1 = turn-taking, 0 = frozen

    Hbar = 2 sum_k f_k (1 - f_k),  f_k = (ps_k - eps)/(1 - 2 eps)

Hbar is the mean Hamming distance between two independent draws from the
cloud, in genotype space, out of 16 bits -- the cloud's radius, where N_eff is
its size.  `ps_k` is what both codes already write out, with the identical
formula, so Hbar is available on both sides of the comparison.

The tables depend only on eps, so they are computed once and cached in
`../Data/mech_e<n>.npz`.

Validated (see `python3 mech.py`): against `selfplay.py`, which does the same
16x16 solve one code at a time by a different route (lstsq on the augmented
system rather than a replaced row), and against DiskM2IN's own `eff` column.
"""
import os
import numpy as np

NC = 65536
NS = 16
MIRROR = np.array([0, 2, 1, 3])
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Data")

# the four two-round occupancies that name a one-sided pair's mode.
# k = 4a + b, a = last round, b = the round before; 1 = CD (I cooperated,
# the co-player defected), 2 = DC.
K_CD_DC, K_DC_CD = 4 * 1 + 2, 4 * 2 + 1      # alternating: the roles swapped
K_CD_CD, K_DC_DC = 4 * 1 + 1, 4 * 2 + 2      # frozen: the roles held


def _build(eps):
    """PI16 for every one of the 65536 codes, in one batched solve."""
    pmax, pmin = 1.0 - eps, eps
    codes = np.arange(NC, dtype=np.int64)
    bits = ((codes[:, None] >> np.arange(NS)[None, :]) & 1).astype(bool)
    pc = np.where(bits, pmax, pmin)                       # (NC,16) P(cooperate)

    # in state k = 4a+b the co-player is in the mirrored state
    KM = np.array([4 * MIRROR[k // 4] + MIRROR[k % 4] for k in range(NS)])
    p1 = pc                                               # mine
    p2 = pc[:, KM]                                        # theirs

    P = np.zeros((NC, NS, NS))
    pr = (p1 * p2, p1 * (1.0 - p2), (1.0 - p1) * p2, (1.0 - p1) * (1.0 - p2))
    for k in range(NS):
        a = k // 4                                        # this round's a is
        for x in range(4):                                # next state's b
            P[:, k, 4 * x + a] += pr[x][:, k]

    # pi P = pi with sum pi = 1: (P^T - I) pi = 0, last row replaced by the
    # normalisation.  Every transition has probability >= eps^2 > 0, so the
    # chain is irreducible and the solution is unique.
    M = np.transpose(P, (0, 2, 1)) - np.eye(NS)[None, :, :]
    M[:, -1, :] = 1.0
    rhs = np.zeros((NC, NS))
    rhs[:, -1] = 1.0
    pi16 = np.linalg.solve(M, rhs)
    resid = np.abs(np.einsum("nij,ni->nj", P, pi16) - pi16).max()
    return pi16, resid


def tables(eps=1e-4, cache=True):
    """{PI16, W, ALT, FRZ, EFFPAY} for every code, cached by eps."""
    tag = "e%d" % round(-np.log10(eps))
    path = os.path.normpath(os.path.join(DATA, "mech_%s.npz" % tag))
    if cache and os.path.exists(path):
        z = np.load(path)
        return {k: z[k] for k in z.files}
    pi16, resid = _build(eps)
    T = dict(PI16=pi16,
             W=pi16.reshape(NC, 4, 4).sum(axis=2),
             ALT=pi16[:, K_CD_DC] + pi16[:, K_DC_CD],
             FRZ=pi16[:, K_CD_CD] + pi16[:, K_DC_DC],
             eps=np.array(eps), resid=np.array(resid))
    if cache:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        np.savez_compressed(path, **T)
    return T


# ------------------------------------------------------------------ per cell
def selfpay(T, cR, cS, cT):
    """each code's own self-play payoff at this game (65536-vector)"""
    W = T["W"]
    return cR * W[:, 0] + cS * W[:, 1] + cT * W[:, 2]


def hbar(ps, eps=1e-4):
    """mean Hamming distance, out of 16 bits, between two independent draws.

    ps_k is the cloud's probability of playing C in state k, which is
    eps + (1-2eps) f_k with f_k the abundance-weighted fraction of the cloud
    whose bit k is set -- both codes write ps with that identical formula.
    Two independent draws differ at bit k with probability 2 f_k (1 - f_k).
    """
    f = (np.asarray(ps, dtype=float) - eps) / (1.0 - 2.0 * eps)
    f = np.clip(f, 0.0, 1.0)
    return float(2.0 * (f * (1.0 - f)).sum())


def cell(g, T, d=None):
    """Every mechanism functional of one cell, from its abundance vector.

    `g` is a game dict as wfload (or the IN loader) returns it: it needs
    `pi`, `cR`, `cS`, `cT`, `emax`, `ps`, and `eff` if the realised payoff is
    to be compared with the self-play one.  `d` is exact.load(), optional; when
    given the eps->0 mass decomposition is added too.
    """
    p = np.asarray(g["pi"], dtype=float)
    s = p.sum()
    if s > 0:
        p = p / s
    W = p @ T["W"]                                   # (4,) the outcome mixture
    alt, frz = float(p @ T["ALT"]), float(p @ T["FRZ"])
    one = float(W[1] + W[2])                         # one-sided mass
    sp = float(p @ selfpay(T, g["cR"], g["cS"], g["cT"]))
    out = dict(W=W, wcc=float(W[0]), wdd=float(W[3]), onesided=one,
               alt=alt, frz=frz,
               tau=(alt / (alt + frz)) if (alt + frz) > 0 else np.nan,
               eff_self=sp / g["emax"], selfpay=sp,
               hbar=hbar(g["ps"], float(T["eps"])) if "ps" in g else np.nan,
               coop=float(W[0] + W[1]))
    if "eff" in g:
        # the ONLY cross-play observable: sf.f measured what the population
        # actually earned against itself as a MIXTURE, where eff_self is what
        # it would earn if every player met its own clone.  The difference is
        # assortment, and IN cannot have it -- its population is monomorphic
        # by construction, so for IN the two coincide identically.
        out["gap"] = float(g["eff"]) - out["eff_self"]
    if d is not None:
        import exact
        meff, riv, ne = exact.masks(d, g["u"], g["v"])
        partner = meff & ne
        out.update(mass_eff=float(p[meff].sum()),
                   mass_riv=float(p[riv].sum()),
                   mass_ne=float(p[ne].sum()),
                   mass_partner=float(p[partner].sum()),
                   mass_fr=float(p[partner & riv].sum()),
                   nfr=int((partner & riv).sum()))
    return out


# --------------------------------------------------------- behaviour classes
def classes(T, dp=2, key="PI16"):
    """Group the 65536 codes by what they DO, not by their 16 bits.

    Two codes are the same BEHAVIOUR if their exact self-play two-round
    occupancy agrees to `dp` decimals.  At dp=3 that is 1853 classes over the
    65536 codes.

    This is not a convenience -- it is the finding.  In genotype space a WF
    cell looks irreproducible: a median of 8 distinct rank-1 strategies over
    its 10 replicates, all 10 different in 75 of the 455 cells, and TVD 0.25
    between two independent runs of the same cell.  Project onto behaviour and
    the same data is sharp: a median of ONE distinct behaviour per cell
    (exactly one in 395 of 455), 64 of those 75 cells naming a single
    behaviour, and a reproducibility ledger of TVD 0.0013 with the leader
    agreeing 19 times out of 19.  The cloud is a NEUTRAL NETWORK -- many
    genotypes, one behaviour -- and the process selects the behaviour.

    `dp` is the only free parameter and the conclusion is insensitive to it
    across the window the data support, but the DEFAULT IS 2 AND NOT 3, which
    was found by building the per-game tables.  At 3 dp the near-pure
    mutual-cooperation behaviour is split in two by the rounding boundary at
    0.9995: two classes whose PI16 differ by 0.0002, separated only by where a
    rare error excursion goes.  That is an artefact of decimal rounding, not a
    difference in behaviour.  At 2 dp they merge, and every statistic improves
    -- the 10 replicates agree unanimously in 429 of 455 cells rather than 395,
    and WF and IN name the same leading behaviour in 315 rather than 307.

    At 2 dp the whole 455-game plane is spanned by exactly THREE leading
    behaviours: near-pure mutual cooperation (274 games, 7639 codes), one-sided
    play that TURN-TAKES (173 games, 1216 codes), and one-sided play that is
    FROZEN into a fixed exploiter/victim pair (8 games, the same 1216 codes).
    The last two have identical last-round frequencies (0, 0.5, 0.5, 0) and are
    separated only by the two-round occupancy -- which is what `tau` measures.

    At 6 dp the definition resolves outcome frequencies to 1e-6, four orders
    below the ledger's own 3e-4 floor on efficiency, and is fitting noise.
    """
    K = np.round(T[key], dp)
    _, inv = np.unique(K, axis=0, return_inverse=True)
    return inv, int(inv.max()) + 1


def project(p, inv, nclass):
    """an abundance vector over codes -> the same mass over behaviours"""
    q = np.bincount(inv, weights=np.asarray(p, dtype=float), minlength=nclass)
    s = q.sum()
    return q / s if s > 0 else q


if __name__ == "__main__":
    import sys
    import selfplay

    eps = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-4
    T = tables(eps, cache=False)
    print("batched solve over %d codes at eps = %g" % (NC, eps))
    print("  max |pi P - pi| = %.2e" % float(T["resid"]))

    # 1. against selfplay.py, which solves the same chain a different way
    worst = 0.0
    for c in (0, 1, 23175, 22663, 55394, 47588, 44844, 53036, 65535, 4063):
        w = selfplay.selfplay_freqs(c, eps)
        worst = max(worst, float(np.abs(w - T["W"][c]).max()))
    print("  vs selfplay.selfplay_freqs on 10 codes: max diff %.2e" % worst)

    # 2. the turn-taking split on the codes the two processes disagree about
    print("\n  code       w_CC     w_CD     w_DC     w_DD      ALT      FRZ"
          "     tau")
    for c, nm in ((55394, "IN's absorbing state"),
                  (47588, "its mirror"),
                  (44844, "a WF far-field leader"),
                  (53036, "another"),
                  (56108, "another"),
                  (23175, "the c=0.65 ray leader"),
                  (65535, "AllC"), (0, "AllD")):
        w, a, f = T["W"][c], T["ALT"][c], T["FRZ"][c]
        t = a / (a + f) if (a + f) > 0 else float("nan")
        print("  %5d  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %6.4f   %s"
              % (c, w[0], w[1], w[2], w[3], a, f, t, nm))
