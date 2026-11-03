#!/usr/bin/env python3
"""
fig2_check.py -- an independent receipt for panel c of Figure 2 (the fraction of the efficient set that is a partner).

    python3 fig2_check.py      ->  prints, and writes cache/fig2_check.txt

Panel c draws the closed form partners_volume.frac_of_efficient of PartnersRivalsM1f/check/working (the module that
PartnersRivalsM1f/Figures/fig_partners.py imports).  Here the same number is estimated by Monte Carlo from the RAW
limit test, with no use of the closed form or of the two-inequality criterion:

  * efficient strategies are drawn uniformly from the efficient set of the game's side of the switch line
    (friendly.md section 1): below, p = (1, p2, p3, p4) with (p2, p3, p4) uniform in the open cube (all in Eff_CC);
    above, (p1, p4) uniform in the open square on each of the two alternating sheets (p2, p3) = (0, 1) and (1, 0)
    (all in Eff_alt), half the draws on each;
  * each draw is confirmed efficient by the raw test too: deficit(1e-9) > -1e-5 or deficit(1e-9) > 0.2 deficit(1e-7)
    (a draw near a trap, e.g. p4 = 1e-4, has deficit O(eps/p4), shrinking 100-fold; on the sheet (1,0) the float
    chain at eps = 1e-9 carries a W_CD - W_DC round-off of ~3e-7, i.e. a deficit of ~1e-6 that does not shrink,
    which is why the absolute floor is 1e-5 and not the 1e-6 of the stability test);
  * it is called stable (limit Nash) by PartnersRivalsM1f/check/m1lib.nash_gap against the 16 deterministic
    deviators, with the two-epsilon classifier of PartnersRivalsM1f/Text/nash.md section 5:
    stable iff gap(1e-9) > -1e-6 or gap(1e-9) > 0.2 gap(1e-7).
The Monte Carlo fraction is compared with the closed form at nine games covering every branch of the closed form.
All imports from other projects are read-only (no bytecode is written).
"""
import os
import sys

sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M1F = "/Users/martin/Documents/CL1/PartnersRivalsM1f"
sys.path.insert(0, os.path.join(M1F, "check"))
sys.path.insert(0, os.path.join(M1F, "check", "working"))
import m1lib as L                                                  # noqa: E402
import partners_volume as PV                                       # noqa: E402

N = 2000
RNG = np.random.default_rng(20260926)
GAMES = [(-3.0, 0.7), (-6.0, 3.0), (0.3, 0.3), (2.5, -4.0), (-2.0, -3.0),      # below the switch line
         (6.0, -3.0), (-4.0, 6.0), (3.0, 3.0), (0.5, 2.0)]                      # above it


def stable(p, u, v):
    g1 = L.nash_gap(p, u, v, 1e-7)
    g2 = L.nash_gap(p, u, v, 1e-9)
    return g2 > -1e-6 or g2 > 0.2 * g1


def efficient(p, u, v):
    d1 = L.efficiency_deficit(p, u, v, 1e-7)
    d2 = L.efficiency_deficit(p, u, v, 1e-9)
    return d2 > -1e-5 or d2 > 0.2 * d1


def draws(u, v):
    if u + v < 1:
        q = RNG.random((N, 3))
        return np.column_stack([np.ones(N), q])
    q = RNG.random((N, 2))
    P = np.empty((N, 4))
    P[:, 0], P[:, 3] = q[:, 0], q[:, 1]
    P[: N // 2, 1], P[: N // 2, 2] = 0.0, 1.0
    P[N // 2:, 1], P[N // 2:, 2] = 1.0, 0.0
    return P


lines = []
worst = 0.0
for u, v in GAMES:
    P = draws(u, v)
    eff = np.array([efficient(p, u, v) for p in P])
    assert eff.all(), "a draw from the efficient set is not efficient at (%g, %g)" % (u, v)
    st = np.array([stable(p, u, v) for p in P])
    mc = st.mean()
    se = np.sqrt(max(mc * (1 - mc), 1.0 / N) / N)
    cf = float(PV.frac_of_efficient(np.array(u), np.array(v)))
    z = (mc - cf) / se
    worst = max(worst, abs(z))
    extra = ""
    if u + v > 1:
        extra = "   sheet (0,1): %.4f   sheet (1,0): %.4f" % (st[: N // 2].mean(), st[N // 2:].mean())
    s = ("(u,v) = (%5.1f,%5.1f)  %s  Monte Carlo %.4f +- %.4f   closed form %.4f   z = %+5.2f%s"
         % (u, v, "below" if u + v < 1 else "above", mc, se, cf, z, extra))
    print(s)
    lines.append(s)
s = "largest |z| over the %d games: %.2f  (N = %d draws per game)" % (len(GAMES), worst, N)
print(s)
lines.append(s)


def say(s):
    print(s)
    lines.append(s)


# ---- panel b: the two 3-dimensional carriers of the stable set, probed by the raw test (nash.md section 2).
#      dim 3 off Q = {u >= 1, v >= 0} and at (1, 0); on Q minus (1, 0) neither facet slab carries a stable point.
say("panel b: share of random points of the open facet slabs that are stable (raw test)")
NB = 600
for u, v in ((3.0, 3.0), (1.0, 2.0), (2.0, 0.0), (1.0, 0.0), (0.5, 3.0), (3.0, -1.0), (-2.0, 2.0)):
    A = np.column_stack([np.ones(NB), RNG.random((NB, 3))])                # {p1 = 1, p4 > 0}
    B = np.column_stack([RNG.random((NB, 3)), np.zeros(NB)])               # {p4 = 0, p1 < 1}
    fa = np.mean([stable(p, u, v) for p in A])
    fb = np.mean([stable(p, u, v) for p in B])
    q = (u >= 1 and v >= 0) and (u, v) != (1.0, 0.0)
    say("  (u,v) = (%4.1f,%4.1f)  %-22s  slab p1 = 1: %.4f   slab p4 = 0: %.4f"
        % (u, v, "in Q minus (1,0)" if q else ("the corner (1,0)" if (u, v) == (1.0, 0.0) else "off Q"), fa, fb))


# ---- panels d and e: the vertex census of the two prisms, by the raw rival margin against the 16 co-players,
#      and efficiency at a wedge-S game (2.5, -4), where the friendly rivals are the efficient members of Riv-
def rival(p, sign):
    m1 = L.rival_margin(p, sign, 1e-7)
    m2 = L.rival_margin(p, sign, 1e-9)
    return m2 > -1e-6 or m2 > 0.2 * m1


say("panels d, e: the binary vertices (raw test)")
for bits in ("0000", "1000", "0010", "1010", "0100", "1100", "1011", "1110", "1111", "1101"):
    p = np.array([float(b) for b in bits])
    say("  %s  competitive for T > S: %-5s  for T < S: %-5s  efficient at (2.5,-4): %s"
        % (bits, rival(p, +1), rival(p, -1), efficient(p, 2.5, -4.0)))
os.makedirs(os.path.join(HERE, "cache"), exist_ok=True)
with open(os.path.join(HERE, "cache", "fig2_check.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
