#!/usr/bin/env python3
"""The SIXTEEN binary strategies' partner field on the same square, for the comparison.

Recomputed here from this project's own machinery -- the 16 facet triples of each vertex
resident and the efficiency predicates of m1lib -- rather than read off
../../PartnersRivalsM1, so that the two shares in the answer are like for like.
"""
import numpy as np

import partners_forms as F
from m1lib import DET, genome, in_eff_alt, in_eff_cc, in_eff_line

SQ, NG = 8.0, 3201


def main():
    V = np.array(DET)
    a, b, c = F._bank(V)                       # (16 deviators, 16 residents)
    g = (np.arange(NG) + 0.5) / NG * (2 * SQ) - SQ
    U, W = (x.ravel() for x in np.meshgrid(g, g, indexing="ij"))
    below, above = U + W < 1, U + W > 1
    cnt = np.zeros(U.size, int)
    who = []
    for j in range(16):
        p = DET[j]
        nash = (c[:, j][:, None] - a[:, j][:, None] * (U - 1)[None, :]
                - b[:, j][:, None] * W[None, :]).min(axis=0) >= -1e-9
        eff = np.where(below, in_eff_cc(p), np.where(above, in_eff_alt(p), in_eff_line(p)))
        part = nash & eff
        cnt += part
        if part.any():
            who.append((genome(p), 100 * part.mean()))
    print("the 16 binary strategies on -8 <= u,v <= 8   (%d x %d cell-centred grid)"
          % (NG, NG))
    for k in range(int(cnt.max()) + 1):
        if np.any(cnt == k):
            print("    %d partner(s) : %8.4f %%" % (k, 100 * np.mean(cnt == k)))
    print("    NO partner    : %8.4f %%   (the published figure's 56.25 %%)"
          % (100 * np.mean(cnt == 0)))
    print("    a partner     : %8.4f %%   (the published figure's 43.75 %%)"
          % (100 * np.mean(cnt > 0)))
    top = cnt.max()
    m3 = (cnt == 3)
    print("    the maximum 3 is attained on %.4f %% and coincides with {u<=1, v<=0}: %s"
          % (100 * m3.mean(), bool(np.array_equal(m3, (U <= 1) & (W <= 0)))))
    print("    counts above 3 occur on %.6f %% of the square (measure-zero rays: 1011 on"
          " v = (u-1)/3 and 1101, as ../../PartnersRivalsM1 records)"
          % (100 * np.mean(cnt > 3)))
    print("    strategies that are a partner somewhere: %s"
          % [(gname, round(sh, 4)) for gname, sh in who])
    print("\n    CONTINUUM: a partner on 100.0000 % of the same square, and the")
    print("    saturation region 'every efficient strategy is a partner' is")
    print("    {u<=1, v<=0} below the line -- EXACTLY the 16's maximum-count region --")
    print("    together with {u>=0, v>=0} above it: 135.5/256 = 52.9297 %.")


if __name__ == "__main__":
    main()
