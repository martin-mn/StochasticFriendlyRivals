"""Remaining exact checks of the memory-one stability / partner statements.

 1. Quadrant u >= 1, v >= 0 minus (1,0): which strata carry 2-dimensional stable sets there
    (SI: 'the largest stable set is the set of alternators ... of dimension 2').
 2. The efficient edges {p_CC = 1, p_DD = 0, max(p_CD,p_DC) = 1} off the segments: partner status
    (not covered by Eq. (3), whose domain is p_DD > 0); compare with Eq. (3) at p_DD = 0.
 3. Efficiency closed forms (main Eq. (2), SI Eq. (S1)) against the exact self-play limit on every stratum.
 4. 'Every game has a partner' / the partner fraction is positive at every game: exact search on a grid
    including the switch line, T = S, u = 0, v = 0, u = 1 and (1,0).
 5. (1,1,0,1) at (3/2,-1/2): ALLC earns 7/6.
 6. Snowdrift games at which no binary memory-one strategy is stable exist.
Run: python3 audit3.py  (uses audit.sample_points; ~1 min on 12 processes)
"""
import itertools
import os
from fractions import Fraction as F
from multiprocessing import Pool

import m1x
import poly
from audit import STRATA, sample_points

Q_GAMES = [(F(1), F(1)), (F(1), F(3)), (F(2), F(0)), (F(5), F(0)), (F(2), F(1)), (F(3), F(2)),
           (F(2), F(1, 2)), (F(3), F(1)), (F(3, 2), F(1, 4)), (F(5), F(7)), (F(1), F(1, 1000)),
           (F(1001, 1000), F(0)), (F(1001, 1000), F(1, 1000)), (F(3, 2), F(1, 2))]


def eff_closed(p, side):
    """Main-text Eq. (2) / SI Eq. (S1). side: -1 below, 0 on, +1 above the switch line."""
    pc, a, b, c = p
    seg = (a, b) in ((0, 1), (1, 0))
    below = pc == 1 and not seg and (c > 0 or max(a, b) == 1)
    above = seg and pc < 1 and c > 0 and (pc, c) != (0, 1)
    if side < 0:
        return below
    if side > 0:
        return above
    return (pc == 1 and c > 0) or (pc == 1 and c == 0 and a > 0 and b > 0 and max(a, b) == 1) or above


def work(st):
    out = {'st': ''.join(st), 'quad': {}, 'effbad': [], 'edge': []}
    for p in sample_points(st, 7):
        V, G = m1x.gains(p)
        Tp = [(A, B, C) for (_, A, B, C, _) in G]
        # 1. quadrant games
        for g in Q_GAMES:
            if m1x.stable(G, *g):
                out['quad'].setdefault(str(g), []).append([str(x) for x in p])
        # 3. efficiency closed forms vs exact
        for side, game in ((-1, (F(0), F(0))), (0, (F(1, 2), F(1, 2))), (1, (F(1), F(1)))):
            ex = m1x.efficient(V, *game)
            # efficiency is a property of V alone on each side
            if ex != eff_closed(p, side):
                out['effbad'].append(([str(x) for x in p], side, [str(x) for x in V], ex))
        # 2. efficient p_DD = 0 edges
        pc, a, b, c = p
        if pc == 1 and c == 0 and (a, b) not in ((0, 1), (1, 0)) and max(a, b) == 1:
            eq3 = [(F(0), c, -(1 - a)), (1 - a, b, -(1 - a))]
            ok, bad, VT, VC = poly.compare(Tp, eq3)
            out['edge'].append(([str(x) for x in p], ok, [[str(x) for x in pt] for pt in VT],
                                [(k, [str(x) for x in pt]) for k, pt in bad[:4]]))
    return out


def partner_search(points_with_G, games):
    res = {}
    for g in games:
        found = None
        for p, V, G in points_with_G:
            if m1x.efficient(V, *g) and m1x.stable(G, *g):
                found = p
                break
        res[g] = found
    return res


if __name__ == '__main__':
    with Pool(int(os.environ.get('CERT_MAX_THREADS') or 12)) as pool:
        outs = pool.map(work, STRATA, chunksize=1)
    print('== 1. quadrant games: strata with stable sample points (stratum: games -> count)')
    for o in outs:
        if o['quad']:
            print(' ', o['st'], {g: len(ps) for g, ps in o['quad'].items()})
    print('== 3. efficiency closed form mismatches:', sum(len(o['effbad']) for o in outs))
    for o in outs:
        for r in o['effbad'][:3]:
            print('  ', o['st'], r)
    print('== 2. efficient p_DD = 0 edges vs Eq. (3) at p_DD = 0')
    for o in outs:
        for r in o['edge']:
            print('  ', o['st'], 'p', r[0], 'Eq3-equal' if r[1] else 'DIFFERENT', 'T-vertices', r[2][:6], r[3])
    # 4. partner at every game
    cands = []
    for st in STRATA:
        for p in sample_points(st, 7):
            if p[0] == 1 or (p[1], p[2]) in ((0, 1), (1, 0)):
                V, G = m1x.gains(p)
                cands.append((p, V, G))
    grid = [F(k, 4) for k in range(-16, 17)] + [F(-10), F(10), F(-100), F(100)]
    games = {(x, y) for x in grid for y in grid} | {(x, 1 - x) for x in grid} | {(x, x - 1) for x in grid} \
        | {(F(1), F(0))}
    res = partner_search(cands, sorted(games))
    miss = [g for g, f in res.items() if f is None]
    print('== 4. games without a partner among', len(cands), 'candidates:', len(miss), 'of', len(res), miss[:10])
    # 5.
    V, G = m1x.gains((F(1), F(1), F(0), F(1)))
    W = dict((n, w) for (n, A, B, C, w) in G)['1111']
    print('== 5. (1,1,0,1) at (3/2,-1/2): self-play', m1x.selfpay(V, F(3, 2), F(-1, 2)), 'ALLC earns',
          m1x.selfpay(W, F(3, 2), F(-1, 2)), 'efficient', m1x.efficient(V, F(3, 2), F(-1, 2)),
          'Eq3 I,II =', 1 * F(-1, 2) - 0, 0 + 0)
    # 6. Snowdrift: 0 < u < 1 < 1 + v
    bins = [m1x.gains(q) for q in m1x.BIN]
    sd = [(F(i, 10), F(j, 10)) for i in range(1, 10) for j in range(1, 60)]
    none = [g for g in sd if not any(m1x.stable(Gq, *g) for (_, Gq) in bins)]
    print('== 6. Snowdrift grid games (0<u<1, v>0) with no stable binary strategy:', len(none), 'of', len(sd),
          'e.g.', none[:5])
