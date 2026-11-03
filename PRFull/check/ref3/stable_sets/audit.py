"""Systematic exact audit of the memory-one stability / partner statements of ms.tex on all 81 strata.

For every sample point p of every stratum (each coordinate 0, open, 1) the exact eps->0 gains of the 16
binary deviators are affine in (u,v), so the true stable set of games is a polygon T(p).  Each claim of
the manuscript that covers p is turned into a polygon Cl(p) on a domain D (e.g. u+v >= 1), and the two
are compared EXACTLY inside the box |u|,|v| <= 1000 by testing every vertex of each against the other
(poly.compare).  In addition every p is tested on a fixed grid of games (incl. the switch line, T=S,
u=0, v=0, u=1, (1,0)) and at its own ZD-equaliser game.

Run: python3 audit.py  (12 processes, ~1 min).  Writes audit_out.json and prints a summary.
"""
import json
import os
import random
import sys
from fractions import Fraction as F
from itertools import product
from multiprocessing import Pool

import m1x
import poly

STRATA = list(product('0o1', repeat=4))
OV = [F(1, 2), F(1, 3), F(3, 4)]
EXT = [F(1, 1000), F(999, 1000)]


def sample_points(st, seed):
    rng = random.Random(seed)
    free = [i for i, c in enumerate(st) if c == 'o']
    base = [F(0) if c == '0' else F(1) if c == '1' else None for c in st]
    pts = set()
    for vals in product(OV, repeat=len(free)):
        q = list(base)
        for i, x in zip(free, vals):
            q[i] = x
        pts.add(tuple(q))
    for vals in product(EXT, repeat=len(free)):
        q = list(base)
        for i, x in zip(free, vals):
            q[i] = x
        pts.add(tuple(q))
    for _ in range(8 if free else 0):
        q = list(base)
        for i in free:
            d = rng.randint(7, 997)
            q[i] = F(rng.randint(1, d - 1), d)
        pts.add(tuple(q))
    return sorted(pts)


def U(a, b, c):
    return (F(a), F(b), F(c))


LINE = [U(1, 1, -1), U(-1, -1, 1)]          # u + v = 1
ABOVE = [U(-1, -1, 1)]                      # u + v >= 1


def claims(p):
    """List of (id, domain, claimed half-planes, kind) covering p.  kind 'stable' or 'partner'."""
    pc, a, b, c = p
    out = []
    exc = (a, b) in ((0, 1), (1, 0))
    # A: Eq. (3) on {p_CC = 1, p_DD > 0} off the segments (proved earlier today; regression check)
    if pc == 1 and c > 0 and not exc:
        out.append(('A_eq3', [], [U(0, c, -(1 - a)), U(1 - a, b, -(1 - a))], 'stable'))
    # B: {p_DD = 0, p_CC < 1} off its edges
    if c == 0 and pc < 1 and not exc:
        out.append(('B_pdd0', [], [U(1 - pc, 0, b), U(1 - a, b, b)], 'stable'))
    # C: alternators above (and on) the switch line
    if (a, b) == (0, 1) and pc < 1 and c > 0 and (pc, c) != (0, 1):
        out.append(('C_alt_above', ABOVE, [U(-c - 1, c - 1, c - 1), U(-pc, pc - 2, pc)], 'stable'))
    # C1: (1,0,1,x), x > 0, on the switch line: partner iff x(u-1) >= -1
    if (pc, a, b) == (1, 0, 1) and c > 0:
        out.append(('C1_seg101x_line', LINE, [U(-c, 0, c - 1)], 'partner'))
        out.append(('C1b_seg101x_altconds', LINE, [U(-c - 1, c - 1, c - 1), U(-pc, pc - 2, pc)], 'stable'))
    # D: frozen alternators above the switch line: stable iff T = S (u - v = 1)
    if (a, b) == (1, 0) and pc < 1 and c > 0 and (pc, c) != (0, 1):
        out.append(('D_frozen_above', ABOVE, [U(1, -1, -1), U(-1, 1, 1)], 'stable'))
    # D1: (1,1,0,x), x > 0, on the switch line: partner only at (1,0)
    if (pc, a, b) == (1, 1, 0) and c > 0:
        out.append(('D1_seg110x_line', LINE, [U(1, 0, -1), U(-1, 0, 1)], 'partner'))
    # E: interior: stable iff equaliser (ZD sheet or constant u/(u+v))
    if all(0 < x < 1 for x in p):
        if pc == a == b == c:
            x = pc
            out.append(('E_interior_const', [], [U(1 - x, -x, 0), U(x - 1, x, 0)], 'stable'))
        else:
            us, vs = zd_game(p)
            out.append(('E_interior_ZD', [], [U(1, 0, -us), U(-1, 0, us), U(0, 1, -vs), U(0, -1, vs)], 'stable'))
    return out


def zd_game(p):
    """The game at which p is a ZD equaliser (co-player payoff vector (1,1+v,u,0)), or None."""
    pc, a, b, c = p
    al = pc - 1 - c
    if al == 0:
        return None
    return (b - c) / al, (a - 1 - c) / al - 1


GRID_U = [F(x) for x in (-5, -3, -2, -1)] + [F(k, 4) for k in range(-3, 9)] + [F(x) for x in (3, 5)]
BASE_GAMES = sorted({(u, v) for u in GRID_U for v in GRID_U}
                    | {(u, 1 - u) for u in GRID_U} | {(u, u - 1) for u in GRID_U}
                    | {(u, -u) for u in GRID_U} | {(F(1), F(0))})


def closed_info(p):
    pc, a, b, c = p
    cl = []
    if pc == 1:
        cl.append('CC')
    if c == 0:
        cl.append('DD')
    if (a, b) == (0, 1):
        cl.append('alt')
    if (a, b) == (1, 0):
        cl.append('frozen')
    if (pc, c) == (0, 1):
        cl.append('CC-DD')
    return cl


def work(st):
    seed = sum('0o1'.index(ch) * 3 ** i for i, ch in enumerate(st)) + 12345
    rows = []
    for p in sample_points(st, seed):
        V, G = m1x.gains(p)
        T = [(A, B, C) for (_, A, B, C, _) in G]
        VT = poly.vertices(T)
        dimT = poly.affine_dim(VT)
        binding = [s for s in range(4) if V[s] > 0 and 0 < p[s] < 1]
        zg = zd_game(p)
        res = {'p': [str(x) for x in p], 'V': [str(x) for x in V], 'closed': closed_info(p),
               'binding': binding, 'dimT': dimT,
               'VT': [[str(x) for x in pt] for pt in VT] if dimT <= 1 else None,
               'zd': [str(x) for x in zg] if zg else None, 'claims': []}
        # stable at its ZD game?
        if zg and abs(zg[0]) <= 1000 and abs(zg[1]) <= 1000:
            res['stable_at_zd'] = m1x.stable(G, *zg)
        for cid, D, Cl, kind in claims(p):
            ok, bad, VTD, VCD = poly.compare(T, Cl, D)
            # grid check as an independent second test
            gbad = []
            for (u, v) in BASE_GAMES + ([zg] if zg else []):
                if not poly.feasible(D, (u, v)):
                    continue
                if m1x.stable(G, u, v) != poly.feasible(Cl, (u, v)):
                    gbad.append((u, v))
            effbad = []
            if kind == 'partner':
                for pt in VCD:
                    if not m1x.efficient(V, *pt):
                        effbad.append(pt)
            item = {'id': cid, 'ok': ok and not gbad, 'n_bad': len(bad), 'n_gridbad': len(gbad),
                    'eff_bad': [[str(x) for x in pt] for pt in effbad]}
            if bad or gbad:
                item['bad'] = [(k, [str(x) for x in pt], [(n, str(g)) for n, g in m1x.gainers(G, *pt)])
                               for k, pt in bad[:6]]
                item['gridbad'] = [([str(x) for x in pt], [(n, str(g)) for n, g in m1x.gainers(G, *pt)])
                                   for pt in gbad[:6]]
                item['VT'] = [[str(x) for x in pt] for pt in VTD[:12]]
            res['claims'].append(item)
        rows.append(res)
    return ''.join(st), rows


if __name__ == '__main__':
    with Pool(int(os.environ.get('CERT_MAX_THREADS') or 12)) as pool:
        out = dict(pool.map(work, STRATA, chunksize=1))
    json.dump(out, open('audit_out.json', 'w'), indent=1)
    # summary
    tally = {}
    for st, rows in out.items():
        for r in rows:
            for c in r['claims']:
                t = tally.setdefault(c['id'], [0, 0, set()])
                t[0] += 1
                if not c['ok'] or c['eff_bad']:
                    t[1] += 1
                    t[2].add(st)
    print('claim                    tested  failed  strata-with-failures')
    for k, (n, f, s) in sorted(tally.items()):
        print('%-24s %6d  %6d  %s' % (k, n, f, sorted(s)))
    print('\nper-stratum: n points, self-play support, binding states, dim of T (0..2, -1 empty) counts')
    for st, rows in out.items():
        sup = sorted({tuple(i for i in range(4) if F(r['V'][i]) > 0) for r in rows})
        bind = sorted({tuple(r['binding']) for r in rows})
        dims = {}
        for r in rows:
            dims[r['dimT']] = dims.get(r['dimT'], 0) + 1
        print(st, len(rows), 'supp', sup, 'bind', bind, 'dimT', dict(sorted(dims.items())),
              'closed', rows[0]['closed'])
