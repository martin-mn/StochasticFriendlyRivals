"""(2) Classification of stability on every stratum of the family p = (0, a, b, 1), all real (u,v).

For each of the 9 strata (open square, 4 open edges, 4 corners):
  * the eps->0 limit of self-play and of each of the 16 binary deviators is computed by
    Markov-chain-tree cofactors with the free coordinate(s) symbolic; lower eps-orders are checked
    to vanish identically and the leading coefficient of the normaliser is certified positive on
    the whole open stratum (cyclelib.limit_on_stratum), so each limit holds at EVERY point;
  * the 16 gains g_sigma = pi(sigma,p) - pi(p,p) are affine in (u,v).  With the claimed
    conditions L1 = a + 2v <= 0 (from ALLD) and L2 = b + 2u - 1 <= 0 (from ALLC), each gain is
    written g = l1 L1 + l2 L2 + mu, and l1, l2 >= 0, mu == 0 are certified on the stratum;
    ALLD and ALLC are certified to have l1 > 0 resp. l2 > 0 (necessity).
  * where that decomposition fails (it does on two corners), the exact stability region in the
    (u,v)-plane is computed as the polygon {g_sigma <= 0 for all 16 sigma}: its vertices and its
    recession cone, exactly.
"""
import sys, json, itertools
sys.path.insert(0, '.')
from cyclelib import *


def rat_sign(f, free):
    """Certified sign (+1, -1, 0) of rational function f on the open box; raises if mixed."""
    f = sp.factor(sp.cancel(f))
    if f == 0:
        return 0, '0'
    num, den = sp.fraction(f)
    sg = 1
    certs = []
    for part, pw in ((num, 1), (den, 1)):
        c, facs = sp.factor_list(part)
        sg *= (1 if c > 0 else -1)
        for g, m in facs:
            try:
                certify_positive(g, free); s = 1
            except ValueError:
                certify_positive(-g, free); s = -1
            sg *= s ** m
            certs.append('%s:%s' % (g, '+' if s > 0 else '-'))
    return sg, '; '.join(certs)


def polygon(gains):
    """Exact region {(u,v): g <= 0 for all g} for affine g.  Returns (vertices, recession rays)."""
    rows = []
    for g in gains:
        g = sp.expand(g)
        cu, cv, c0 = g.coeff(u), g.coeff(v), g.subs({u: 0, v: 0})
        if cu == 0 and cv == 0:
            if c0 > 0:
                return 'empty', []
            continue
        rows.append((cu, cv, c0))
    verts = set()
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(rows, 2):
        det = a1 * b2 - a2 * b1
        if det == 0:
            continue
        U = (-c1 * b2 + c2 * b1) / det; Vv = (-a1 * c2 + a2 * c1) / det
        if all(r[0] * U + r[1] * Vv + r[2] <= 0 for r in rows):
            verts.add((U, Vv))
    # recession cone: directions d with r[0]d_u + r[1]d_v <= 0 for all rows; test 8 + extreme dirs
    cand = set()
    for (a1, b1, _) in rows:
        cand.add((b1, -a1)); cand.add((-b1, a1))
    rays = sorted({d for d in cand if all(r[0] * d[0] + r[1] * d[1] <= 0 for r in rows)})
    return sorted(verts), rays


report = {}
ok_all = True
for name, (p, free) in STRATA.items():
    p = tuple(sp.Integer(x) if isinstance(x, int) else x for x in p)
    aa, bb = p[1], p[2]
    L1 = aa + 2 * v
    L2 = bb + 2 * u - 1
    Wp, kp, cp = limit_on_stratum(p, p, free)
    pip = payoff(Wp)
    rec = {'p': str(p), 'selfplay_W': [str(x) for x in Wp], 'selfplay_eps_order': kp,
           'selfplay_cert': cp, 'deviators': {}}
    gains = []
    decomp_ok = True
    for sig in BIN:
        W, k, c = limit_on_stratum(tuple(map(sp.Integer, sig)), p, free)
        g = sp.factor(sp.cancel(payoff(W) - pip))
        gains.append(g)
        ge = sp.expand(sp.cancel(g))
        l2 = sp.cancel(sp.diff(g, u) / 2)
        l1 = sp.cancel(sp.diff(g, v) / 2)
        mu = sp.cancel(g - l1 * L1 - l2 * L2)
        s1, c1 = rat_sign(l1, free) if not sp.diff(l1, u) and not sp.diff(l1, v) else (None, 'dep')
        s2, c2 = rat_sign(l2, free) if not sp.diff(l2, u) and not sp.diff(l2, v) else (None, 'dep')
        dok = (mu == 0) and s1 in (0, 1) and s2 in (0, 1)
        if sig == (0, 0, 0, 0):
            dok = dok and s1 == 1
        if sig == (1, 1, 1, 1):
            dok = dok and s2 == 1
        decomp_ok &= dok
        rec['deviators'][label(sig)] = {'W': [str(x) for x in W], 'eps_order': k, 'D_cert': c,
                                       'gain': str(g), 'l1': str(sp.factor(l1)), 'l2': str(sp.factor(l2)),
                                       'mu': str(sp.factor(mu)), 'sign_l1': s1, 'sign_l2': s2, 'ok': dok}
    rec['claimed_conditions_exact'] = decomp_ok
    print('==', name, '| self-play W =', Wp, '(eps-order of normaliser %d)' % kp)
    if decomp_ok:
        print('   stable  <=>  %s <= 0  and  %s <= 0   [all 16 gains = l1*L1 + l2*L2, l1,l2 >= 0 certified, mu == 0;'
              ' ALLD l1>0, ALLC l2>0]' % (L1, L2))
        for sig in BIN:
            d = rec['deviators'][label(sig)]
            print('     %s  l1 = %-28s l2 = %s' % (label(sig), d['l1'], d['l2']))
    else:
        if free:
            ok_all = False
            print('   DECOMPOSITION FAILS on a positive-dimensional stratum')
        verts, rays = polygon(gains)
        rec['true_region_vertices'] = [[str(x) for x in vv] for vv in verts]
        rec['true_region_rays'] = [[str(x) for x in r] for r in rays]
        print('   claimed region: %s <= 0, %s <= 0 ; TRUE region vertices %s, recession rays %s'
              % (L1, L2, verts, rays))
        for sig, g in zip(BIN, gains):
            print('     %s gain = %s' % (label(sig), g))
    report[name] = rec
json.dump(report, open('prove.json', 'w'), indent=1)
print('\npositive-dimensional strata all match the claim:', ok_all)
