"""Self-play of p(x) after one error at <CC,CC>, eps = 0 afterwards (no further error), exact.
Tracks: rounds, times A (the erring player) and B are exploited, offers (first CC), and checks balance."""
from fractions import Fraction as F
from lo import pstar, M
NAMES = ['CC', 'CD', 'DC', 'DD']

def run(x, depth=200):
    p = pstar(x)
    # state: (sA, exploitedA, exploitedB, firstCCs); sA = A's view; B's view is M[sA]
    # round 1: A defects by error at <CC,CC>, B cooperates: outcome DC for A
    from collections import defaultdict
    dist = {(4 * 2 + 0, 0, 1, 0): F(1)}       # A at <DC,CC>; B exploited once (the error)
    done = defaultdict(F)                       # (eA, eB, firstCCs) -> prob at restoration of <CC,CC>
    ER = F(0); rounds = F(0); mass = F(1)
    for t in range(depth):
        nd = defaultdict(F)
        for (s, eA, eB, k), pr in dist.items():
            a = s // 4
            pa, pb = F(p[s]), F(p[M[s]])
            for ca, qa in ((1, pa), (0, 1 - pa)):
                for cb, qb in ((1, pb), (0, 1 - pb)):
                    w = pr * qa * qb
                    if w == 0: continue
                    o = 2 * (1 - ca) + (1 - cb)
                    ns = 4 * o + a
                    eA2 = eA + (o == 1); eB2 = eB + (o == 2)
                    k2 = k + (o == 0 and a != 0)       # a CC round not preceded by CC
                    rounds += w
                    if ns == 0:
                        done[(eA2, eB2, k2)] += w
                    else:
                        nd[(ns, eA2, eB2, k2)] += w
        dist = nd
    left = sum(dist.values())
    return done, left, rounds

for x in [F(1, 4), F(1, 7), F(5, 13)]:
    done, left, rounds = run(x)
    assert all(eA == eB for (eA, eB, k) in done), 'unbalanced path'
    EA = sum(pr * eA for (eA, eB, k), pr in done.items())
    Ek = sum(pr * k for (eA, eB, k), pr in done.items())
    print('x = %s: mass not yet recovered after 200 rounds %.2e; E[times each exploited, counting the error] ~ %.10f;'
          ' E[CC rounds not preceded by CC] ~ %.10f; E[rounds after the error round] ~ %.10f' % (x, float(left), float(EA), float(Ek), float(rounds)))
    # exact geometric values: success per decision = x/2
    print('   exact: 1 + (1 - x/2)/(x/2) =', 1 + (1 - x / 2) / (x / 2), '; offers before recovery = x * 2/x =', 2)
