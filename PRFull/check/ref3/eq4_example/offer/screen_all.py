"""Exact margins of p(x) against all 65536 binary co-players (leading-order exact arithmetic).
usage: python3 screen_all.py NUM DEN  -> writes margins_NUM_DEN.txt (code, margin) and prints a summary"""
import os, sys, time
from fractions import Fraction as F
from multiprocessing import Pool
from lo import code2q, pstar, margin

X = F(int(sys.argv[1]), int(sys.argv[2]))
P = pstar(X)


def work(rng):
    return [(c, margin(P, code2q(c))) for c in range(*rng)]


if __name__ == '__main__':
    t = time.time()
    chunks = [(i, i + 512) for i in range(0, 65536, 512)]
    with Pool(int(os.environ.get('CERT_MAX_THREADS') or 10)) as pool:
        res = [r for part in pool.map(work, chunks) for r in part]
    res.sort()
    with open('margins_%s_%s.txt' % (X.numerator, X.denominator), 'w') as f:
        for c, m in res:
            f.write('%d %s\n' % (c, m))
    ms = [m for _, m in res]
    mn = min(ms)
    neg = [c for c, m in res if m < 0]
    zero = [c for c, m in res if m == 0]
    print('x = %s: %d co-players, min margin %s attained by %d; negative: %d; zero: %d; %.1fs'
          % (X, len(res), mn, ms.count(mn), len(neg), len(zero), time.time() - t))
    if neg[:20]:
        print('  negative codes (first 20):', neg[:20])
    pos = sorted(m for m in ms if m > 0)
    print('  smallest positive margin:', pos[0], '=', float(pos[0]),
          ' attained by', [c for c, m in res if m == pos[0]][:10])
