"""All 81 strata of the cube (each p_i at 0, interior, or 1): is the limit-Nash set
FULL-DIMENSIONAL inside the stratum, at each of 64 games?

For a stratum of dimension k the Nash set has dimension k iff some relative-interior point
satisfies every non-identically-vanishing row STRICTLY (the rows are affine in (u,v), so the
strict set is open).  Rows are computed once per sampled point in exact arithmetic and then
evaluated at every game for free."""
from fractions import Fraction as F
from itertools import product
import random, pickle, sys
from multiprocessing import Pool
from limlib import nash_rows

NS = {0: 1, 1: 40, 2: 120, 3: 300, 4: 600}
STRATA = list(product([0, 'i', 1], repeat=4))


def sample(spec, rng):
    return tuple(F(rng.randint(1, 199), 200) if c == 'i' else F(c) for c in spec)


def work(job):
    si, spec = job
    k = sum(1 for c in spec if c == 'i')
    rng = random.Random(1000 + si)
    pts = []
    for _ in range(NS[k]):
        p = sample(spec, rng)
        rows, Wp = nash_rows(p)
        pts.append((p, rows, tuple(Wp)))
    return si, k, pts


if __name__ == '__main__':
    with Pool(12) as pool:
        res = pool.map(work, list(enumerate(STRATA)))
    with open('census.pkl', 'wb') as f:
        pickle.dump(res, f)
    print('done', sum(len(r[2]) for r in res), 'points')
