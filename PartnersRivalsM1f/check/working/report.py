from fractions import Fraction as F
from itertools import product
import pickle

res = pickle.load(open('census.pkl', 'rb'))
STRATA = list(product([0, 'i', 1], repeat=4))
GU = [F(-3), F(-1), F(-1,4), F(0), F(1,4), F(1), F(3,2), F(3)]
GAMES = [(u, v) for u in GU for v in GU]


def name(spec):
    return ''.join('*' if c == 'i' else str(c) for c in spec)


full = {}          # (si, game) -> True if stratum is full-dimensional there
anyweak = {}
for si, k, pts in res:
    Z = set(range(16))
    for p, rows, Wp in pts:
        Z &= {j for j, r in enumerate(rows) if r == (0, 0, 0)}
    for g in GAMES:
        u, v = g
        f = w = False
        for p, rows, Wp in pts:
            st = all((r[0]*(u-1) + r[1]*v < r[2]) if k else (r[0]*(u-1) + r[1]*v <= r[2])
                     for j, r in enumerate(rows) if j not in Z)
            wk = all(r[0]*(u-1) + r[1]*v <= r[2] for r in rows)
            f |= st
            w |= wk
        full[(si, g)] = f
        anyweak[(si, g)] = w

# max dimension attained by a full-dimensional stratum piece
print('per game: max k over strata that are full-dimensional there, and the carriers')
hdr = {}
for g in GAMES:
    best = max(k for si, k, _ in res if full[(si, g)])
    carr = sorted(name(STRATA[si]) for si, k, _ in res if full[(si, g)] and k == best)
    hdr[g] = (best, carr)
# print a compact table
print('      v=' + '  '.join(f'{str(v):>6s}' for v in GU))
for u in GU:
    print(f'u={str(u):>5s} ' + '  '.join(f'{hdr[(u,v)][0]:>6d}' for v in GU))
print()
print('carriers of the maximum, by game:')
seen = {}
for g in GAMES:
    seen.setdefault((hdr[g][0], tuple(hdr[g][1])), []).append(g)
for (b, c), gs in sorted(seen.items()):
    print(f'  dim {b}  carriers {list(c)}')
    print(f'        at {len(gs)} games, e.g. {gs[:4]}')
print()
print('every stratum that is ever full-dimensional, with the games where it is:')
for si, k, _ in res:
    gs = [g for g in GAMES if full[(si, g)]]
    if gs:
        print(f'  {name(STRATA[si]):5s} k={k}  at {len(gs):3d}/64 games')
