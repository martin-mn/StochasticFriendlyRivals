# `check/wsep/`: the two families of friendly rivals in wedge W are separated

This kit answers internal review, round 1, item 2. It is a computer-assisted proof that no friendly rival of the wedge W
(T > S, below the switch line) has p6 = p_{CD,DC} or p8 = p_{DC,CC} strictly between 0 and 1, and none has
p6 ≠ p8. The same machinery also fixes p9, p11, p14 on each family, and p4, p5, p13, p15 on all of W.
Everything is computed in the ε → 0 limit before any inequality is tested. Nothing is sampled. No tolerance is used.

## Statement certified

Let p be a friendly rival of W. Then:
- p0 = 1, p4 = 0, p13 = 0, p15 = 0 and p5 < 1;
- and exactly one of the following holds:
  - (i) p6 = p8 = 0 and p9 = p11 = p14 = 1;
  - (ii) p6 = p8 = 1, 0 < p9 < 1, p11 < 1 and p14 < 1.

Both cases occur: the binary friendly rival 19079 is in case (i), and the exact point eq:certW is in case (ii).
Both are re-certified by `lemma_check.py`.
- The friendly rivals of W therefore lie in the disjoint closed faces {p6 = p8 = 0} and {p6 = p8 = 1}, and meet both, so the set is **disconnected**.
- Case (ii) contains no binary strategy, and the closure of case (ii) contains no binary friendly rival.
- Whether each family is itself connected is **not** established. Of the 126 binary-anchored strata, 18 are not faces of its 5-cube.

## Why a stratum is ruled out (all exact)

A stratum sets each coordinate to 0, to 1, or to free in (0,1). One of the following rules it out:

1. **Not efficient.** Efficiency below the switch line means w_CC = 1 in self-play. By L1 this is decided by the support of the limit, which is an integer exponent computation and constant on the open stratum. p0 = 1 is necessary: if the limit support lies in the four states whose last outcome is CC, stationarity at ε = 0 forces it to be {CC,CC} with p0² = 1. The `c0` run confirms that none of the 28,697,814 strata with p0 < 1 is efficient.
2. **Beaten by ALLD.** ALLD cooperates with probability exactly ε, so w_DC → 0 and the margin is −w_CD ≤ 0. It is < 0 on the whole open stratum iff a state with last outcome CD lies in the limit support. That is again exact and constant on the stratum (L1).
3. **Beaten by one binary co-player q, by L3 + L4** (`cert/certify.py`). Every corner c of the closed stratum cube must have either E_c(q) > E(q), so that N_q(c) = 0, or a binary margin certified ≤ 0. The margin counts as ≤ 0 when it is strictly negative by the outward-rounded interval engine, a structural zero, or zero by integrality, since D*·margin is an integer. At least one corner must be strictly negative. In addition, the centre margin must be certified negative by an interval and is also computed exactly in rationals. Since N_q = Σ_c N_q(c) Π x^c (1−x)^(1−c) has positive weights on the open cube, N_q < 0 there, so q out-earns every point of the open stratum wherever T > S.

The only pruning used is p0 = 1. The SI's ALLD closed form (eq:si-alld) is not used as a lemma. It is cross-checked instead, with 0 disagreements on every efficient stratum of every class.

## Files

| file | role |
|---|---|
| `wsep.c` | Driver. It `#include`s `../cert/certeng.c` unchanged. It enumerates a class, decides efficiency and ALLD exactly, and searches for an excluding q: first a move-to-front pool, then all 65,536 co-players ordered by centre margin. |
| `cover.c` | For a list of co-players, which ones exclude each stratum. Used for the greedy set cover. |
| `verify.py` | Independent re-check of **every** certificate row with the x86_64 library path. It computes the exact rational centre margin (which must lie in the C enclosure, with equal E), recomputes the corner table with `en.corner_scan` (no early exit), and compares the counts. On 4000 random strata per class it re-decides efficiency and ALLD exactly. On 1500 certificates per class it tests 4 random rational interior points each, some 1/1000 from the faces; all must be strictly negative. |
| `exact_corners.py` | Fully rational re-check of sampled corner tables: 3000 certificates, 126,211 corners. |
| `corner_pairs.py`, `exact_table.py` | Fully rational re-check of **every** corner of **every** certificate. There are 42.45 M corner evaluations but only 583,817 distinct (binary corner, q) pairs. Each pair is evaluated once with `ex.lead_exact`, and then every row is checked against the resulting table. |
| `cover_prep.py`, `cover.c`, `cover_greedy.py` | Greedy set cover of the separation strata by the co-players found. |
| `lemma_check.py` | Exact L3 corner rule plus multilinear interpolation on 400 (stratum, q) pairs, a finite-ε convergence check, and the two positive controls. |
| `stats.py` | Writes `out/summary.json`. |
| `run_all.sh` | Reproduces everything, on 3 cores under `nice 10`. |
| `out/<class>.tsv` | One row per efficient stratum not excluded by ALLD: `c pattern d E q how nneg nzs nzi nhigh centre_lo centre_hi ntried`, where `how` is P (pool) or F (full search). |
| `out/<class>.exact.tsv` | The same rows with the exact centre margin, and `ok` from `verify.py`. |
| `out/extra/` | The p4, p5, p13, p15 runs and the p11/p14 family runs. |
| `ctrl_W364.txt`, `ctrl_W364.out` | Negative control: none of the 364 known friendly-rival strata of W is excluded. |

## Classes

All classes have p0 = 1. Templates use `.` for any value and `[..]` for a set.

| class | condition | strata |
|---|---|---|
| c6 | p6 free | 3^14 |
| c8 | p8 free, p6 ∈ {0,1} | 2·3^13 |
| c68 | (p6,p8) ∈ {(0,1),(1,0)} | 2·3^13 |
| c9a | p6 = p8 = 0, p9 ∈ {0, free} | 2·3^12 |
| c9b | p6 = p8 = 1, p9 ∈ {0,1} | 2·3^12 |
| x4 | p4 ∈ {1, free} | |
| x13 | p13 ∈ {1, free} | |
| x5 | p5 = 1 | |
| x15 | p15 ∈ {1, free} | |
| f1p11 | p6 = p8 = 0, p11 ∈ {0, free} | |
| f1p14 | p6 = p8 = 0, p14 ∈ {0, free} | |
| f2p11 | p6 = p8 = 1, p11 = 1 | |
| f2p14 | p6 = p8 = 1, p14 = 1 | |

Classes c6, c8 and c68 together are every stratum with (p6, p8) ∉ {(0,0), (1,1)}.

Every certificate row is verified three times: by the C interval engine, by `verify.py`, and exactly by `exact_table.py`.

Expected output of `stats.py`, for use as a regression check: every class ends with `NOT_EXCLUDED 0` and `lemma_violations 0`, and `verify.py` reports `rows_failing 0`.

```
./compile.sh && cc -O2 -ffp-contract=off -fno-fast-math -o cover cover.c -lm -lpthread
./run_all.sh          # about 5 min of C + python on 3 cores; the greedy cover adds about 8 min
```
