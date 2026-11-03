# `check/cert/`: certificates for memory-two friendly rivals (ε → 0)

This library answers the internal reviewer's request for certificates in place of sampling: boxes of positive
volume, whole open strata, and strata excluded outright. Everything is computed in the ε → 0 limit, before
any inequality is tested.

```python
import sys; sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/check')
from cert import engine as en, exact as ex, certify as ce
```

To build, run `./compile.sh`. On macOS this makes a universal `libcert.dylib`; on Linux/Cannon it makes `libcert.so`.
To test, run `cd PRFull/check && nice -n 10 python3 cert/tests.py`. The full suite (sections `a b c d e t`) takes about 4 minutes on 3 cores.

Local rule: C threads and exact-arithmetic pools are capped at `engine.MAX_THREADS = 3`. On Cannon, override this with `CERT_MAX_THREADS`.

## Conventions

These match PartnersRivalsM2f.
- `p_s` is the probability of cooperating at state s = 4x + y, where x is the last outcome and y the one before. The outcomes are CC=0, CD=1, DC=2, DD=3, with the focal player's action first.
- Errors act as p → ε + (1−2ε)p on both players.
- A binary co-player is a genome g in 0..65535 whose bit s is q_s (1 = C).
- A pattern is `'1*1*001*1*1**0*0'`, where character s is p_s. It can also be an int8 array of 0, 1 and 2 (= free), or a census code (base-3 digits 0 = 0, 1 = free, 2 = one).
- The margin is sign·(w_DC − w_CD), with sign = +1 where T > S.
- `WEDGES` maps each wedge to its sign and required efficiency class: W → (+1,'C'), S → (−1,'C'), N → (+1,'A'), E → (−1,'A'), ray+ → (+1,'Z'), ray− → (−1,'Z'). Class C means w_CC = 1, A means w_CD + w_DC = 1, and Z means w_DD = 0.

## The structural lemma

Status: verified; see `certify.py` for the statement and tests (c) for the checks. Fix a binary co-player q and a stratum. Then:
- The limit is w_k = A\*_k(x)/D\*(x).
- A\*_k is the ε^E coefficient of the tree sum rooted at k, where E is the least tree order. It is multilinear and nonnegative in the basis {x_s, 1−x_s}.
- The margin therefore has the sign of the multilinear polynomial N_q = D\*·margin.

Four consequences follow:
- **L1.** The support is constant on an open stratum.
- **L2.** On a box inside the open stratum, the minimum is attained at a vertex.
- **L3.** At a corner, N_q(c) = D\*(c)·(the corner's own binary margin) if E_c = E. If E_c > E, then N_q(c) = D\*(c) = 0. E_c < E never occurs.
- **L4.** If N_q ≤ 0 at every corner and N_q < 0 at one interior point, then N_q < 0 on the whole open cube. In the other direction, N_q ≥ 0 at every corner gives margin ≥ 0 on the whole open stratum.

How E is obtained: the GTH pivots multiply to A_0(ε), by the matrix-tree theorem. So E_0 = Σ_k exp(pivot_k), and E = E_0 + min_k exp(π_k/π_0).

**Caveats found**
1. **The "1/37 gap" is false.** The claim that the smallest non-zero |w_DC − w_CD| over binary pairs is 1/37 (internal review notes, not deposited; an earlier draft of the paper; M2f methods.tex:76) does not hold. The pair p = 15204, q = 21477 gives exactly 1/963, and 637 of 10⁵ random pairs fall below 1/37. The library uses **integrality** instead:
   - For a binary pair every edge has leading coefficient 1, so A\*_k is the number of minimum-order spanning trees rooted at k.
   - D\*·(w_DC − w_CD) is therefore an integer, and an enclosure with D\*_hi·max|m| < 1 proves the margin is 0.
   - The de Bruijn graph B(4,2) has 4¹³ spanning arborescences per root, 2³⁰ in all. So every non-zero binary margin is at least 2⁻³⁰ ≈ 9.3·10⁻¹⁰, which is *below* the census tolerance of 10⁻⁹. The largest D\* in the 10⁵ sample was 1428.
2. **Self-play is not multilinear.** Both x_k and x_mirror(k) enter, and x_s can appear twice in one tree. Efficiency is still exact, by L1, because it depends only on the support.
3. **The co-player quantifier.** Every certificate here quantifies over the 65 536 binary co-players. That this covers all of [0,1]¹⁶ is the paper's box-vertex lemma. It also follows from the same lemma in q, since the co-player factor out of state k uses q_mirror(k) and mirror is a bijection.

## API

| function | returns | what it proves / computes |
|---|---|---|
| `en.f_eval(p, q)` | w16, W, Ek, E, D, supp | Float limit. **E, Ek and supp are exact integers**; D ≈ D\*. |
| `en.f_scan(p, qlist=None, nthreads)` | W (n×4), E, supp, D | The same over all 65 536 (or listed) co-players. |
| `en.i_eval(plo, phi, q)` | wlo/whi, Wlo/Whi, Ek, E, Dlo/Dhi, supp | Outward-rounded enclosures, valid for **every** x in the box. Free coordinates must lie in (0,1). |
| `en.i_scan(plo, phi, sign, qlist, nthreads)` | cls, mlo, mhi, E, supp | Per co-player: cls = +1 certified > 0, −1 certified < 0, 0 structural zero (holds on the whole open stratum), 2 unresolved. |
| `en.i_scanW(...)` | Wlo, Whi (n×4), E, supp | Enclosures of all four frequencies. |
| `en.corner_scan(pat, sign, qlist, nthreads)` | E, cnt (n×7), pos_ok, neg_ok | L3 corner table. Columns: pos, neg, zero_struct, zero_integral, E_c>E, E_c<E, unresolved. pos_ok means N_q ≥ 0 on the whole closed cube. |
| `en.corner_detail(pat, sign, g)` | per-corner Ec, cls, mlo, mhi, Dhi | One co-player. |
| `en.vertex_scan(vlo, vhi, sign, qlist, ...)` | counts per q, unresolved (q, vertex) pairs | Box vertices for L2. |
| `en.binary_pair_exact(pg, qg)` | (A[16] ints, D) | **Exact** binary-pair limit w_k = A_k/D from the interval engine, with each integer isolated. |
| `en.efficiency_class(pat)` | 'C'/'A'/'Z'/'-' | **Exact** efficiency class of the open stratum (L1). |
| `ex.lead_exact(p, q)` | w16, W, E, Ek, D, supp | Exact rationals (gmpy2.mpq, returned as Fraction). |
| `ex.margin`, `ex.freqs`, `ex.selfplay` | | Exact convenience wrappers. |
| `ex.finite_eps(p, q, eps)` | (π, A_k) | Exact finite-ε chain via matrix-tree determinants. Used as a test reference. |
| `ce.certify_point(p, wedge)` | counts, min_margin, selfplay, certified | For a rational point, the **exact** sign of the margin against every q: an interval, a structural zero, or an exact rational on the unresolved ones. Also the exact self-play and efficiency. |
| `ce.certify_box(pat, lo, hi, wedge)` | certified, per-stage counts, volume | **The whole closed box ⊂ friendly rivals.** Each q passes by (a) the L3/L4 corner test, (b) a positive interval over the box, or (c) L2: every vertex ≥ 0 by interval or exact. A negative vertex is a counterexample. |
| `ce.certify_ratio_box(c, h, wedge)` or `lo=, hi=` | certified, U/L, rho_min | **The box c ± h ⊂ friendly rivals** by the ratio route. Each q passes if w_unf is structurally 0, or if w_fav(c)/w_unf(c) ≥ U/L (interval, then exact), or if q is an exact identity. There is no vertex enumeration, so this is the route for d = 15. |
| `ce.certify_open_stratum_all_q(pat, wedge)` | certified, corner_counts | **Every point of the open stratum** is competitive (L3 + L4, positive direction), plus exact efficiency. Refuses when the estimate exceeds 14 min: the cost is 2^d × 65 536 corner evaluations. |
| `ce.exclude_stratum(pat, wedge, candidates=pool)` | excluded, q, centre margin, corner counts | **No point of the open stratum is competitive.** One q has every corner ≤ 0 (L3) and a centre margin < 0, by a strict interval or exact arithmetic (L4). |
| `ce.is_identity(pat, sign, q)` | bool | Whether N_q ≡ 0 on the stratum (every corner value is 0). |

## `witnesses/`

The 25 exact witnesses, written by `make_witnesses.py` and listed in `index.json`. Two of their `source` fields date
from before the paper took its present form. `S_SI2` is the boundary point c of S (stratum `1***1***1**10***`) that
SI Table 2 names and SI section "Certificates" discusses ("A boundary point"), not the point in row S of the table,
c15 = (1,3/4,3/4,1/2, 7/8,3/4,3/4,1/2, 7/8,3/4,3/4,7/8, 1/8,1/2,1/2,1/2), which has no witness file and is certified
with its boxes in `results/S_new15_box.json` (`results/run_Snew15.py`). `S15_box_centre`, marked there as not yet in
the paper, is the centre of the S box of SI Box 2 whose half-width is h = 1/1024 (`results/S_box.json`). The files
are kept as they were generated.

## Timings

Martin's Mac, x86_64 python under Rosetta:

| call | time |
|---|---|
| `f_eval` | 17 µs |
| `i_eval` | 25 µs |
| `f_scan`, 65 536 co-players | 83 ms on 1 thread, 29 ms on 3 |
| `i_scan`, 65 536 co-players | 161 ms on 1 thread, 55 ms on 3 |
| `corner_scan` | 0.78 µs per (q, corner) on 3 threads, 2.3 µs on 1; native arm64 is 1.67 µs on 1 |
| `lead_exact` | 0.32 ms with mpq, 1.07 ms with Fraction |

Certificate runs:

| certificate | time |
|---|---|
| `certify_point` on the W witness | 1.5 s |
| `certify_point` on the N witness | 5.5 s |
| `certify_box` W ± 1/16 (d = 7) | 8 s |
| `certify_box` N ± 1/12 (d = 8) | 15 s |
| `certify_ratio_box` S (d = 15) | 0.08 s |
| `certify_open_stratum_all_q` at d = 5 | 1.7 s |
| `exclude_stratum` | 66 ms cold, 2.8 ms from a pool |
