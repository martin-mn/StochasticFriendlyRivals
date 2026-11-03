# `check/ubound/`: upper bounds on the dimension of the friendly-rival set (internal review, round 1, item 1b)

Everything here is in the ε → 0 limit, taken before any inequality is tested. Nothing is sampled.

## What is proved

For each locus below, every open support stratum of the stated dimensions that satisfies the locus's efficiency requirement is **excluded**. Excluded means one binary co-player q gives sign·(w_DC − w_CD) < 0 at *every* point of the open stratum.

| locus | sign of T − S | efficiency (exact class) | excluded strata | count |
|---|---|---|---|---|
| wedge W | + | C (w_CC = 1) | all of dimension ≥ 8 | 795 531 |
| wedge N | + | A (w_CD + w_DC = 1) | all of dimension ≥ 9 | 107 694 |
| wedge E | − | A | all of dimension ≥ 9 | 107 694 |
| ray u+v=1, u<1 | + | C, A or Z (w_DD = 0) | all of dimension ≥ 9 | 444 722 |
| W, face p₆ = p₈ = 0 | + | C | all of dimension ≥ 6 with p₆, p₈ pinned at 0 | 178 760 |

Twelve co-players suffice for all of them (`qstar.txt`). Combined with the certified boxes in `../cert/results/`:
- W has dimension exactly 7 (W_box).
- N and E have dimension exactly 8 (N_box, E_box).
- The T > S ray has dimension exactly 8: the N box is class A with sign +, so it lies in the ray's set.
- The face p₆ = p₈ = 0 has dimension exactly 5: the whole open 5-cube `1*1*0*0101*1*010` is certified (cube5_stratum).

## Method (`ub.c`)

**The exclusion criterion.** N_q = D*·margin is multilinear on the closed cube of a stratum. For x in the open cube it is a combination of the 2^d corner values N_q(c) with strictly positive weights.

By lemma L3:
- N_q(c) = D*_c·m(c, q), with D*_c ≥ 1, when E_c(q) = E_σ(q);
- N_q(c) = 0 when E_c(q) > E_σ(q);
- E_σ = min_c E_c.

Therefore q excludes σ exactly when every corner at the minimum exponent has m ≤ 0 and at least one has m < 0. This needs no interior evaluation and no gap constant (no "1/37").

**The computation.**
1. For each q, compute a row: E_c and the sign of m(c, q) for all 65 536 binary c. These come from the outward-rounded interval engine, with integrality deciding exact zeros.
2. Run a ternary DP over all 3¹⁶ strata, one child per lowest free coordinate. Each stratum gets (min E, OR of the flags at min E).
3. Cost: 0.18 s per co-player on one core.
4. The co-players were chosen adaptively. For strata not yet excluded, the float margin at the stratum's centre ranks the candidates. The ranking only proposes candidates; the proof is the DP. 182 co-players were tried and a greedy set cover kept 12.

**Checks, all passed.**
- `ub verify`: every one of the 1 189 899 (stratum, sign, q) triples was re-derived by a separate route, direct enumeration of all 2^d corners (545 639 040 corner evaluations). For every triple:
  - no E_c < E_σ, and E_σ (taken at the centre) = min E_c;
  - the corner criterion holds;
  - the centre margin enclosure is strictly negative.
- `exact_rows.py`: the 12 rows (786 432 binary pairs) were recomputed in exact rational arithmetic (gmpy2). They are byte-identical to the interval rows, so the certificate is free of floating point.
- `exact_check.py`:
  - exact margins < 0 at the centre and at a random k/97 point of 2000 target strata;
  - 20 000 binary pairs exact;
  - 4000 efficiency classes exact at random k/89 points.
- `census_check.py`: the 1 524 791 census strata (M2f `st_eff.txt`) have the same exact class, and none of the census friendly-rival strata is excluded.

## By-product: the census is complete at the stratum level

With 1681 co-players (`results/pool_c2_q.txt`, found adaptively, all dimensions; `./run_all.sh all`), the efficient strata **not** excluded are:
- **W: 364.** These are exactly the corrected census set (st_eff's 362 plus g = 5648342 and 5648369). Split by (p₆, p₈): 126 at (0,0) and 238 at (1,1), none elsewhere. The 50 off-diagonal strata left after the first pool were each excluded by an exhaustive single-co-player search (`offdiag.py`, `runs/offdiag.json`).
- **N: 5022, and E: 5022.** The two dimension histograms are identical. The N set is st_eff's N set ∪ mirror(E).
- **Ray T>S: 7726.** This is st_eff's 7721 plus 5 strata, matching the corrected census count.

So the census counts, which were lower bounds, are also upper bounds.

Consequence for W:
- The friendly-rival set lies in the two disjoint closed faces {p₆ = p₈ = 0} and {p₆ = p₈ = 1}.
- The two families are therefore separated.

All 14 241 667 all-dimension exclusions (every excluded efficient stratum of W, N, E and the ray, with its co-player) passed the independent `ub verify` re-check (`results/verify_alldims_c2.out`, 2 x 390 s).

## Timings (Mac, 3 threads, nice 10, machine shared)

| step | time |
|---|---|
| classes of all 3^16 strata | 15 s |
| one co-player (row + DP + marking) | 0.18 s on 1 core |
| 12-co-player certificate | 0.7 s |
| adaptive search to 182 co-players | 22 s |
| `ub verify` of the 1 189 899 triples | 325-367 s |
| exact rows of the 12 co-players | 77 s |
| 1681-co-player by-product | 102 s, verify 2 x 390 s |

## Files

| file | contents |
|---|---|
| `ub.c`, `certeng.c` | modes `classes`, `cover`, `verify`, `sets`, `row`, `rowcmp`, `check`. `certeng.c` is a verbatim copy of `../cert/certeng.c` |
| `status.py`, `greedy.py`, `setcover.py`, `mktriples.py`, `exact_rows.py`, `exact_check.py`, `census_check.py`, `analysis.py`, `extra.py`, `offdiag.py`, `export.py`, `counts.py` | drivers and checks |
| `qstar.txt` | the 12 certificate co-players |
| `results/ub_certificate.tsv.gz` | one line per excluded (stratum, sign) with its co-player |
| `results/certificate_summary.txt` | the co-players and the count tables |
| `results/*.out` | logs |
| `run_all.sh` | reproduces the certificate from scratch, in about 10 min on 3 threads (`all`: plus the by-product, about 15 min more) |
| `results/pool_c2_q.txt` | the 1681 co-players of the by-product |

The bulky intermediate arrays (per-stratum co-players, triples, verification outputs, 1.3 GB) are not kept here: `run_all.sh` regenerates them. `cls.bin` (the classes, 43 MB) and `dim.npy` are regenerated by steps 1 and 2.

Stratum encoding in the binaries: index = Σ_j dig_j 3^j, where dig 0 means pinned 0, 1 means pinned 1 and 2 means free. This differs from the census code; see `census_check.py:conv`.
