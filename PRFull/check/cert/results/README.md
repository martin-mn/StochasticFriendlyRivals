# `check/cert/results/`: positive-volume certificates (internal review, round 1, point P2 and items 1–2)

Every statement below is in the ε → 0 limit, taken before any inequality is tested. Each one is quantified over all 65 536 binary memory-two co-players. That this covers every memory-two co-player is the paper's box-vertex lemma, or equally the structural lemma applied to q. Nothing here is sampled: the float scans in the files are sanity checks only.

The library is `../` (`cert.engine`, `cert.exact`, `cert.certify`). The driver is `p2lib.py`.

Run everything from `PRFull/check` as `nice -n 10 python3 cert/results/run_<x>.py`. Each run uses at most 3 threads or processes, and each is under 11 minutes on Martin's Mac.

## Proof routes

Fix q and a stratum. Then N_q = D\*·sign·(w_DC − w_CD) is multilinear in the free coordinates, with D\* > 0 on the open stratum (`../certify.py` gives the lemma, `../tests.py` (c) the checks).

- **Vertex route (L2)**, `p2lib.box_certificate`. For a closed box inside the open stratum, min N_q is attained at a vertex. Every one of the 2^d × 65 536 vertex margins is settled in one of three ways:
  - a strictly positive outward-rounded interval enclosure;
  - a structural zero, meaning w_CD = w_DC = 0 by the ε-exponents alone;
  - an **exact** mpq rational, whenever the enclosure contains 0.

  With `exact_all=False` (used only for E), the zero co-players are instead settled by the L3/L4 corner identity test, N_q ≡ 0.
- **Corner route (L3/L4)**, `certify_open_stratum_all_q`. If N_q(c) ≥ 0 at every corner of the closed stratum cube, then the margin is ≥ 0 on the whole open stratum. The corner values are binary-pair margins, and they are exactly zero by integrality: D\*·margin is an integer.
- **Ratio route**, `certify_ratio_box`. Each A\*_k is a nonnegative sum of monomials in {x_s, 1 − x_s}. So on the box, w_fav/w_unf ≥ ρ(c)·L/U. This is the route for d = 15.
- **Hybrid**, `p2lib.hybrid_box_certificate`. Each co-player is settled by the ratio route where it suffices, and otherwise by the vertex route over all 2^15 vertices.
- **Efficiency is exact (L1).** The self-play support is fixed by the exponents, which depend only on the pattern, so the class is constant on the open stratum. The exact self-play is also computed at the centre and at every vertex.
- **Lemma check.** In every vertex-route file, N_q at random rational interior points equals exactly its multilinear interpolation from the vertices.

## Results

| file | region (closed box unless stated) | d | struct. zero | positive by interval | resolved exactly | smallest positive vertex margin (exact) | time |
|---|---|---|---|---|---|---|---|
| `W_box` | W witness ± 1/16 | 7 | 33 260 | 27 540 | 4 736 (606 208 pairs, all exactly 0) | 736/691793 ≈ 1.06e-3 | 97 s |
| `W_box_large` | p1, p3 ∈ [1/256, 255/256], p7 ∈ [1/256, 95/128], p9, p11, p14 ∈ [1/256, 5/16], p12 ∈ [1/256, 73/128]; volume 0.0121, about 2.5·10⁴ × the ±1/16 box | 7 | 33 260 | 27 540 | 4 736 (all 0) | 11776/30087256599 ≈ 3.9e-7 | 87 s + 341 s search |
| `N_box` | N witness ± 1/12 | 8 | 13 436 | 34 512 | 17 588 (4 502 528 pairs, all exactly 0) | 43/10760 ≈ 4.0e-3 | 515 s (first run 600 s) |
| `E_box` | E witness ± 1/12 (mirror of N) | 8 | 13 436 | 34 512 | 0 exact; 17 588 identities by the corner test | 43/10760 | 45 s |
| `cube5_stratum` | the **whole open** stratum `1*1*0*0101*1*010` | 5 | 17 088 | – | all 65 536 have N_q ≡ 0 (48 448 non-structural identities) | margin ≡ 0 | 3.5 s (+34 s faces) |
| `S_box` | c = (32,27,15,25,29,28,9,6,21,14,27,15,6,12,30,22)/32 ± 1/1024 | 15 | 25 159 | 36 667 by ratio + 3 710 with w_DC ≡ 0 | 0 | ρ_min = 419698/307347 vs U/L = 1.172756 | 9 s |
| `S_new11_box` | c′ ± 3/16, in the old S witness's stratum | 11 | 24 122 | 41 414 | 0 | 85/13974532 ≈ 6.1e-6 | 226 s |
| `S_new15_box` | c15 ± 1/128 (hybrid); c15 ± 1/256 (ratio only) | 15 | 25 159 | 24 070 by ratio + 3 710 with w_DC ≡ 0, and 12 597 by all 32 768 vertices | 0 | ρ_min = 9/5 | 362 s |

The minimum over all (vertex, co-player) pairs is exactly 0 in every file. It is attained by the structural zeros and by the identities, the co-players whose margin is identically 0 on the stratum. The centre of each W, N and E box is the SI Table 2 witness, so those three witnesses are **interior** points of the friendly-rival set.

## The 5-cube

`cube5_stratum` decides all 3⁵ strata of the closed cube:
- 108 strata, those with p1 ∈ {\*, 1}, p5 ∈ {\*, 0} and p3, p10, p12 arbitrary, are certified friendly-rival strata;
- the 81 strata with p1 = 0 are inefficient, which is exact;
- the 54 strata with p5 = 1 and p1 > 0 are each excluded by co-player 16519.

So the friendly rivals in the closed cube are **exactly** the points with p1 > 0 and p5 < 1. This proves the sentence "appears to hold friendly rivals only" of an earlier draft of the paper.

## The S witness of SI Table 2 is a boundary point (`S_witness_boundary`)

The witness is (1,½,½,½, 1,½,½,½, 1,½,½,1, 0,½,½,½).
- 6 728 co-players tie with it exactly without being structural zeros, and **none** of these ties is an identity.
- N_q is multilinear, hence harmonic. By the strong minimum principle, each of these q therefore has points with margin < 0 in every neighbourhood of the witness.
- Explicit exterior points: p2 − δ against q = 1142 gives an exact margin of −3/692 at δ = 1/100 and −3.9e-13 at δ = 2⁻⁴⁰. The six points of the internal review notes (p1, p2, p5, p6, p9 or p10 lowered by 1/100) are beaten by 6728, 6728, 1088, 1288, 1928 and 1928 co-players respectively, with exact minima between −1/235 and −1/451.

**Proposed replacement row for SI Table 2.** Both candidates are certified; `S_new15` is recommended because it lies in the maximal stratum.

- `S_new15`: c15 = (1,¾,¾,½, ⅞,¾,¾,½, ⅞,¾,¾,⅞, ⅛,½,½,½), stratum `1***************`, dimension 15. Self-play (1,0,0,0). Minimum margin 0, attained by 25 159 co-players (all structural). Certified box ± 1/128.
- `S_new11`: c′ = (1,¾,¾,½, 1,¾,¾,½, 1,¾,¾,1, 0,½,½,½), stratum `1***1***1**10***`, dimension 11. Self-play (1,0,0,0). Minimum margin 0, attained by 24 122 co-players (all structural). Certified box ± 3/16. This is the minimal change: six entries go from ½ to ¾.

Both are defensible by the SI's float linear programme (M2f `m2lib.occupation_lp`, value 0). That check is in double precision and is not exact.

## Caveats

- The "1/37 gap" is false (see `../README.md`), and it is not used anywhere here.
- `W_box_large` was found by a greedy float-and-interval screen and then certified from scratch. The three faces p9, p11, p14 = 5/16 cannot be moved outward, even by 1/256.
- `N_box` was rerun after a memory-only refactor of `p2lib` (the enclosures are now streamed). `W_box` and `E_box` were rerun too, and they are identical apart from timings.
