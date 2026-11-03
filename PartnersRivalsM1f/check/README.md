# `check/` — exact ε → 0 checks over the memory-one continuum

These programs check the memory-one results of the paper (see the paper's SI section on memory one)
over the continuous strategy space p = (p₁,p₂,p₃,p₄) ∈ [0,1]⁴. Where possible the limit ε → 0 is
taken exactly, with no threshold. Where floats are used, the two possible outcomes are separated by
many orders of magnitude, and the separation is reported.

The checks come in four groups, called items in the scripts: rivals (item 1), efficiency and
friendly rivals (item 2), partners (item 3) and the dimension of the limit-Nash set (item 4). The
same checks over the sixteen binary strategies are in `../../PartnersRivalsM1/check/`. The Cannon
kit used for the large runs is `../cannon/`.

## Contents

| file | what it is |
|---|---|
| `m1lib.py` | shared harness, imported by every script |
| `sweep_rivals.py`, `sweep_rivals.out` | item 1, float sweep of the rival sets |
| `exact_rivals.py` | item 1, exact Fractions with the leading ε-exponent (no saved output) |
| `strata_rivals.py` | item 1, all 81 strata by exact leading-order arithmetic (no saved output) |
| `symbolic_rivals.py`, `symbolic_rivals.out` | item 1, sympy with ε and pᵢ symbolic, Bernstein sign certificates |
| `sweep_efficiency.py`, `sweep_efficiency.out` | item 2, float sweep of the efficient and friendly-rival sets |
| `exact_efficiency.py`, `exact_efficiency.out` | item 2, exact Fractions |
| `sweep_partners.py`, `sweep_partners.out` | item 3, float check of the partner sets against the 16-deviator test |
| `exact_partners.py`, `exact_partners.out` | item 3, exact Fractions |
| `nash_field.py`, `nash_field.out` | item 4, dimension of the limit-Nash set at each game |
| `working/` | development scratch: derivations, intermediate scans and one-off probes; not receipts |

Each script prints its report to standard output; a `.out` file is the saved report of a run. The
closed forms each script tests are written out in its docstring.

## Conventions

All conventions are fixed in `m1lib.py` and are shared with `../../PartnersRivalsM1`.

- Payoffs are normalised to (R,S,T,P) = (1,u,1+v,0).
- A strategy p = (p₁,p₂,p₃,p₄) gives the probability of cooperating after CC, CD, DC and DD. The
  focal player's own action is written first. States are numbered CC=0, CD=1, DC=2, DD=3.
- Error model: each probability p is replaced by ε + (1−2ε)p. The strategy itself stays a point of
  the closed cube, and ε is only a floor on the error rate.
- The sixteen binary strategies are the vertices of the cube. Their genome strings list p₁…p₄ in
  order, 1 for C (so TFT is `1010`, WSLS is `1001`, ALLD is `0000`).
- The pair chain has four states; w = (w_CC, w_CD, w_DC, w_DD) is its stationary distribution, and W
  is the limit of w as ε → 0. The payoff difference of the two players is π(p,q) − π(q,p) = (T−S)(w_DC − w_CD).
- Co-players and deviators are the 16 binary strategies. Why that suffices is in the paper's SI
  section on co-players and deviators of any memory; segments 2 and 3 of the Cannon run below test
  it numerically against all memory-two co-players.

`m1lib.py` provides the pair chain in float and in exact rationals, payoffs, `rival_margin`,
`nash_gap`, `efficiency_deficit`, the sixteen vertices (`DET`, `DET_EXACT`), and the predicates
`in_riv_plus`, `in_riv_minus`, `in_equaliser`, `in_eff_cc`, `in_eff_alt`, `in_eff_line` (with the
helper `_alt_corner`). Its exact route is checked in two ways: Akin's identity has residual exactly 0
in Fractions, and TFT in self-play gives exactly (¼,¼,¼,¼).

## Requirements and running

Python 3 with numpy. `nash_field.py` also needs scipy, and `symbolic_rivals.py` needs sympy. Several
scripts use a pool of 12 worker processes. Run every script from this directory.
`exact_efficiency.py` imports `exact_rivals.py`, and `sweep_efficiency.py` imports
`sweep_rivals.py`. `sweep_partners.py` imports `partners_lib` from `working/`, so run it with
`PYTHONPATH=working`.

---

# Item 1 — the rival sets

```bash
python3 sweep_rivals.py      # float: grid, random interior points, dense sample of facets and edges
python3 exact_rivals.py      # exact Fractions: 456 976 pairs, leading ε-exponent computed exactly
python3 strata_rivals.py     # all 81 strata of the cube, with explicit witnesses
python3 symbolic_rivals.py   # sympy, ε a symbol: Bernstein sign certificates, cell by cell
```

## Results

| route | scope | result |
|---|---|---|
| `sweep_rivals.py` | ε = 10⁻¹⁰ and 10⁻¹². A 21⁴ = 194 481 grid; 200 000 uniform interior points; 69 152 points on all 8 facets, all 24 two-dimensional faces, the edges and thin slabs near them | 0 disagreements. Grid: 4830 rivals predicted and 4830 computed, for each sign of T−S. Facet sample: 11 626 and 11 626; it lies on the faces p₄ = 0 and p₁ = 1. The lowest margin of a predicted rival is −2.0×10⁻⁹ and the highest margin of a non-rival −1.7×10⁻², so any tolerance from 10⁻⁸ to 10⁻³ gives the same verdicts. Report: `sweep_rivals.out` |
| `exact_rivals.py` | 28 561 rational strategies × 16 co-players = 456 976 pairs; exact Fractions; leading ε-exponent read off exactly, not fitted | 0 disagreements. Akin's residual is the zero polynomial in ε on every pair. Riv⁺ ∩ Riv⁻ is exactly the 12 grid points (1,s,1−s,0) with s < 1; s = 1 is excluded. Leading exponents 0, 1, 2, 3 occur 396 224, 58 170, 848 and 12 times; 1 722 margins are identically zero; no exponent is negative. The co-players {ALLD, `0001`} alone reproduce Riv⁺, and {ALLC, `0111`} alone reproduce Riv⁻; ALLD alone gives 2 028 points, 858 too many |
| `strata_rivals.py` | all 81 strata (each pᵢ at 0, interior or 1) × 16 co-players | the closed form holds on every stratum; 132 explicit witnesses |
| `symbolic_rivals.py` | sympy; ε and all four pᵢ symbolic. The cube is cut into cells (each pᵢ fixed at 0 or 1, or free); the cell's parametrisation is substituted first and the limit taken afterwards; each sign is certified by non-negative Bernstein coefficients on the unit box | All certified, in 49 s. Akin's residual is exactly 0 against all 16 co-players with every pᵢ free. On all 54 cells with p₄ > 0, ALLD gives the exact limit −p₄/(1−p₂+p₄) < 0. The rival inequality is certified on every cell of {p₄=0, p₂+p₃≤1, p₂<1}, and a refuting co-player is found on every cell piece where p₂+p₃ > 1. The mirror identity holds for all 16. Report: `symbolic_rivals.out` |
| Cannon array **45636422** | 512 tasks, all COMPLETED | 73 003 425 522 strategies, 0 disagreements, no `DISAGREE` line; see below |

## The Cannon run

The kit is `../cannon/` (source `m1f_pairs.c`, `compile.sh`, `task.sh`, `riv1.sbatch`,
`riv1_probe.sbatch`; see its README). The full array **45636422** (`../cannon/JOBID`) ran on
`sapphire`: 512 tasks, 72 core-hours. It was sized from the probe **45630088**
(`../cannon/JOBID_PROBE`, partition `test`), which measured 44 ns per four-state limit on the grid,
78 ns on the random sample and 1.62 µs per sixteen-state limit.

The array index is `segment × 128 + shard`. The sizes used are those exported in `riv1.sbatch`
(`N1=500 N2=10000000000 N3=30 N4=1000000`). The defaults in `task.sh` and the segment table in
`../cannon/README.md` are smaller.

| segment | mode | sweep | co-players |
|---|---|---|---|
| 0 | `m1` | exhaustive grid of spacing 1/500, 501⁴ = 6.30×10¹⁰ strategies | the 16 memory-one vertices |
| 1 | `rand` | 10¹⁰ random points drawn from all 81 strata, off the grid | the 16 memory-one vertices |
| 2 | `m2` | grid of spacing 1/30, 923 521 strategies | all 65 536 memory-two strategies |
| 3 | `rand2` | 10⁶ random points from all 81 strata | all 65 536 memory-two strategies |

Segments 2 and 3 test the restriction to memory-one co-players: no memory-two co-player does better
against a memory-one resident than the best memory-one co-player does.

**Receipts.** `sacct -j 45636422 -X` reports 512 COMPLETED and nothing else. The output files are
`../cannon/out/riv1_45636422_<task>.out` and `.err`. Each `.out` ends in one `SUMMARY` line. Summed
over the 512 lines: tested 73 003 425 522, `bad_plus` 0, `bad_minus` 0, Riv⁺ 1 063 129 628, Riv⁻
1 063 119 060. No `.out` file contains a `DISAGREE` line. Each `.err` file holds only the
`/usr/bin/time -v` report, with exit status 0. To recompute the totals:

```bash
cd ../cannon/out
grep -h '^SUMMARY' riv1_45636422_*.out | awk '{for(i=2;i<=NF;i++){split($i,a,"=");s[a[1]]+=a[2]};n++}
  END{printf "%d tasks tested=%.0f bad_plus=%d bad_minus=%d riv+=%.0f riv-=%.0f\n",
      n,s["tested"],s["bad_plus"],s["bad_minus"],s["rivals_plus"],s["rivals_minus"]}'
grep -l DISAGREE riv1_45636422_*.out | wc -l       # 0
```

## Tolerances

The exact-rational routes use no tolerance, and Akin's residual is checked as a polynomial identity.
The float and Cannon routes call a limit zero when its absolute value is below 10⁻¹². The Cannon
engine computes limits by a leading-order GTH reduction, which has no subtractions, so each computed
limit is either exactly zero or far from it. On the float grid the measured gap is seven orders of
magnitude.

One place needs care: the edge E = {(p₁,1,0,0)}. At ε = 0 the chain on E has two recurrent classes,
so the limit depends on the O(ε) transitions between them and has to be computed at finite ε. The
ALLD margin on E is −1/(3−p₁), which lies in [−½,−⅓].

---

# Item 2 — efficiency and friendly rivals

```bash
python3 sweep_efficiency.py   # float sweep of the three efficiency forms and the friendly-rival sets
python3 exact_efficiency.py   # exact Fractions, no threshold
```

The efficiency predicates are `in_eff_cc`, `in_eff_alt` and `in_eff_line` in `m1lib.py` (helper
`_alt_corner`); the sweeps compare against them.

| route | scope | result |
|---|---|---|
| exact Fractions | grids with 5, 9 and 11 values per coordinate (5⁴, 9⁴, 11⁴); the order in ε tested by solving at both ε = 10⁻³ and ε = 10⁻⁶ | 0 mismatches. Strategies with W_CC = 1, with alternation, and with W_DD = 0: 99, 30, 137 on 5⁴; 647, 126, 789 on 9⁴; 1209, 198, 1427 on 11⁴. w_CD = w_DC exactly in all 1250 solves |
| `sweep_efficiency.py`, `exact_efficiency.py` | as in their docstrings | reports in `sweep_efficiency.out` and `exact_efficiency.out` |
| Cannon array **45658567** | 256 tasks, all COMPLETED. A grid of spacing 1/500 (6.30×10¹⁰) and 10¹¹ points drawn from all 81 strata. The closed-form predicates are not used: each point's self-play limit is classified, both rival margins are tested, and friendly rivals are counted per wedge W, N, E, S of the (u,v) plane | 163 001 502 001 strategies. FR_W = 0, FR_N = 0, FR_E = 0, FR_S = 8 543 035 103 |

The same array also counts efficient strategies without using the efficiency closed forms. The
closed forms predict effALT = 4 800 499 998 and effCC = 15 965 500 499; the array measured
4 800 499 428 and 15 965 745 936, relative differences 1.2×10⁻⁷ and 1.5×10⁻⁵, within sampling noise.

**Cannon details.** Segments 4 (`fr`, grid) and 5 (`frand`, random) of `../cannon/task.sh`, run by
`../cannon/fr1.sbatch` (`N5=500 N6=100000000000`). Array **45658567** is in `../cannon/JOBID_FR`. It
was sized from the probe **45656466** (`../cannon/JOBID_FR_PROBE`, `fr1_probe.sbatch`), which
measured 0.806 µs and 1.51 µs per point. Each `../cannon/out/fr1_45658567_<task>.out` ends in one
`FR-SUMMARY` line; the totals above are the sums of those 256 lines, and every `.err` file holds only
the `/usr/bin/time -v` report with exit status 0.

## Floating-point conditioning

The float pair chain of `m1lib.py` (`stationary`, a least-squares solve) loses accuracy for
ε ≤ 10⁻¹¹. At ε = 10⁻¹¹ it gave errors of about 5×10⁻⁶ and missed 13 alternators (out of 30) on
the 5⁴ grid, all at a single corner. The p₂ ↔ p₃ symmetry of self-play, which the bad answer broke, exposed
the problem. With that solver, use ε = 10⁻⁷ or 10⁻⁹; for decisive checks use exact Fractions or
`../cannon/m1f_pairs.c`. (`sweep_rivals.py` and `sweep_efficiency.py` do not use the least-squares
solve: they use a subtraction-free Markov-chain-tree computation, which keeps full relative accuracy
at ε = 10⁻¹².)

---

# Item 3 — Nash equilibria and partners

```bash
PYTHONPATH=working python3 sweep_partners.py   # float: closed forms vs the raw 16-deviator test, over a grid of games
python3 exact_partners.py                      # exact Fractions, no threshold
```

| route | scope | result |
|---|---|---|
| `sweep_partners.py` | 109 Eff_CC and 136 Eff_alt candidates × 3819 games (the grid of spacing 1/2 on [−8,8]², plus points on and ±0.01 off every boundary of the closed forms), about 936 000 resident–game pairs | 0 disagreements with the closed forms below, on and above the line u+v = 1; 0 differences in classification between ε = 10⁻⁹ and 10⁻⁷ |
| exact Fractions | 260 tests comparing the six-inequality and two-inequality forms; 939 over all 17 remaining boundary strata of Eff_CC; 503 random interior points; 302 on the (1,0) sheet and 338 on the (0,1) sheet, with T = S over-sampled; 736 against Akin's ratio form where it applies | 0 mismatches |
| Nash existence | the 1089 games of the grid of spacing 1/2; witness p = (c,c,c,c), with c the cooperation probability of the symmetric equilibrium of the one-shot game | 0 failures; worst gap −8.0×10⁻⁹ at ε = 10⁻⁹ |
| graded field | closed forms against deterministic brute-force quadrature | largest difference 2.8×10⁻⁵ below the line, 1.5×10⁻⁴ above |

Reports: `sweep_partners.out`, `exact_partners.out`. `sweep_partners.py` decides limit-Nash
membership from two values of ε: a point counts as limit-Nash if gap(10⁻⁹) ≥ −10⁻⁵ or
|gap(10⁻⁹)| < 0.1·|gap(10⁻⁷)|.

There is no Cannon run for item 3, and no `JOBID` file for it.

---

# Item 4 — the dimension field of the limit-Nash set

```bash
python3 nash_field.py       # dimension field, carriers, existence; ε = 10⁻⁷ and 10⁻⁹
```

`nash_field.py` computes, for each game (u,v), the dimension of the limit-Nash set N(u,v) ⊂ [0,1]⁴
and tests which strata carry its three-dimensional pieces. Q denotes the closed quadrant
{u ≥ 1, v ≥ 0}.

- Games: 322, namely the 289 games of the unit grid on [−8,8]², plus 33 probes: on both boundary
  rays, at distance 10⁻³ on either side of them, and at the corner (1,0).
- Carriers: each claimed carrier gets a Chebyshev ball, found by linear programming inside its
  stratum, and the ball is then tested against the raw 16-deviator test. Non-carriers are probed by
  uniform sampling and by Nelder–Mead maximisation of the Nash gap from many starting points.
- When a certificate ball is too thin for the float chain, the test is repeated in exact Fractions.
  This was needed at 2 of the 322 games, (0.999, 3) and (0.999, 8).

Result (`nash_field.out`, run time 114 s): the measured dimension equals the claimed one at all 322
games, and every failure and mismatch count in sections 1–8 is 0.
`../Figures/fig_nash.py` repeats the test of both facet criteria by exact linear programming on a
65 × 65 grid of games, with 0 mismatches, and also checks the saturation of the live alternating
sheet on all of Q.

**The limit test.** For a limit-Nash strategy the gap at finite ε is of order ε, but its constant
becomes large in boundary layers (for example, p₄ of size O(ε) on {p₁=1}). A fixed tolerance at a
single ε therefore misclassifies thin regions. `nash_field.py` evaluates the gap at ε = 10⁻⁷ and
10⁻⁹ and calls a point limit-Nash if gap(10⁻⁹) > −10⁻⁶ or gap(10⁻⁹) > 0.2·gap(10⁻⁷). Between these
two ε the gap of a limit-Nash point shrinks by a factor of 100, while a genuine failure does not
shrink. Section 9 of `nash_field.out` compares the two rules on 1200 strategy–game pairs: a fixed
tolerance of 10⁻⁴ at each ε separately disagrees once, at a boundary-layer point that the two-ε rule
calls Nash at both; 188 points are classified Nash by the ratio branch, with ratios from 5.2 to 100.

There is no Cannon run for item 4.
