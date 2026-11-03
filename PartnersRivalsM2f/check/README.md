# `check/` — local checks and receipts for the memory-two continuum

This directory holds the validation of the memory-two engine, the exact certificates of single points, small local
analyses, and the local copies of the Cannon output with the scripts that collected and summarised it. Each local
program runs on a 12-core machine in at most about 15 minutes (one exception is noted at `boxcover.py`). The Slurm
kits are in `../cannon/`.

Many programs name their inputs by absolute paths on the author's machine (`/Users/martin/Documents/CL1/...`).
`PRFull/release/relocate.py` points them at the unpacked release.

Not in the release: `working/` (development scratch), apart from the programs and output of five checks that
`working/README.md` lists; `derivations/` (internal notes); and `__pycache__/`. No program here reads them.

## Conventions

* Payoffs (R, S, T, P) = (1, u, 1 + v, 0). T = S on the line u − v = 1; the switch line is u + v = 1.
* Outcome code o = 2·(focal player defects) + (co-player defects): CC = 0, CD = 1, DC = 2, DD = 3, the focal
  player's own action first.
* Memory-two state s = 4·(last outcome) + (outcome before it), s = 0, …, 15. A strategy is p ∈ [0,1]¹⁶, with p[s] the
  probability of cooperating after state s.
* Binary genome: a 16-bit integer whose bit s is 1 for C in state s (ALLC = 65535, ALLD = 0).
* Errors: both players play ε + (1 − 2ε)p; every result is the limit ε → 0, taken at leading order in a GTH
  (Grassmann–Taksar–Heyman) elimination.
* W = (W_CC, W_CD, W_DC, W_DD) are the limit frequencies of the four outcomes. Efficient: W_CC = 1 below the switch
  line, W_CD = W_DC = ½ above it. Rival (competitive) for the sign of T − S: the minimum of
  sign(T − S)·(W_DC − W_CD) over the co-players is ≥ 0 (tolerance 1e-9 in floating point). The scans run over the
  65536 binary co-players; the paper's SI section "Co-players and deviators of any memory" explains why that suffices.
* Wedges: W (u + v < 1, T > S), S (u + v < 1, T < S), N (u + v > 1, T > S), E (u + v > 1, T < S), and the two rays
  of the switch line (T > S and T < S).
* Stratum: the 16 coordinates each marked 0 (p_j = 0), interior (0 < p_j < 1) or 1 (p_j = 1). Its index is
  g = Σ d_j·3^j with d_j = 0, 1, 2 in that order; as a pattern it is printed with `0`, `*`, `1`, position j = state j.
* Census class of a stratum: `C` (W_CC = 1), `A` (perfect alternation), `Z` (W_DD = 0 only), `-` (not efficient).

## The engine

| file | what |
|---|---|
| `../cannon/m2f.c` | the ε → 0 engine in C: leading-order GTH for a pair p, q ∈ [0,1]¹⁶ |
| `m2f.py` | ctypes driver for `../cannon/libm2f.dylib`: `freqs(p, q)`, `freqs16`, `selfplay(p)`, `rival_margin(p, sign)` (minimum over the 65536 binary co-players and the co-player that attains it), `batch`, `batch_wit`, `g2p`, `p2g` |
| `exact.py` | the same algorithm with `Fraction` coefficients and integer exponents: exact rational limits, no floating point |
| `m2lib.py` | the conventions above as code, and `occupation_lp(p, sense)`: a linear programme over the long-run state–action frequencies of the MDP that the co-player controls, giving the exact minimum of sense·(W_DC − W_CD) over every co-player, with no bound on its memory and from any start; `is_rival(p, sense)` |

The release ships no compiled code. Build the library before using `m2f.py`:

```bash
cd ../cannon
cc -O3 -arch arm64 -arch x86_64 -fPIC -shared -o libm2f.dylib m2f.c m2f_py.c -lm
```

The universal build serves an arm64 or an x86_64 Python. On Linux, drop the `-arch` flags and keep the output name:
`m2f.py` loads that file by its path.

### Validation of the engine

* Self-play limits of all 65536 binary genomes: byte-identical to those of an independent program for binary
  strategies (not in this release).
* Rival sets of the binary genomes, one per sign of T − S, 2640 genomes each: equal as sets, not only in size, to
  `rivP` and `rivM` of `PartnersRivals/Figures/m2_masks.npz`. `binary_margin.npz` holds the T > S margins.
* `exact.py` against the C engine: agreement on 200 random binary pairs and on continuum pairs.
* Memory-one strategies embedded as p[s] = π[s >> 2] (they look at the last outcome only) give the classes of the
  memory-one engine of `../../PartnersRivalsM1f`. In 20 000 stratified draws the friendly-rival counts by wedge are
  W 0, S 2590, N 0, E 0. TFT has self-play limit (¼, ¼, ¼, ¼) and is a rival for both signs.

### Tolerance

Over all 65536 binary strategies against all 65536 binary co-players, a margin W_DC − W_CD that is not zero is at
least 1/37 ≈ 0.027 in absolute value, and the margin of every rival is exactly 0. The tolerance 1e-9 therefore never
decides a close case there. The continuum points below are certified in exact rationals.

## Exact certificates

Each script uses `exact.py`, scans all 65536 binary co-players with 12 processes and takes well under a minute.

| script | point p = (p₀, …, p₁₅) | stratum, wedge | exact output |
|---|---|---|---|
| `certify7.py` | (1, ½, 1, ½, 0, 0, 1, ½, 1, ¼, 1, ¼, ¼, 0, ¼, 0) | `1*1*001*1*1**0*0`, W | self-play (1, 0, 0, 0); minimum of W_DC − W_CD is 0, at 37 996 co-players |
| `certify8N.py` | (⅔, 1, 0, ⅓, ⅓, ⅓, 0, 0, ⅔, 1, ⅔, ⅓, ⅓, 0, 1, 0) | `*10***00*1***010`, N | self-play (0, ½, ½, 0); minimum of W_DC − W_CD is 0, at 31 024 co-players |
| `certify11S.py` | (1, ½, ½, ½, 1, ½, ½, ½, 1, ½, ½, 1, 0, ½, ½, ½) | `1***1***1**10***`, S | self-play (1, 0, 0, 0); minimum of W_CD − W_DC is 0, at 30 850 co-players |
| `certify_segment.py` | p(t) = (0, 0, 0, t, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0) | N | self-play (0, ½, ½, 0) at t = 1/10, 1/4, 1/2, 3/4, 9/10; (⅙, ⅓, ⅓, ⅙) at t = 0; (²⁄₉, ⅓, ⅓, ¹⁄₉) at t = 1; minimum of W_DC − W_CD at t = ½ is 0 |

The W and N points are those of the paper's SI Table 2; the S point is the simpler one named in that table's legend.
See also the paper's SI section "Certificates". The docstring of `certify8N.py` repeats that of `certify7.py`; the code
certifies the N point above. `PRFull/check/certify8N.py` is a copy with a corrected docstring.

## Other local programs

| file | what it computes |
|---|---|
| `validate_binary.py` | runs `occupation_lp` on all 65536 binary genomes for both signs (12 processes), compares with `PartnersRivals/Figures/m2_masks.npz` and writes `binary_lp.npz`. Output: 2144 genomes per sign, equal as sets to `defP`/`defM` of that file and not to its 2640 `rivP`/`rivM`, because the LP admits co-players of any memory |
| `witness.py` | greedy set cover over binary co-players; writes `witnesses.npy`, 19 co-players such that each of the 62 896 binary genomes that is not a rival for T > S has a negative margin against at least one. Used only as a screen |
| `boxvertex.py` | numerical test that the worst memory-two co-player of a fixed p is binary: 24 rivals (the 8 binary friendly rivals of W, the point of `certify7.py`, 15 random rivals in its stratum), each against 40 Nelder–Mead searches and 4000 random q ∈ [0,1]¹⁶. Output: no q falls below the binary minimum by more than 1.7 × 10⁻¹⁶ |
| `explore_W.py` | moves one coordinate at a time of each of the 8 binary friendly rivals of W (`111*000101*1*010`) over a grid in [0,1] and records efficiency and rivalry |
| `boxcover.py` | `python3 boxcover.py <tag> '<glob of stratum files>'`: minimum cover of a set of strata by boxes, products of subsets of {0, interior, 1} per coordinate (prime boxes, then an integer programme with `scipy.optimize.milp`). Do not run it on wedge N or S locally: the prime step grows exponentially and did not finish in 15 minutes on the 4927 strata of `strata/N_strata.txt` |
| `rebuild_refs.py` | orders a LaTeX reference list by first citation; its input, the project's manuscript draft, is not in this release |

| data file | contents |
|---|---|
| `binary_lp.npz` | `rivP`, `rivM` (bool), `valP`, `valM` (LP values), all 65536 genomes; written by `validate_binary.py` |
| `binary_margin.npz` | `margP`, `argP`: per binary genome, the minimum of W_DC − W_CD over the 65536 binary co-players and the co-player that attains it; read by `witness.py` |
| `binary_eff.npz` | `W` (65536 × 4): self-play limits of the binary genomes; `effLINE`: W_DD = 0. Read by `PRFull/Figures/fig1data.py` |
| `witnesses.npy` | the 19 co-players chosen by `witness.py` |
| `boxcover_W.out` | `boxcover.py` on `strata/W_strata.txt`: 214 strata, 85 prime boxes, a minimum cover of 41 boxes, each drawn as a 4 × 4 grid (rows: last outcome; columns: outcome before) |
| `strata/` | an early local search, 4 draws per stratum, one stratum index per line: `W_strata.txt` (214), `N_strata.txt` (4927), `ray_strata.txt` (2222, the ray T > S). Superseded by `st_eff.txt` |

## The census and the sweeps

Cannon arrays (IDs from the `JOBID` files; ledger lines in `PRFull/release/accounting/`):

| array | JOBID file | tasks | Slurm script | what |
|---|---|---|---|---|
| 45751662 | `../cannon/JOBID` | 512 | `../cannon/census.sbatch` | the stratum census of all 3¹⁶ strata |
| 45751663 | `../cannon/JOBID_SWEEP` | 256 | `../cannon/sweep.sbatch` | predicate-free sweep, 4 × 10⁶ random strategies per task |
| 45752239 | `../cannon/JOBID_DIM` | 384 | `../cannon/dim.sbatch` | maximum-dimension probe: six (class, sign) blocks × 64 stripes |
| 45831015 | `../cannon/alpha/JOBID_a333` | 512 | `../cannon/alpha/run_333.sbatch` | rarity sweep at α = 0.333, 3 190 000 strategies per task |

| file | what |
|---|---|
| `collect.sh` | needs `ssh cannon` and the remote kit `~/CL1/PartnersRivalsM2f/ca1`. Checks with `sacct`, task by task, that all 512 census tasks are COMPLETED and stops otherwise; copies the lines of `out/st_*.txt` whose class is not `-` into `st_eff.txt`; runs `census.py` into `census.out`; sums the `sweep n=` lines of `out/sw_*.err` into `sweep.out` |
| `st_eff.txt` | the census lines of the efficient strata with at least one confirmed rival, 1 524 791 lines `g class nP nM K`: stratum index, census class, number of draws confirmed as rivals for T > S and for T < S, number of draws |
| `census.py` | `python3 census.py st_eff.txt`: per locus (four wedges, two rays) the number of strata, the dimension histogram and the maximal strata; the binary-invisible strata (no binary friendly rival in the closed cell); and two exact consistency checks: (1) the C/D mirror s ↦ 15 − s must map wedge N onto wedge E; (2) the 0-dimensional strata must be exactly the binary friendly rivals, 8 / 1519 / 80 / 80 in W / S / N / E |
| `census.out` | the output of `census.py`. Check 2 passes in all four wedges. Check 1 holds only for a complete census: here N has 5020 strata and E 5022, the two extra strata of E are absent from mirror(N), and the output takes the union, 5022, as the lower bound |
| `sweep.out` | totals over the 256 tasks of 45751663 (1.024 × 10⁹ strategies): `effCC`, `effALT`, `effZ`, `rivP`, `rivM`, `FR_W`, `FR_S`, `FR_N`, `FR_E`, `rayP`, `rayM` |
| `dim.out` | 45752239 per (class, sign) block: the number of strata of dimension ≥ 5 with a confirmed rival, the largest dimension and the histogram. Its first line is the per-task `sacct` check (empty: all 384 tasks clean). `PRFull/Figures/sitable1.py` recomputes the six blocks from `../cannon/ca1_out/dim_*` |
| `sweep_a333.tsv` | one tab-separated line per task of 45831015: `ALPHA am n effCC effALT effZ rivP rivM FR_W FR_S FR_N FR_E rayP rayM` (α in thousandths, number of strategies, then counts; `FR_*` friendly rivals by wedge, `ray*` on the two rays). 512 × 3 190 000 = 1.633 × 10⁹ strategies; the paper's SI Table 3 |

The counts in `census.out` are those of the census alone. The paper's counts also take in strata found by the dimension
probe, by re-censuses and by independent recomputations, so some of them are larger (Methods; SI section "The
memory-two census" and SI Table 1). `PRFull/Figures/sitable1.py` recomputes every cell of SI Table 1 from `st_eff.txt`
and `../cannon/ca1_out/dim_*` (the census and the dimension probe together give 364 strata in W and 5022 in N).

In the release, the census output that stayed on Cannon is in `../cannon/ca1_out/` (`pt_*` confirmed points at full
double precision, `dim_*`, `sw_*`, `probe_45749304_*`), with the md5 of each file in `../cannon/ca1_out.md5`. The
stratum lines of all 512 census tasks are in `PRFull/release/census_class/census_st_45751662.tgz`.

## The evolution over the memory-two cube

Local copies of the output of `../cannon/wfc/wfc.c` (Slurm scripts, task files and `JOBID_*` files in
`../cannon/wfc/`). Each array was checked with `sacct` task by task.

| directory | array (`JOBID` file) | tasks | contents |
|---|---|---|---|
| `wfc_a333/` | 48366947 (`JOBID_a333`), `tasks_a333.txt` | 512 | N = 1000 run at α = 1/3: `r_1o3_<game>_1.txt`, ten `REP` lines each (5120 in all); `n_1o3_<game>.txt`, one `NEU` line of 20 000 draws |
| `wfc_a333/` | 48391602 (`JOBID_neu`), `tasks_neu.txt` | 84 | neutral refinement at 21 games: `nr_1o3_<game>_<chunk>.txt`, one `NEU` line of 250 000 draws each |
| `wfc_f2/` | 48724650 (`JOBID_f2`), `tasks_f2.txt` | 512 | N = 100 run (`wfc2`) at α = 1/3: `r_1o3_<game>_1.txt`, 5120 `REP` lines; no neutral draws (the analysis uses those of `wfc_a333/`) |
| `wfc_probe48329275/` | 48329275 (`JOBID_probe`), `tasks_probe.txt` | 88 | α = 1/2 at 8 games × 10 replicates (`r_500_<game>_<rep>.txt`), α = 0.333 at 8 games (`r_333_*`, `n_333_30.txt`). `sacct` marks tasks 82–88 FAILED: the move of a strategy dump failed after their `REP` lines were complete |
| `wfc_probe48343431/` | 48343431 (`JOBID_probe2`), `tasks_probe2.txt` | 28 | α = 1/2 at two games of wedge W (lines 30 and 38) × 10 replicates, for the exact counter of atom 111; α = 1/3 at 8 games (`r_1o3_*`, `n_1o3_*`) |
| `wfc_probe48723624/` | 48723624 (`JOBID_probe_f2`), `tasks_probe_f2.txt` | 48 | N = 100: α = 1/2 at 4 games × 10 replicates (`r_500_*`), α = 1/3 at 8 games (`r_1o3_*`) |

File names: `r_<α>_<game>_<first replicate>.txt`, `n_<α>_<game>.txt`, `nr_<α>_<game>_<chunk>.txt`; the α tag is `1o3`
(exactly 1/3), `333` (0.333) or `500` (1/2); `<game>` is the line of `../cannon/wfc/games.dat` (game id 1999 + line).
The columns, from the header of `wfc.c`:

```
REP gid ipt u v Emax pay ef alpha irep nsamp ndist a000 a001 a010 a011 a100 a101 a110 a111
    nhit nmis nnashscan nrivscan t_sim t_class inW ex111 nfr111 ncand111 total
NEU gid ipt u v Emax alpha M c000 c001 c010 c011 c100 c101 c110 c111 nnashscan nrivscan t
```

`aXYZ` counts the samples with efficient = X, Nash = Y, competitive = Z (`a101` must be 0); `cXYZ` counts the same
atoms among the M draws from the mutation law. In wedge W (`inW` = 1), `ex111` is the exact number of
individual-generations held by friendly rivals out of `total`.

| script | what |
|---|---|
| `wfc_binall.py <dir>` | compares the classifier of `wfc.c`, run as `./wfc binall <task> 500` on all 65536 binary genomes (files `B<task>.bin`, not kept here), with the exact atom table `PartnersRivals/Cannon/atoms512.bin`; prints the mismatches per game |
| `wfc_probe.py <dir> [F1\|F2]` | α = 1/2 check against the binary runs, whose per-replicate data are `PartnersRivals/Data/f1_e4.rep` (F1, N = 1000) and `f2_e4.rep` (F2, N = 100): per game the mean and replicate s.d. of efficiency and atom shares, Welch's t, and the most abundant atom. `python3 wfc_probe.py wfc_probe48329275`; `python3 wfc_probe.py wfc_probe48723624 F2` |
| `wfc_stats.py [F1\|F2]` | the summary statistics of the evolutionary runs on the 512 games for three strategy sets: the memory-one cube (`PartnersRivalsM1f`), the binary memory-two genomes (`PartnersRivals1`) and the memory-two cube (`wfc_a333/` or `wfc_f2/`): games by quadrant and wedge, shares, efficiency, most abundant and most enriched atoms, replicate agreement, cost. It reads `../Figures/figwfc.py`, `PartnersRivals1/FinalFigures/wfdata.py` and `PartnersRivalsM1f/Figures/cwdata.py` |

See the paper's SI section "The evolutionary maps" and the Computational Supplement; the kit itself is described in
`../cannon/wfc/README.md`.
