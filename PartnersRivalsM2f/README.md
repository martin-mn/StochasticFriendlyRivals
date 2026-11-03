# PartnersRivalsM2f: the memory-two strategy cube

Programs and data for memory-two strategies with continuous coordinates, p ∈ [0,1]¹⁶, taken in the limit of
vanishing implementation error. The 65536 binary memory-two strategies are the corners of this cube. The binary
case is the project `../PartnersRivals`; the memory-one cube is `../PartnersRivalsM1f`.

The paper uses this directory in the SI sections "Memory two", "The memory-two census" (SI Table 1),
"Certificates" and "The evolutionary maps", and in the Computational Supplement, sections "The limit computation
and its tolerances" and "The census, its campaigns and the rarity sweep".

## Conventions

They are those of the paper's Methods, in short:

- Game: (R,S,T,P) = (1,u,1+v,0). The switch line is u+v = 1; the line T = S is u−v = 1.
- Wedges: W = {u+v<1, T>S}, S = {u+v<1, T<S}, N = {u+v>1, T>S}, E = {u+v>1, T<S}. The switch line splits into
  two rays, u < 1 (T>S) and u > 1 (T<S).
- Outcomes are coded with the player's own action first: CC = 0, CD = 1, DC = 2, DD = 3.
- A state is s = 4·(last outcome) + (outcome before it), s = 0..15. p_s is the probability of C in state s.
- A binary genome is a 16-bit integer whose bit s is 1 for C (ALLC = 65535, ALLD = 0).
- Error model: both players play ε + (1−2ε)p in place of p. Every quantity is the leading term as ε → 0,
  computed exactly by GTH elimination on (coefficient, power of ε) pairs. There is no threshold in ε and no
  extrapolation.
- W = (W_CC, W_CD, W_DC, W_DD) is the limit frequency of the four outcomes.
- Self-play class, read with a tolerance of 10⁻¹²: C if W_CC = 1, A if W_CD + W_DC = 1 (in self-play
  W_CD = W_DC = ½), Z if W_DD = 0 without C or A, `-` otherwise. Efficiency means class C below the switch
  line, class A above it, and W_DD = 0 on it (classes C, A and Z).
- Rival margin: sign·(W_DC − W_CD) minimised over the 65536 binary co-players, with sign +1 for T > S and
  −1 for T < S. A margin ≥ −10⁻⁹ counts as a rival. The paper's SI section "Co-players and deviators of any
  memory" explains why the binary co-players suffice.
- Support stratum: each p_j is pinned to 0, pinned to 1, or free in (0,1). The stratum index is
  g = Σ d_j·3^j with d_j = 0 (p_j = 0), 1 (free), 2 (p_j = 1). A pattern string lists p_0 … p_15 from left to
  right, with `*` for a free coordinate, e.g. `1*1*001*1*1**0*0`. There are 3¹⁶ = 43 046 721 strata.

## Contents

| path | what |
|---|---|
| `cannon/` | the census kit: engine, census, Slurm scripts, JOBID files; `cannon/README.md` |
| `cannon/alpha/` | the kit of the rarity sweep at α = 0.333 and its α-grid probe |
| `cannon/wfc/` | the kit of the evolutionary runs over the cube; `cannon/wfc/README.md` |
| `cannon/archive/` | `cannon_2026-09-26.tgz`: the raw output of the N = 1000 evolutionary run and its probes, with the sampled strategies, and the alpha kit with its output, as they were on Cannon; `cannon_out_f2_2026-09-26.tgz`: the raw output of the N = 100 run, with the sampled strategies (see `cannon/wfc/README.md`) |
| `cannon/ca1_out/` | added by the release: the census output that was kept only on Cannon (below) |
| `cannon/ca1_out.md5` | md5 of each file of `ca1_out/`, taken on Cannon |
| `check/` | local programs, validation and the collected results; `check/README.md` |
| `check/working/` | from the derivation pass's scratch, which the release otherwise leaves out: the programs and output of five checks (`t2_fw.c`, `effcc_task1/effcc.c` with its logs, `t2_altform2.c`, `verify_task1/`, and `task2/` with the kit of the re-censuses, `t2.c` and `t2b.c`); `check/working/README.md` |
| `check/strata/` | early local lists of rival strata (wedge W, wedge N, ray T>S), 4 draws per stratum; superseded by the census |
| `check/wfc_a333/`, `check/wfc_f2/` | output of the evolutionary runs at N = 1000 and N = 100 |
| `check/wfc_probe48329275/`, `check/wfc_probe48343431/`, `check/wfc_probe48723624/` | output of the three evolutionary probes |
| `Figures/figwfc.py` | loader and plots of the evolutionary output; imported by `PRFull/Figures/fig45.py` |

### `cannon/`

| file | what |
|---|---|
| `m2f.c` | the engine: `freqs(p,q,W)` (four outcome frequencies), `freqs16` (sixteen state frequencies), `rival_margin(p,sign,&arg,earlyquit)` (margin and arg-min co-player) |
| `census.c` | modes `strata` (all 3¹⁶ strata) and `sweep` (random strategies, no closed form used) |
| `witnesses.h` | 157 binary co-players used as a screen before the full scan (see Checks) |
| `dimprobe.c` | for one (class, sign) pair, every stratum with at least DMIN free coordinates; reports the largest number of free coordinates with a confirmed rival |
| `confirm.c` | full 65536-co-player confirmation for a list of strata |
| `onestratum.c` | the per-stratum logic of `dimprobe` on one stratum, verbose |
| `strata.c` | self-play class of all 3¹⁶ strata at two interior draws each (local) |
| `strata2.c` | the same plus a screen against 19 co-players (local) |
| `wit2.c` | arg-min co-players over random continuum strategies, used to grow the witness set |
| `m2f_main.c` | command-line driver: binary self-play limits, binary rival test, one binary pair |
| `m2f_py.c` | the C side of the ctypes driver `check/m2f.py` |
| `compile.sh`, `task.sh` | build and array task of the census and the sweep |
| `census_probe.sbatch`, `census.sbatch`, `sweep.sbatch`, `dim.sbatch` | Slurm scripts (below) |
| `JOBID_PROBE`, `JOBID`, `JOBID_SWEEP`, `JOBID_DIM` | the array job IDs, written at submission |

### `check/`

| file | what |
|---|---|
| `m2f.py` | ctypes driver for `../cannon/libm2f.dylib` |
| `exact.py` | the engine again in exact `Fraction` arithmetic |
| `m2lib.py` | `occupation_lp(p, sense)`: the exact minimum of sense·(W_DC − W_CD) over the occupation-measure polytope of the Markov decision process whose controller is the co-player |
| `validate_binary.py` | runs that LP on all 65536 binary genomes; writes `binary_lp.npz` |
| `witness.py` | greedy set cover of the binary non-rivals by co-players; writes `witnesses.npy` |
| `boxvertex.py` | numerical search for a continuum co-player that does worse than the best binary one |
| `certify7.py`, `certify8N.py`, `certify11S.py`, `certify_segment.py` | exact rational certificates of single points (see Checks) |
| `explore_W.py` | starts from each of the eight binary friendly rivals of wedge W, moves one coordinate at a time over a grid, and marks efficiency and rivalry |
| `boxcover.py`, `boxcover_W.out` | minimum cover of a set of strata by three-valued boxes, and its output on the early wedge-W list |
| `collect.sh` | collects the census and the sweep from Cannon (below) |
| `census.py`, `census.out` | analysis of the census and its output |
| `st_eff.txt` | the census lines of efficient strata |
| `sweep.out`, `dim.out` | totals of the sweep and of the dimension probe |
| `sweep_a333.tsv` | the 512 task lines of the rarity sweep at α = 0.333 |
| `binary_eff.npz` | `W` (65536 × 4, self-play limit of each binary genome) and `effLINE` (W_DD = 0) |
| `binary_margin.npz` | `margP`, `argP`: rival margin for T > S of each binary genome and its arg-min co-player |
| `binary_lp.npz` | `rivP`, `rivM` (LP verdicts), `valP`, `valM` (LP minima) |
| `wfc_binall.py`, `wfc_probe.py`, `wfc_stats.py` | checks and statistics of the evolutionary runs (`cannon/wfc/README.md`) |
| `rebuild_refs.py` | a bibliography helper for the manuscript sources, which are not in this release; no result depends on it |

## Building

    cd cannon && ./compile.sh            # census, with cc -O3 -march=native

`compile.sh` builds only `census`. Each other program in `cannon/` is a single C file that calls the engine, so it
is compiled together with `m2f.c` in the same way. `check/m2f.py` loads a shared library. Build it as below;
the universal binary also serves an x86_64 Python on an arm64 Mac:

    cd cannon && cc -O3 -arch arm64 -arch x86_64 -fPIC -shared -o libm2f.dylib m2f.c m2f_py.c -lm

`cannon/alpha/compile.sh` and `cannon/wfc/compile.sh` build their own programs in their directories.

## Running

| command | output |
|---|---|
| `census strata <stripe> <nstripes> <Keff> <st_out> <pt_out> [Kother]` | strata g ≡ stripe (mod nstripes); Keff draws per efficient stratum, Kother (default 8) otherwise |
| `census sweep <n> <seed> <out>` | n random strategies: half uniform on the cube, half with each coordinate 0 or 1 with probability 0.36 each and uniform otherwise |
| `alpha/census psweep <n> <seed> <alpha_permille> <out>` | n strategies, each coordinate 0 with probability α, 1 with probability α, uniform otherwise |
| `dimprobe <dmin> <K> <sign> <C\|A\|Z> <stripe> <nstripes>` | one line `g nfree class` per stratum with a confirmed rival |
| `confirm <stratumlist> <sign> <K> <stripe> <nstripes>` | margin and class of K draws per listed stratum |
| `onestratum <g> <K> <sign>` | verbose trace for one stratum |
| `strata <stripe> <nstripes> [out]`, `strata2 <stripe> <nstripes> [K] [out]` | class lines (`C g`, `A g`, `Z g`), or `g class nP nM` for strata with screened survivors |
| `wit2 <n> <seed>` | `argP margP argM margM` per random strategy |
| `m2f self\|rival\|pair ...` | binary self-play limits, binary rival margins, or one binary pair |

Interior coordinates are drawn from a seeded xorshift generator; in the census the seed of each stratum is a
fixed function of g, so every stratum's draws can be repeated.

## Output formats

- `st_*.txt` and `check/st_eff.txt`: one line `g cls nP nM K` per stratum with at least one confirmed rival:
  stratum index, self-play class at the first draw, number of draws confirmed as rivals for T > S and for T < S,
  and the number of draws. `st_eff.txt` keeps the lines of class C, A or Z (1 524 791 lines).
- `pt_*.txt`: `P g cls p_0 … p_15` (T > S) or `M g cls p_0 … p_15` (T < S) at full double precision (`%.17g`),
  up to four points per efficient stratum and sign.
- `sw_*.txt`: the points of the sweep that are efficient below the switch line and rivals for T > S. Each sweep
  task writes its totals as one line on stderr (`sweep n=...: effCC ... effALT ... effZ ... rivP ... rivM ... |
  FR_W ... FR_S ... FR_N ... FR_E ... | rayP ... rayM ...`). `check/sweep.out` sums them over the tasks.
- `dim_<class>_<sign>_<stripe>.txt`: the `dimprobe` lines of one task. `check/dim.out` gives, per (class, sign)
  block, the number of strata with at least five free coordinates, the largest number found and a histogram.
- `psweep` writes one tab-separated line to stdout:
  `ALPHA am n effCC effALT effZ rivP rivM FR_W FR_S FR_N FR_E rayP rayM`, with α = am/1000.
  `check/sweep_a333.tsv` holds the 512 lines of array 45831015.
- FR_W and FR_S count draws of class C that are rivals for T > S and for T < S; FR_N and FR_E do the same for
  class A. rayP and rayM count draws of class C, A or Z that are rivals for T > S and for T < S.

## Cluster runs

All ran on Harvard's Cannon cluster, one core per task. The census kit was in `~/CL1/PartnersRivalsM2f/ca1/`,
the alpha kit in `~/PartnersRivalsM2f_alpha`.

| script | partition | tasks | JOBID file | array | what |
|---|---|---|---|---|---|
| `census_probe.sbatch` | test | 4 | `JOBID_PROBE` | 45749304 | cost per unit: 4 of 20 000 stripes |
| `census.sbatch` | sapphire, 12 h | 512 | `JOBID` | 45751662 | all 3¹⁶ strata, Keff = 64, Kother = 8 |
| `sweep.sbatch` | sapphire, 12 h | 256 | `JOBID_SWEEP` | 45751663 | 256 × 4 × 10⁶ = 1.024 × 10⁹ random strategies, seed 1000 + task |
| `dim.sbatch` | sapphire, 12 h | 384 | `JOBID_DIM` | 45752239 | `dimprobe 5 512` for six (class, sign) pairs × 64 stripes |
| `alpha/probe.sbatch` | sapphire, 1 h | 8 | `alpha/JOBID_probe` | 45830542 | α-grid sizing probe, 10⁵ draws per task at α = 0, 0.1, 0.2, 0.3, 0.333, 0.36, 0.4, 0.5 (from the submission ledger) |
| `alpha/run_333.sbatch` | sapphire, 8 h | 512 | `alpha/JOBID_a333` | 45831015 | 512 × 3 190 000 draws at α = 0.333, seed 100000 + 1000·am + task |

In the census, the screen against `witnesses.h` runs first and every survivor is confirmed against all 65536
co-players. The six `dim.sbatch` blocks are classes C, A and Z, each with sign +1 and −1: wedges W, S, N, E and
the two switch-line rays.

The census was sized from its probe, not from a local timing: 210–270 s per 1/20 000 of the strata on Cannon,
about 1330 core-hours in all. `cannon/alpha/run_250.sbatch` and `run_450.sbatch` (α = 0.25 and 0.45) have no
JOBID file. The arrays of the evolutionary runs are listed in `cannon/wfc/README.md`.

The paper's rarity sweep (SI Table 3) is a rerun at exactly α = 1/3 in `PRFull/cannon/sweep13`, because
`psweep` reads α in thousandths. `PRFull/release/reproduce_figures.sh` compares that rerun with
`check/sweep_a333.tsv`.

### Collection

`check/collect.sh` runs `sacct` on the census array and stops unless every task is COMPLETED. It then copies the
census lines of class C, A or Z into `check/st_eff.txt`, runs `census.py` on them (output `census.out`), and sums
the sweep's stderr lines into `sweep.out`. `cannon/alpha/collect.sh <JOBIDfile>` does the same check for an alpha
array and sums its `ALPHA` lines per α.

Where the release keeps the rest of the output:

- `PRFull/release/census_class/census_st_45751662.tgz`: the stratum lines of all 512 census tasks, also those
  of class `-`, with the task logs and an md5 list.
- `cannon/ca1_out/`: `pt_*.txt` (512 files), `dim_*` (the 384 `dimprobe` outputs and their Slurm logs),
  `sw_*` (the 256 sweep outputs and their Slurm logs) and `probe_*` (the Slurm logs of the census probe),
  2440 files in all. Each matches `cannon/ca1_out.md5`.

## Checks

Job states. Every task of 45749304, 45751662, 45751663, 45752239, 45830542 and 45831015 is COMPLETED in `sacct`
(`PRFull/release/accounting/sacct_table7.psv`).

The engine:

- the self-play limits of all 65536 binary genomes agree byte for byte with the binary engine
  `PartnersRivals/check/pairs.c`;
- the two binary rival sets agree with that project's element by element (2640 and 2640), not only in count;
- `exact.py` agrees with the C engine on 200 random binary pairs and on continuum pairs.

The co-players:

- `validate_binary.py`: the LP marks 2144 binary strategies, which are the `defP`/`defM` masks of
  `PartnersRivals/Figures/m2_masks.npz`, not the 2640 rivals of `rivP`/`rivM` in the same file.
- `witness.py`: 19 co-players (`witnesses.npy`, the first 19 entries of `witnesses.h`) already give a negative
  margin to each of the 62 896 binary strategies that are not rivals for T > S. The other 138 entries of
  `witnesses.h` are the most frequent arg-min co-players that `wit2` found over 22 125 random continuum
  strategies. A point that fails any of them is not a rival, so the screen never rejects a rival.
- `boxvertex.py` tests 24 rivals. For each it tries 4000 random co-players q ∈ [0,1]¹⁶ and 40 Nelder–Mead
  runs from random starts. None gets more than 1.7 × 10⁻¹⁶ below the binary minimum.

The census (`check/census.out`):

- Check 2 passes: in each of W, S, N and E the strata with no free coordinate coincide, as sets, with the
  binary friendly rivals of `m2_masks.npz` (8, 1519, 80 and 80 of them).
- Check 1 compares wedge N with the C/D mirror image of wedge E. It finds 2 strata of E whose images the census
  missed in N; with them N has 5022 strata, as E has.
- The counts in `census.out` come from the census as run, with a finite number of draws per stratum, so they are
  lower bounds. The final counts are in the paper's SI section "The memory-two census" and SI Table 1, which
  `PRFull/Figures/sitable1.py` recomputes from the census and the dimension probe.

Exact certificates (Python `Fraction`, no floating point; each scans all 65536 co-players on 12 processes):

| script | point | prints |
|---|---|---|
| `certify7.py` | (1, ½, 1, ½, 0, 0, 1, ½, 1, ¼, 1, ¼, ¼, 0, ¼, 0), stratum `1*1*001*1*1**0*0` | self-play limit, W_CC = 1, the T > S margin, its arg-min and the number of co-players at margin 0 |
| `certify8N.py` | (⅔, 1, 0, ⅓, ⅓, ⅓, 0, 0, ⅔, 1, ⅔, ⅓, ⅓, 0, 1, 0) | the same, with W_CD + W_DC = 1 in place of W_CC = 1 (its docstring repeats that of `certify7.py`; the code checks this point) |
| `certify11S.py` | (1, ½, ½, ½, 1, ½, ½, ½, 1, ½, ½, 1, 0, ½, ½, ½), stratum `1***1***1**10***` | the same for T < S and W_CC = 1 |
| `certify_segment.py` | p(t) = (0,0,0,t, 0,0,0,0, 0,1,0,1, 1,0,1,0) | self-play limits at t = 0, 1 and at five interior t; the T > S margin at t = ½ |

`boxcover.py` grows prime boxes in exponential time. On the 4927 wedge-N strata of `check/strata/` it did not finish
in 15 minutes; use it on the wedge-W list only.
