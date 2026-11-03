# PartnersRivalsM1f — programs and data for stochastic memory-one strategies

This directory holds the programs and data for memory-one strategies p = (p₁,p₂,p₃,p₄) in the cube [0,1]⁴.
The paper uses them in three places:

* SI section "Memory one": the checks of the closed forms (`check/`) and the two Cannon verification arrays
  (`cannon/`);
* Figure 2: its scripts in `PRFull/Figures/` (`fig2.py`, `fig2_check.py`, `fig2_check_f.py`) import
  `check/m1lib.py`, `check/working/partners_volume.py` or both;
* the memory-one row of Figures 4 and 5 (run F1) and of SI Figures 1 and 2 (run F2): the evolutionary runs in
  `cannon/cw/`, read through `Figures/cwdata.py` and `Figures/fig4c.py` by `PRFull/Figures/fig45.py`.

The Computational Supplement lists the programs in its Table 1 and the cost of every Cannon array in its Table 2.

## Conventions

* **Game.** (R,S,T,P) = (1,u,1+v,0). Hence T − S = 1 + v − u and 2R − (S+T) = 1 − u − v.
* **Wedges.** The lines u − v = 1 (T = S) and u + v = 1 (2R = S + T) cut the (u,v) plane into four wedges, named as
  in the code: W = {u+v < 1, u−v < 1}, S = {u+v < 1, u−v > 1}, N = {u+v > 1, u−v < 1}, E = {u+v > 1, u−v > 1}.
* **States.** CC = 0, CD = 1, DC = 2, DD = 3, with the focal player's move first. pᵢ is the probability of C after
  state i. The co-player sees the same state with the two moves swapped.
* **Binary strategies** are the 16 vertices of the cube, written p₁p₂p₃p₄ with 1 for C: `1010` is TFT, `0000` ALLD,
  `1111` ALLC. A binary memory-two co-player answers state j = 4·(last outcome) + (outcome before); bit j of its
  16-bit genome is 1 for C.
* **Error model.** A strategy p is played as ε + (1 − 2ε)p: each intended move is flipped with probability ε. Every
  property is the ε → 0 limit of the pair chain's stationary distribution. The limit is taken first, and the weak
  inequality is applied to it.
* **Names.** The code says "Nash" and "rival"; the paper says "stable" and "competitive".
* **Atoms.** Three bits E N R (efficient, Nash, rival), index 4E + 2N + R, written as the code `ENR` (e.g. `110`).
  The code `101` never occurs; the packer and the reader refuse data in which it does.
* **Strata.** Each coordinate is 0, 1 or interior (codes 0, 1, 2), so there are 81 strata. Stratum index =
  code(p₁) + 3·code(p₂) + 9·code(p₃) + 27·code(p₄).

## Layout

| path | contents |
|---|---|
| `check/` | Python checks of the memory-one closed forms (float, exact rational, symbolic), their reports (`*.out`), and the shared library `m1lib.py`; see `check/README.md` |
| `check/working/` | scratch derivations and one-off probes. They are not receipts. `partners_volume.py` (the graded partner measure) is imported by the Figure 2 scripts |
| `cannon/` | the kit of the two verification arrays: `m1f_pairs.c`, `compile.sh`, `task.sh`, four Slurm scripts, four `JOBID*` files; see `cannon/README.md` |
| `cannon/out/` | stdout and stderr of every task of those four arrays |
| `cannon/cw/` | the kit of the evolutionary runs: sources, build and Slurm scripts, task files, job IDs, per-cell output, packed tables; see `cannon/cw/README.md` and `cannon/cw/CHANGES.md` |
| `Figures/` | in this release only `cwdata.py` (reader of the packed tables; its docstring defines their formats) and `fig4c.py` (its `continuous()` builds the memory-one row of Figures 4 and 5) |

No compiled program is deposited; the build scripts make them.

## `check/`: local checks

```bash
cd check
python3 sweep_rivals.py       # float: rival sets on a grid, at random interior points, and on faces and edges
python3 exact_rivals.py       # exact rationals: the same on a rational grid, with the leading eps exponent
python3 strata_rivals.py      # all 81 strata against the 16 binary co-players, with witnesses
python3 symbolic_rivals.py    # sympy, eps a symbol: sign certificates cell by cell (Bernstein)
python3 sweep_efficiency.py   # float: the efficient sets and the friendly rivals per wedge
python3 exact_efficiency.py   # exact rationals: the same
python3 sweep_partners.py     # float: the Nash and partner criteria against the raw 16-deviator test
python3 exact_partners.py     # exact rationals: the same
python3 nash_field.py         # the dimension of the Nash set over a grid of games
```

Each script prints a report and writes no file; where a `.out` file of the same name exists, it is that report.
`strata_rivals.py` and `symbolic_rivals.py` need sympy; `nash_field.py` needs scipy. `m1lib.py` holds the pair
chain in floats and in exact fractions, payoffs, the rival, Nash and efficiency tests, and the closed-form
predicates (`in_riv_plus`, `in_riv_minus`, `in_eff_cc`, `in_eff_alt`, `in_eff_line`).

Numerical notes:

* Each float script prints the values of ε it uses. The float pair chain can lose accuracy at very small ε (errors
  of about 5e-6 were seen at ε = 1e-11 in the efficiency checks). Use exact fractions or `cannon/m1f_pairs.c` for
  anything decisive.
* `nash_field.py` decides limit-Nash from two values of ε: Nash if gap(1e-9) > −1e-6, or if gap(1e-9) > 0.2·gap(1e-7).

## `cannon/`: verification arrays

`m1f_pairs.c` computes the exact ε → 0 limit of the pair chain for a continuous memory-one strategy against binary
co-players. It uses GTH state reduction in leading-order (coefficient, exponent) arithmetic. GTH never subtracts,
so the leading coefficients carry no cancellation error; the code treats a value below 1e-12 as zero.

```bash
./compile.sh                                   # cc -O3, libm only
./m1f_pairs m1|m2 <n> <task> <ntasks>          # grid of spacing 1/n; co-players: 16 memory-one | 65 536 memory-two
./m1f_pairs rand|rand2 <N> <task> <ntasks>     # N random points drawn from all 81 strata; same co-players
./m1f_pairs fr|frand <n|N> <task> <ntasks>     # friendly rivals per wedge, on a grid | at random points
```

* Modes `m1`, `rand`, `m2`, `rand2` compare the computed rival margins with the closed-form rival sets (`pred_plus`,
  `pred_minus` in the source; the same predicates as `in_riv_plus`, `in_riv_minus` in `check/m1lib.py`). Each task
  prints one `SUMMARY` line and one `DISAGREE` line for each point where the two differ.
* Modes `fr`, `frand` decide efficiency from the computed self-play limit and rivalry from the computed margins
  against the 16 vertices, without any closed form. They count efficient strategies and friendly rivals per wedge.
  Each task prints one `FR-SUMMARY` line and one `IMPOSSIBILITY-HIT` line for each friendly rival found in W, N or E.
* Random points use splitmix64 seeded from the task index, so every sample can be regenerated.

`task.sh` maps an array index to segment = id / SHARDS and shard = id % SHARDS. The segments are 0 `m1`, 1 `rand`,
2 `m2`, 3 `rand2`, 4 `fr`, 5 `frand`. The Slurm scripts set the sizes:

| script | array | partition | what | job ID (file) |
|---|---|---|---|---|
| `riv1_probe.sbatch` | 0-3 | test | rate probe of segments 0-3 at reduced size | 45630088 (`JOBID_PROBE`) |
| `riv1.sbatch` | 0-511 | sapphire | segments 0-3, 128 shards each: the 1/500 grid (501⁴ points) and 10¹⁰ random points against the 16 vertices; the 1/30 grid (31⁴ points) and 10⁶ random points, each against the 65 536 binary memory-two strategies | 45636422 (`JOBID`) |
| `fr1_probe.sbatch` | 0-1 | test | rate probe of segments 4 and 5 | 45656466 (`JOBID_FR_PROBE`) |
| `fr1.sbatch` | 0-255 | sapphire | segments 4 and 5, 128 shards each: the 1/500 grid and 10¹¹ random points | 45658567 (`JOBID_FR`) |

Output: `cannon/out/<name>_<jobid>_<task>.out` and `.err`, with `<name>` = `probe`, `riv1`, `frprobe` or `fr1`.
The `.err` files hold the `/usr/bin/time -v` report of each task.

Receipts:

* `sacct` (in the release: `PRFull/release/accounting/sacct_table7.psv`): 45636422 has 512 of 512 tasks COMPLETED and
  45658567 has 256 of 256. Of the probes, 45656466 has 2 of 2 COMPLETED; one of the four tasks of 45630088 hit its
  time limit (TIMEOUT).
* The `out/` files were copied from Cannon and matched the Cannon copies by md5.
* Summing the task lines gives, for 45636422: 512 tasks, 73 003 425 522 strategies tested, `bad_plus` 0,
  `bad_minus` 0, no `DISAGREE` line. For 45658567: 256 tasks, 163 001 502 001 strategies, FR_W = FR_N = FR_E = 0,
  FR_S = 8 543 035 103, no `IMPOSSIBILITY-HIT` line. To recompute:

```bash
cd cannon/out
cat riv1_45636422_*.out | awk '/^SUMMARY/ {for (i = 2; i <= NF; i++) {split($i, a, "="); s[a[1]] += a[2]}; n++}
  /^DISAGREE/ {d++} END {print n, "tasks; tested", s["tested"], "bad_plus", s["bad_plus"], "bad_minus", s["bad_minus"], "DISAGREE lines", d+0}'
cat fr1_45658567_*.out | awk '/^FR-SUMMARY/ {for (i = 2; i <= NF; i++) {split($i, a, "="); s[a[1]] += a[2]}; n++}
  /^IMPOSSIBILITY-HIT/ {h++} END {print n, "tasks; tested", s["tested"], "FR_W", s["FR_W"], "FR_N", s["FR_N"], "FR_E", s["FR_E"], "FR_S", s["FR_S"], "hits", h+0}'
```

## `cannon/cw/`: evolutionary runs

`cw.f` is a Wright–Fisher program for memory-one strategies anywhere in the cube. It runs the process of the
binary program `sf1.f` of the project DiskM1WF (byte copy `ref_sf1.f`; `CHANGES.md` lists every change). Each
lineage carries its own 4-vector. Mutants come from a kernel that sets each coordinate independently to ε, to
1 − ε, or to a uniform interior value, with probability 1/3 each, so every stratum has mass 1/81. When a lineage
dies, and at the end of each replicate, its abundance is credited to its exact ε → 0 atom and to its stratum.

| file | role |
|---|---|
| `cw.f` | the program of run F1 |
| `cw2.f`, `mkcw2.py` | run F2: `cw2.f` is generated from `cw.f` by `mkcw2.py`; do not edit it |
| `ref_sf1.f`, `ref_sf2.f` | byte copies of the binary programs `sf1.f` and `sf2.f`, for the audit build |
| `m1cls.c`, `m1cls.h` | the classifier: atom and stratum of a strategy at a game. E and R from the closed forms, decided exactly; N from the 16-deviator test with the leading-order GTH of `../m1f_pairs.c` |
| `kappa.c` | the kernel's mass on each atom and stratum at each game (M = 10⁷ draws, seed 41000000 + isl) |
| `mt.f`, `mtb.f`, `payf2.f`, `payf3.f`, `games.dat` | byte copies of the binary run's random-number, payoff and game files |
| `compile.sh` | builds every target (below) |
| `task.sh` | runs one cell, idempotently: `bash task.sh <isl> <rundir> <exe>` |
| `probe.slurm`, `prod.slurm`, `kappa.slurm`, `probe2.slurm`, `prod2.slurm` | the Slurm scripts |
| `submit.sh`, `verify.sh` | submission (writes `JOBID_<tag>`, `IDX_<tag>` and the ledger line) and per-task completion from `sacct` |
| `IDX_<tag>` | the task files: the array indices of each submission |
| `pack.py` | collects a run directory and `kappa/` into the packed tables |
| `check_v4.py` | invariant checks on any output directory |
| `test_cls.py`, `test_cls.c`, `test_fcall.f`, `test_fbind.f90` | tests of the classifier and of the Fortran–C link |
| `quarantine/m1cls_stub.c` | an old test stub; the production targets refuse it |
| `build/sf1t.f`, `build/sf2t.f` | `ref_sf1.f` and `ref_sf2.f` with a short run length, written by `compile.sh test` / `test2` |

The two runs (ε = 10⁻⁴, 10 replicates per game, the second half of each replicate sampled):

| run | N | β | µ | generations | seeds | program | output | figures |
|---|---|---|---|---|---|---|---|---|
| F1 | 1000 | 100 | 10⁻² | 10⁷ | 40000000 + (isl−1)·10 + irep | `cw.x` | `run1/` | Figures 4 and 5 |
| F2 | 100 | 3 | 10⁻⁴ | 10⁸ | 42000000 + (isl−1)·10 + irep | `cw2.x` | `run2/` | SI Figures 1 and 2 |

The games are the 512 rows of `games.dat`, the same games as the binary runs. Cell ipt (0-511) is task
isl = 3·ipt + 3, so each production array is `--array=3-1536:3`.

Build:

```bash
bash compile.sh            # cw.x, cwb.x (BINMODE audit build), kappa.x
bash compile.sh probe      # cwp.x: cw.x with 10^6 generations, for the F1 rate probe
bash compile.sh f2         # cw2.x, cwb2.x
bash compile.sh probe2     # cwp2.x: cw2.x with 10^7 generations, for the F2 rate probe
```

On Cannon, after `module load intel/25.2.1-fasrc01`, Fortran: `ifx -O3 -xHost -fp-model precise`; C:
`icx -O2 -xHost -fp-model precise -ffp-contract=off`. Without ifx, Fortran:
`gfortran -O2 -ffp-contract=off -fno-fast-math`; C: `gcc -O2 -ffp-contract=off`. Never build `m1cls.c` with fast-math: its exact sums need strict IEEE arithmetic.

Arrays (job names `PartnersRivalsM1f.cw*`, comment `PartnersRivalsM1f`):

| script | array | partition | output | job ID (file) | sacct |
|---|---|---|---|---|---|
| `probe.slurm` | 3, 1527, 1530, 1533, 1536 | test | `probe/` | 48340155 (`JOBID_probe`) | 5 of 5 COMPLETED |
| `kappa.slurm` | one task, all 512 cells | sapphire | `kappa/` | 48340156 (`JOBID_kappa`) | 1 of 1 COMPLETED |
| `prod.slurm` | 3-1536:3 | sapphire | `run1/` | 48366968 (`JOBID_prod`, also `JOBID`) | 512 of 512 COMPLETED, 407.6 core-h |
| `probe2.slurm` | 3, 1527, 1530, 1533, 1536 | test | `probe2/` | 48373158 (`JOBID_probe2`) | 5 of 5 COMPLETED |
| `prod2.slurm` | 3-1536:3 | sapphire | `run2/` | 48408751 (`JOBID_prod2`) | 512 of 512 COMPLETED, 404.5 core-h |

Slurm logs are in `log/` (`o_*` and `e_*` for `probe`, `prod`, `kappa`, `f2probe`, `f2prod`). F2 uses the kernel
masses of F1: they depend only on the game and on ε.

Output per cell (in `run1/`, `run2/`, `probe/`, `probe2/`; `games.dat` is staged beside them):

* `<isl>`: one line per replicate: `gid ie u v Emax pay ef meanlive solves_per_gen share0 … share7`, the shares
  indexed 4E + 2N + R.
* `h<isl>`: a header (cell, game, run parameters, the atoms of the 16 corners), then the exact counts: 8 atom lines,
  81 stratum lines, 16 corner lines and `# total`. A cell is complete only when this tail is present.
* `w<isl>`: classifier diagnostics per replicate (lineages classified, classifier calls, near ties, cache hit rate).
* `kappa/k<isl>`: one line of 98 fields: `isl gid ie eps u v M`, 8 atom fractions, 81 stratum fractions, and the
  fractions of the two kinds of near tie.

Collection:

```bash
python3 pack.py run1 kappa packed              # F1
python3 pack.py run2 kappa packed --run=F2     # F2; kappa_e4.tsv is shared with F1
```

`pack.py` checks every cell: the task decode against `games.dat`, the run parameters, 8 + 81 + 16 count lines, atom
and stratum totals equal to replicates × (generations/2) × N, atom `101` zero, 10 replicate lines whose shares sum to
1, and the kappa lines. It exits with status 1 on any failure. It writes to `packed/`:

| file | one row per | columns |
|---|---|---|
| `cw1_e4.tsv`, `cw2_e4.tsv` | cell | `ipt gid u v emax pay ef s000 s100 s010 s001 s110 s011 s111 s101 meanlive` |
| `cw1_e4.rep.tsv`, `cw2_e4.rep.tsv` | replicate | `ipt irep pay ef meanlive solves_per_gen s000 … s101` |
| `cw1_e4.diag.tsv`, `cw2_e4.diag.tsv` | cell | `ipt ncls ncall nneartie neartie_abund neartie_share nhw_max exacttie_share` |
| `kappa_e4.tsv` | cell | `ipt k000 … k101 z00 … z80` (atom and stratum masses of the kernel) |

Here s is the share of the sampled population in an atom, ef = pay/emax, and emax = max(1, (1+u+v)/2).
`Figures/cwdata.py` defines these formats in full. Re-running both `pack.py` commands on the deposited `run1/`,
`run2/` and `kappa/` (2026-10-17) printed "all 512 cells complete and consistent" for each run and gave all seven
files of `packed/` byte for byte. The classifier's near-tie flag (a non-zero Nash margin within
10⁻⁹·max(1, |u|, |1+v|) of zero) was raised for 341 lineages in F1 and 4 in F2, at most 2.1·10⁻⁸ and 10⁻¹⁰ of
any cell's abundance.

The directories `run1/`, `run2/`, `kappa/`, `probe/`, `probe2/`, `log/`, `vcheck/` and `vcheck2/` were copied from
Cannon and matched the Cannon copies by md5.

The statistics of the memory-one row of Figures 4 and 5 and of SI Figures 1 and 2 are printed by
`PRFull/Figures/fig45.py` into `PRFull/Figures/fig45_stats.txt`.

Validation outputs and tests (not data):

* `vcheck/`, `vcheck2/`: runs on Cannon for F1 and F2. In `cb/` and `sf/`, the BINMODE build (the 16 binary
  strategies through `cw.f`'s machinery) and the binary reference program, at cells 243 and 1203; their outputs are
  byte-identical. In `vcheck/c/`, the corner atoms of the linked program at all 512 cells (2 generations).
* `run_local/`, `kappa_local/`, `packed_local/`, `run_local2/`, `packed_local2/` and the two `run_local*.log`: a short
  local run of the whole chain (10⁴ generations or fewer, M = 10⁶). `cwdata.py` refuses these files as production
  data.

`cannon/cw/README.md` gives the full protocol, the validation tables and the cost measurements.
