# `cannon/wfc`: evolution over the memory-two cube [0,1]^16

`wfc.c` simulates a Wright–Fisher process with pairwise-comparison imitation over continuous memory-two strategies. It
classifies every sampled resident exactly, in the limit ε → 0, with the engine `m2f.c`. Its output is the
memory-two stochastic row of the paper's Figures 4 and 5 (N = 1000) and of SI Figures 1 and 2 (N = 100), which
`PRFull/Figures/fig45.py` draws through `../../Figures/figwfc.py`.

One source, two builds (`compile.sh` makes both):

| binary | run | N | β | μ | ε | generations per replicate | samples per replicate |
|---|---|---|---|---|---|---|---|
| `wfc` (default build) | F1 | 1000 | 100 | 10⁻² | 10⁻⁴ | 10⁷ | 1000 times × 4 individuals |
| `wfc2` (`-DNPOP=100 -DBETAV=3.0 -DMUV=1e-4`) | F2 | 100 | 3 | 10⁻⁴ | 10⁻⁴ | 10⁸ | 4000 times × 1 individual |

Each run has 10 replicates per game at 512 games. The second half of each replicate is sampled. The process is that
of the binary runs F1 and F2 of the companion paper on binary strategies (Fortran `sf.f`, project `DiskM2WF`, not part
of this release; see also arXiv:2608.06147), with three changes: the mutation law, strategies stored as vectors, and
abundance estimated by sampling (below).

## Files

| file | content |
|---|---|
| `wfc.c` | the simulation and the modes `run`, `neutral`, `binall` |
| `m2f.c` | exact ε → 0 limit of the stationary distribution of a pair of strategies in [0,1]^16 (leading-order GTH reduction) |
| `witnesses.h` | 157 binary co-players tried first in the rivalry screen (survivors get the full 65536 scan) |
| `games.dat` | the 512 games, one per line (format below) |
| `compile.sh` | builds `wfc` and `wfc2` (the binaries are not deposited) |
| `task.sh` | one `wfc run` task per array element, and optionally a `wfc neutral` |
| `task_neu.sh` | one `wfc neutral` chunk per array element |
| `probe.sbatch`, `probe2.sbatch`, `probe_f2.sbatch` | probe arrays |
| `run_a333.sbatch`, `run_f2.sbatch`, `run_neu.sbatch` | production arrays |
| `tasks_*.txt` | task files, one line per array element, first line a `#` comment naming the columns |
| `JOBID_*` | the array job ID of each array, from `sbatch --parsable` |
| `S1.bin`, `S1_1.bin`, `S1_3.bin` | strategy dumps of test runs of task 1 (record format below) |

## Conventions

* **Game.** Payoffs (R, S, T, P) = (1, u, 1 + v, 0). `games.dat` columns: game ID `gid` (2000 + `ipt`), an unused
  integer, R, S, T, u, v, tag `V_p<ipt>_<quadrant><wedge>` (quadrant `HA`, `PD`, `SD`, `SH`; wedge `W`, `N`, `E`, `S`).
  Line t of the file is task t, with `ipt` = t − 1. The program requires R = 1. Emax = max(R, (S + T)/2).
* **States.** j = 4·(last outcome) + (outcome before). Outcome o = 2·(own defection) + (co-player's defection):
  CC = 0, CD = 1, DC = 2, DD = 3, own action first. p[j] is the probability to cooperate in state j. A binary genome
  has bit j = 1 for C (ALLC = 65535, ALLD = 0).
* **Errors.** Each player plays ε + (1 − 2ε)p. The process uses ε = 10⁻⁴; the classification is the limit ε → 0.
* **Mutation.** Each of the 16 coordinates is drawn independently: 0 with probability α, 1 with probability α, uniform on
  (0, 1) otherwise. The initial population is drawn the same way. α = 1/2 is uniform mutation over the 65536 binary
  genomes, the kernel of the binary runs. The production runs use α = 1/3. Pass α as an exact fraction (`1/3`) or in
  permille (`500`).
* **Atoms.** Code = 4·efficient + 2·Nash + competitive, written `aXYZ` (X efficient, Y Nash, Z competitive). Atom 101
  must be empty; the program stops if it is not.
* **Classification** (ε → 0, tolerance 1e-9):
  * efficient: in self-play the limit has W_CC = 1 below the line u + v = 1, and W_CD = W_DC = 1/2 above it;
  * competitive: sign(T − S)·(W_DC − W_CD) ≥ 0 against each of the 65536 binary co-players;
  * Nash: π(q, p) ≤ π(p, p) for each of the 65536 binary co-players q.

  Binary co-players suffice: the extreme co-player over the cube sits at a vertex of the box. For atom 111 the Nash scan is skipped, except for the first
  20 such strategies of each task, which are scanned and must come out Nash (the program stops otherwise).

## Building and running

    bash compile.sh                      # wfc (F1) and wfc2 (F2)

    ./wfc run <task> <alpha> <itend> <NT> <NI> [nrep] [seedbase] [rep0] [dumpfile]
    ./wfc neutral <task> <alpha> <M> [seed]
    ./wfc binall <task> 500

All modes read `games.dat` in the working directory; `<task>` is a line number of it (1-based).

**`run`** simulates replicates rep0 … rep0 + nrep − 1 (default 1 … 10) for `itend` generations. Replicate r of task t
has the seed seedbase + 10(t − 1) + r (default seedbase 9000000). In the second half, NT sample times are spread evenly
with a random offset, and NI individuals are drawn at each. The efficiency is accumulated exactly over the half, not
sampled. A replicate's result does not depend on how many replicates share its task: strategies carry never-reused
serials, and the pair cache (2^20 slots) only decides whether a value is recomputed.

In wedge W (u + v < 1 and u − v < 1) atom 111 is counted exactly instead of sampled. A mutant with p₀ = 1 and p₁₅ = 0
is tested in self-play at birth. Its individual-generations are counted, and the rivalry test at its death decides
whether they count for 111.

Output on stdout, one line per replicate (F2 files begin with a `#` line naming the build values; the F1 files in
`../../check/wfc_a333/` have no such line):

    REP gid ipt u v Emax pay ef alpha irep nsamp ndist a000 a001 a010 a011 a100 a101 a110 a111
        nhit nmis nnashscan nrivscan t_sim t_class inW ex111 nfr111 ncand111 total

* `pay`: mean payoff over the sampled half; `ef` = pay/Emax.
* `nsamp`: samples (NT·NI); `ndist`: distinct sampled strategies; `aXYZ`: samples in each atom.
* `nhit`, `nmis`: pair-cache hits and misses; `nnashscan`, `nrivscan`: full scans done; `t_sim`, `t_class`: seconds.
* `inW` = 1 in wedge W. There `ex111` is the exact number of individual-generations held by atom 111, from `nfr111`
  distinct strategies out of `ncand111` candidates tested, and `total` = (itend/2)·N.

The strategy dump (`dumpfile`, default `S<task>_<rep0>.bin`) holds every distinct sampled strategy of every
replicate, as records of `int32 irep, int32 count, int32 atom, double p[16]` (140 bytes).

**`neutral`** draws M strategies from the mutation law and classifies each: the α-measure of the atoms at that game
(default seed 777000000 + task). One line:

    NEU gid ipt u v Emax alpha M c000 c001 c010 c011 c100 c101 c110 c111 nnashscan nrivscan t

**`binall`** classifies all 65536 binary genomes at the game and writes `B<task>.bin` (one byte per genome, the atom
code) and a line `BIN gid ipt c000 … c111 t`. It is the check against the exact binary atom table.

## Task scripts, task files and Slurm scripts

`task.sh` reads line `$SLURM_ARRAY_TASK_ID` of `$TASKFILE` (comment lines skipped):
`alpha game rep0 nrep neutralM seedbase`. It runs `wfc run` and writes `$OUTDIR/r_<a>_<game>_<rep0>.txt` and
`$OUTDIR/S_<a>_<game>_<rep0>.bin`, where `<a>` is α with `/` replaced by `o` (`1o3`) or the permille value. If
neutralM > 0 it also writes `$OUTDIR/n_<a>_<game>.txt`. Environment: `ITEND`, and optionally `WFC` (default `./wfc`),
`OUTDIR` (default `out`), `NT` and `NI` (default 1000 and 4).

`task_neu.sh` reads `alpha game M seed chunk` and writes `out/nr_<a>_<game>_<chunk>.txt` (`wfc neutral`).

| Slurm script | task file | array | binary, environment | content |
|---|---|---|---|---|
| `probe.sbatch` | `tasks_probe.txt` | 1–88 | `wfc`, ITEND = 10⁷ | α = 1/2 at 8 games × 10 replicates; α = 0.333 at 8 games × 1 replicate, neutral 20 000 |
| `probe2.sbatch` | `tasks_probe2.txt` | 1–28 | `wfc`, ITEND = 10⁷ | α = 1/2 at games 30 and 38 (wedge W) × 10; α = 1/3 at 8 games × 1, neutral 20 000 |
| `run_a333.sbatch` | `tasks_a333.txt` | 1–512 | `wfc`, ITEND = 10⁷ | F1: α = 1/3, 512 games × 10 replicates, neutral 20 000 per game, seedbase 9000000 |
| `run_neu.sbatch` | `tasks_neu.txt` | 1–84 | `wfc` | 21 games × 4 chunks × 250 000 neutral draws (10⁶ per game) |
| `probe_f2.sbatch` | `tasks_probe_f2.txt` | 1–48 | `wfc2`, ITEND = 10⁸, OUTDIR = `out_f2`, NT = 4000, NI = 1 | α = 1/2 at games 30, 165, 172, 184 × 10; α = 1/3 at 8 games × 1 |
| `run_f2.sbatch` | `tasks_f2.txt` | 1–512 | `wfc2`, ITEND = 10⁸, OUTDIR = `out_f2`, NT = 4000, NI = 1 | F2: α = 1/3, 512 games × 10 replicates, seedbase 9100000, no neutral draws |

All arrays ran on the `sapphire` partition, one CPU and 2 GB per task. F2 uses the α-measures of F1 (the same α and the
same games).

## Arrays

Every array was checked per task with `sacct`; the submission ledger lines are in
`PRFull/release/accounting/jobs_tsv_lines.tsv` and the `sacct` records in `PRFull/release/accounting/sacct_table7.psv`.

| JOBID | file | tasks | what | state |
|---|---|---|---|---|
| 48329275 | `JOBID_probe` | 88 | first probe: α = 1/2 validation, α = 0.333 cost | 81 COMPLETED; tasks 82–88 FAILED at the `mv` of a strategy dump whose name lacked α (a name collision in the clean-up step); their REP lines are complete |
| 48343431 | `JOBID_probe2` | 28 | second probe: exact 111 count at α = 1/2, cost at α = 1/3 | 28/28 COMPLETED |
| 48366947 | `JOBID_a333` | 512 | run F1 | 512/512 COMPLETED, 5120 REP lines, empty stderr, 1:08–1:57 h per task |
| 48391602 | `JOBID_neu` | 84 | 10⁶ further neutral draws at 21 games | 84/84 COMPLETED |
| 48723624 | `JOBID_probe_f2` | 48 | F2 probe: α = 1/2 validation, α = 1/3 cost (444–514 s per replicate, classification included) | 48/48 COMPLETED, empty stderr |
| 48724650 | `JOBID_f2` | 512 | run F2 | 512/512 COMPLETED, 5120 REP lines, empty stderr, 1:00–1:49 h per task, 692 core-hours |

The 21 games of 48391602 are those where 20 000 draws left the most enriched atom undetermined at 95 %.

Changes made after the first probe and before the second: `u01()` returns 52-bit values strictly inside
(0, 1) (it could return 1.0, so an index could go out of bounds); α is read as an exact fraction (it was permille
only); every output name carries α; atom 111 in wedge W is counted exactly instead of sampled.

## Validation

| check | result |
|---|---|
| `binall` at 9 games covering all 4 wedges, against the binary atom table `PartnersRivals/Cannon/atoms512.bin` (`../../check/wfc_binall.py <dir>`) | 0 disagreements in 589 824 genome-games |
| α = 1/2 against the binary run F1, 8 games × 10 replicates (reference `PartnersRivals/Data/f1_e4.rep`), array 48329275 (`../../check/wfc_probe.py <dir>`) | efficiency: Welch \|t\| ≤ 1.68 at all 8 games; mean difference ≤ 3.2e-4 at 7; at the 8th 0.9845 here, 0.9788 in F1 (one F1 replicate at 0.928). Atom shares: largest \|Welch t\| 2.27 of 64. Most abundant atom the same at 8/8 |
| exact 111 count, α = 1/2, wedge-W games 2029 and 2037, array 48343431, against F1's exact histogram | 3.73e-5 here, 3.18e-5 in F1 (t = +1.0); 4.63e-5 here, 3.89e-5 in F1 (t = +1.2) |
| replicate 1 of 48366947 against 48343431, same seeds, 8 games | bit-identical, 8/8 |
| a replicate run alone against the same replicate inside a multi-replicate task | identical |
| `wfc2`, α = 1/2, against binary run F2 (`PartnersRivals/Data/f2_e4.rep`), games 2029, 2164, 2171, 2183 × 10 replicates, array 48723624 (`../../check/wfc_probe.py <dir> F2`) | most abundant atom the same at 4/4; largest \|Welch t\| 2.38 of 27 non-degenerate statistics (3 above 2, about 1.7 expected); efficiency here/F2: 0.9342/0.9319, 0.9789/0.9750, 0.9985/0.9985, 0.9571/0.9474 |
| replicate 1 of 48724650 against 48723624 at the 8 probe games | bit-identical apart from the timing fields |

## Where the output is

On Cannon the kit ran in `~/PartnersRivalsM2f_wfc` (F1 and the probes write to `out/`, F2 to `out_f2/`). Copies of the
per-replicate output:

* `../../check/wfc_a333/`: run F1 (`r_1o3_<task>_1.txt`, `n_1o3_<task>.txt`) and the neutral refinement
  (`nr_1o3_<task>_<chunk>.txt`);
* `../../check/wfc_f2/`: run F2 (`r_1o3_<task>_1.txt`);
* `../../check/wfc_probe48329275/`, `../../check/wfc_probe48343431/`, `../../check/wfc_probe48723624/`: the probes.

The raw output with the strategy dumps and the Slurm logs is archived in `../archive/`:
`cannon_2026-09-26.tgz` (7208 files: `out/`, `out_probe1/`, `out_probe2/`, `out_f2_probe/` and two earlier kits of this
project) and `cannon_out_f2_2026-09-26.tgz` (`out_f2/` of run F2).

## Analysis

    bash compile.sh
    jid=$(sbatch --parsable run_a333.sbatch); echo $jid > JOBID_a333
    sacct -j $(cat JOBID_a333) --format=JobID,State,ExitCode -X --noheader | awk '$2 != "COMPLETED"'   # must print nothing
    python3 ../../Figures/figwfc.py <dir with r_1o3_* n_1o3_* nr_1o3_*>                                   # Fig4c, Fig5c
    python3 ../../Figures/figwfc.py ../../check/wfc_f2 --run F2 --neutral ../../check/wfc_a333            # SIFig4c, SIFig5c
    python3 ../../check/wfc_stats.py F1        # or F2

* Atom share: fraction of the samples in the atom (in wedge W, atom 111 from `ex111/total`). Enrichment: share divided
  by the atom's α-measure at the game, from `wfc neutral` (20 000 draws per game, plus the 10⁶ of array 48391602 at 21
  games).
* The α-measure of atom 111 depends only on the wedge. `figwfc.py` takes the friendly-rival rate per wedge measured at
  α = 0.333 and scales it to α = 1/3 by the ratio of the stratum census (`../../check/st_eff.txt`);
  `PRFull/Figures/fig45.py` replaces these with the rates measured at α = 1/3 in `PRFull/cannon/sweep13`.
* Error of a share: samples within a replicate are correlated, so the error is taken across replicates.
  `figwfc.py` (`undetermined_rep`) calls the most enriched atom undetermined at 95 % from a Student t interval
  (9 degrees of freedom) over the 10 replicate shares and a Poisson interval on the α-measure; such games are drawn
  dotted.
* `../../check/wfc_stats.py` recomputes from the data the statistics of the paper's main-text paragraph "What
  evolution selects". The paper's SI section on the evolutionary maps describes the runs.
