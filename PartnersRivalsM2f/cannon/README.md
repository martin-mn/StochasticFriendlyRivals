# Cannon kit: the memory-two continuum census

This directory holds the exact ε → 0 engine for memory-two strategies in [0,1]¹⁶ and the Slurm kit that ran the
stratum census, the predicate-free sweep and the dimension probe on Harvard's Cannon cluster. The remote copy of the
kit was `~/CL1/PartnersRivalsM2f/ca1/`. The paper's Computational Supplement describes these campaigns (section "The
census, its campaigns and the rarity sweep") and lists their cost (section "Cost", Table 2).

## Contents

| path | what |
|---|---|
| `m2f.c` | the engine: limit outcome frequencies of a pair of strategies, and the rival margin |
| `witnesses.h` | the 157 binary co-players of the screen (`NWIT`, `WITG`) |
| `census.c` | the census (mode `strata`) and the predicate-free sweep (mode `sweep`) |
| `dimprobe.c` | the maximum-dimension probe |
| `compile.sh`, `task.sh` | build of `census`; one array element of `census strata` or `census sweep` |
| `census_probe.sbatch`, `census.sbatch`, `sweep.sbatch`, `dim.sbatch` | the four Slurm arrays (below) |
| `JOBID_PROBE`, `JOBID`, `JOBID_SWEEP`, `JOBID_DIM` | their array job IDs |
| `m2f_main.c` | command-line driver on binary genomes: `m2f self`, `m2f rival`, `m2f pair` |
| `m2f_py.c` | C entry points for `../check/m2f.py` (ctypes); built with `m2f.c` into `libm2f.dylib` |
| `confirm.c` | full 65 536-co-player confirmation for a list of strata |
| `strata.c` | the self-play class of every stratum, from two draws each |
| `strata2.c` | the class plus a screen against the first 19 witnesses only |
| `onestratum.c` | the per-stratum logic of `dimprobe` for one stratum, with verbose output |
| `wit2.c` | the arg-min co-player and the margin, for both signs, over random strategies |
| `census`, `confirm`, `dimprobe`, `m2f`, `onestratum`, `strata`, `strata2`, `wit2` | macOS arm64 binaries of the above |
| `libm2f.dylib` | macOS universal (arm64, x86_64) shared library |
| `alpha/` | the α-stratified sweep kit (below) |
| `wfc/` | the memory-two evolutionary kit, with its own `README.md` |
| `archive/` | two archives of Cannon output (below) |
| `ca1_out/`, `ca1_out.md5` | release only: the raw output copied from Cannon, and its md5 list (below) |

## Build

On Cannon: `./compile.sh` (it runs `cc -O3 -march=native -o census census.c m2f.c -lm`). Its note: `-march=native`
is safe because every partition used here (`sapphire`, `shared`, `test`) compiles and runs on the same node. The
other programs are one source file each, linked with `m2f.c` in the same way, for example

```bash
cc -O3 -march=native -o dimprobe dimprobe.c m2f.c -lm
cc -O3 -march=native -o m2f m2f_main.c m2f.c -lm
cc -O3 -arch arm64 -arch x86_64 -fPIC -shared -o libm2f.dylib m2f.c m2f_py.c -lm   # macOS, for ../check/m2f.py
```

## Conventions

* A strategy is p ∈ [0,1]¹⁶; p[j] is the probability to cooperate in state j.
* State j = 4·(last outcome) + (outcome before). An outcome is o = 2·(own defection) + (co-player's defection), own
  action first: CC = 0, CD = 1, DC = 2, DD = 3.
* A binary strategy is a 16-bit genome; bit j is 1 for C. ALLC = 65535, ALLD = 0.
* Errors: both players play p as ε + (1 − 2ε)p, independently. The engine computes the limit ε → 0 exactly.
* Payoffs are normalised as (R,S,T,P) = (1,u,1+v,0). The switch line is u + v = 1. T > S means u − v < 1 (sign +1);
  T < S means u − v > 1 (sign −1).
* Self-play class, from the limit frequencies W = (W_CC, W_CD, W_DC, W_DD) with tolerance 1e-12: `C` if W_CC = 1
  (efficient below the switch line), `A` if W_CD + W_DC = 1 (alternation, efficient above it), `Z` if W_DD = 0 only
  (efficient on it), `-` otherwise.
* Rival margin: the minimum over co-players of W_DC − W_CD (sign +1) or of W_CD − W_DC (sign −1). A strategy is a
  rival for that sign when the margin is ≥ −1e-9.
* Wedges: W = class C and sign +1, S = C and −1, N = A and +1, E = A and −1. The two rays of the switch line count
  class C, A or Z (W_DD = 0) with sign +1 or −1.
* Stratum index g ∈ [0, 3¹⁶): base-3 digits d_j, j = 0 least significant. d_j = 0 means p_j = 0, d_j = 2 means
  p_j = 1, d_j = 1 means p_j is free. The census and the dimension probe draw free coordinates uniformly on
  (0.02, 0.98), from a random stream seeded by g, so every stratum's draws are reproducible.

## The engine, `m2f.c`

It builds the 16-state chain of a focal strategy p against a co-player q and reduces it by GTH state elimination.
Each entry is held as a leading coefficient and a power of ε: a probability x > 0 gives (x, 0) for C, and x = 0 gives
(1, 1); D likewise with 1 − x. Only positive quantities are added, multiplied and divided, so the leading terms are
exact, with no cut-off.

* `freqs(p, q, W)`: the four limit outcome frequencies from the focal side; `freqs16` gives all 16 states.
* `rival_margin(p, sign, &arg, earlyquit)`: the minimum of sign·(W_DC − W_CD) over the 65 536 binary co-players, with
  the arg-min genome in `arg`. It stops once the minimum falls below `earlyquit`. Why the binary co-players suffice
  is in the paper's SI section "Co-players and deviators of any memory".

`m2f.c` generalises the binary engine `PartnersRivals/check/pairs.c` of the earlier binary-strategy project (not in
this release) to the continuum. Its validation is in `../check/README.md`.

## Programs and their output

`census strata <stripe> <nstripes> <Keff> <st_out> <pt_out> [Kother]` visits the strata g ≡ stripe (mod nstripes).
It classifies each stratum by the self-play of one draw. It then takes Keff draws if the class is C, A or Z, and
Kother (default 8) if it is `-`. For each draw and sign, the draw is first screened against the 157 witnesses; a draw
with a negative margin against any of them is not a rival. Survivors go through `rival_margin` over all 65 536.

* `st_out`: one line `g cls nP nM K` per stratum with at least one confirmed rival draw (nP, nM: confirmed draws for
  sign +1 and −1, out of K).
* `pt_out`: lines `P g cls p0 … p15` and `M g cls p0 … p15` at full double precision, at most 4 per stratum and sign,
  and only for strata of class C, A or Z.
* stderr: one summary line per task (strata emitted, number of full scans).

`census sweep <n> <seed> <out>` draws n strategies: half uniform on the cube, half with each coordinate 0 with
probability 0.36, 1 with probability 0.36, else uniform. It classifies each one and tests both signs, with no use of
any closed form.

* `out`: one line `W p0 … p15` for each draw in wedge W.
* stderr: one line `sweep n=…: effCC effALT effZ rivP rivM | FR_W FR_S FR_N FR_E | rayP rayM` with the counts.

`dimprobe <dmin> <K> <sign> <class> <stripe> <nstripes>` visits the strata with at least dmin free coordinates whose
first draw has the given class. It takes up to K draws per stratum and stops at the first confirmed rival. stdout:
one line `g nfree cls` per stratum with a rival; stderr: the count and the largest number of free coordinates.

The small tools:

* `confirm <stratumlist> <sign> <K> <stripe> <nstripes>`: for each listed stratum, K draws on (0.07, 0.93), each
  scanned against all 65 536; stdout `g r W_CC W_CD W_DC W_DD p0 … p15` per confirmed rival.
* `strata <stripe> <nstripes> [out]`: lines `C g`, `A g`, `Z g`; two draws per stratum on (0.12, 0.88), and any
  disagreement between them goes to stderr.
* `strata2 <stripe> <nstripes> [K=4] [out]`: lines `g cls sP sM`, the draws that pass the 19-witness screen.
* `onestratum <g> <K> <sign>`: the stratum's digits, its class and how many of K draws pass the screen and are
  confirmed.
* `wit2 <n> <seed>`: stdout `argP marginP argM marginM` for n random strategies.
* `m2f self <stripe> <nstripes>`: `g W_CC W_CD W_DC W_DD` for binary self-play; `m2f rival <stripe> <nstripes>
  <sign>`: `g margin argmin`; `m2f pair <g1> <g2>`: the four frequencies of one binary pair.

## What was submitted

| script | array | what |
|---|---|---|
| `census_probe.sbatch` | 4 tasks, `test`, 2 h | `census strata` on stripes 0–3 of 20 000, KEFF = 64, KOTH = 8: the per-unit rate |
| `census.sbatch` | 512 tasks, `sapphire`, 12 h | `census strata` on all 3¹⁶ = 43 046 721 strata, 512 stripes, KEFF = 64, KOTH = 8; writes `out/st_<i>.txt`, `out/pt_<i>.txt` |
| `sweep.sbatch` | 256 tasks, `sapphire`, 12 h | `census sweep` with n = 4 × 10⁶ and seed = task + 1000: 1.024 × 10⁹ strategies in all; writes `out/sw_<i>.txt` |
| `dim.sbatch` | 384 tasks, `sapphire`, 12 h | `dimprobe 5 512 <sign> <class> <stripe> 64` for the six (class, sign) pairs (W, S, N, E and the two rays) × 64 stripes; task i runs pair i / 64 on stripe i mod 64; writes `out/dim_<class>_<sign>_<stripe>.txt` |

`task.sh` takes the mode (`strata` or `sweep`); the array index picks the stripe, and `NSTRIPES`, `KEFF`, `KOTH`,
`NSWEEP` come from the batch script. Slurm logs go to `out/<prefix>_<jobid>_<task>.out` and `.err`.

Sizing used the probe's own rate: 210–270 s per 1/20 000 of the strata on Cannon (184 s on the local machine). That
gives about 1330 core-hours for the census, about 2.6 h for each of 512 tasks.

## Job IDs

Each ID was written to its file at submission and appended to the cross-project ledger `~/jobs.tsv` on Cannon in the
same step. The release's `accounting/jobs_tsv_lines.tsv` has the ledger lines.

| file | job | tasks |
|---|---|---|
| `JOBID_PROBE` | 45749304 | 4, the census probe |
| `JOBID` | 45751662 | 512, the census |
| `JOBID_SWEEP` | 45751663 | 256, the sweep |
| `JOBID_DIM` | 45752239 | 384, the dimension probe |

All tasks of the four arrays COMPLETED; the per-task `sacct` record is in the release's
`accounting/sacct_table7.psv`. To check a campaign on Cannon:

```bash
squeue -j $(cat JOBID)
sacct  -j $(cat JOBID) --format=JobID,State,ExitCode -X --noheader | awk '$2 != "COMPLETED"'   # empty = clean
```

An empty `squeue` only means nothing is running; completion comes from `sacct`, per task.

## Collection and receipts

* `../check/collect.sh` refuses to go on unless `sacct` shows all 512 census tasks COMPLETED. It then copies the
  census lines of class C, A or Z to `../check/st_eff.txt`, runs `../check/census.py` on them (output
  `../check/census.out`), and sums the sweep's stderr lines into `../check/sweep.out`.
* `../check/dim.out` summarises the dimension probe per (class, sign) pair.
* The stratum output of all 512 census tasks, as Cannon wrote it, is in the release's
  `PRFull/release/census_class/census_st_45751662.tgz`, with the md5 of each file.
* `ca1_out/` (release only) holds the rest of the raw output, copied read-only from Cannon's `ca1/out`:
  `pt_<i>.txt` (512); `sw_<i>.txt` (256) with the Slurm logs `sw_45751663_<i>.out`/`.err`;
  `dim_<class>_<sign>_<stripe>.txt` (384) with the logs `dim_45752239_<i>.out`/`.err`; and the census probe logs
  `probe_45749304_<i>.out`/`.err`; 2440 files in all. `ca1_out.md5` is the md5 of each file taken on Cannon, and every
  copied file matches it.

## The witness set

`witnesses.h` lists 157 binary co-players. The first 19 are a greedy set cover of the 62 896 binary non-rivals: each
binary strategy that is not a rival has a negative margin against one of them. The other 138 are the co-players that
most often gave the minimum margin in a sample of 22 125 random strategies from the continuum. A draw that fails
against any witness is certainly not a rival, so the screen never drops a rival; every draw that passes is still
scanned against all 65 536.

## `alpha/`: the α-stratified sweep

Remote kit `~/PartnersRivalsM2f_alpha`. Same `m2f.c` and `witnesses.h` as here. Its `census.c` keeps mode `strata`
and replaces `sweep` with `psweep <n> <seed> <alpha_permille> <out>`: each coordinate is 0 with probability α, 1 with
probability α, and uniform on (0,1) otherwise.

* `compile.sh` builds `census`; `task.sh` runs `census psweep $NSWEEP <seed> $AM out/pt_<AM>_<i>.txt`, with seed
  100000 + 1000·AM + task index, and writes stdout to `out/a_<AM>_<i>.tsv` and stderr to `out/a_<AM>_<i>.err`.
* stdout is one tab-separated line `ALPHA am n effCC effALT effZ rivP rivM FR_W FR_S FR_N FR_E rayP rayM`; stderr
  repeats the counts with the mean fraction of pinned coordinates; the point file has `W p0 … p15` lines.
* `probe.sbatch`: 8 tasks, n = 10⁵ each (job 45830542, `JOBID_probe`; its α grid is in its ledger line).
  `run_333.sbatch`: 512 tasks at α = 0.333, n = 3 190 000 each (job 45831015, `JOBID_a333`). `run_250.sbatch` and
  `run_450.sbatch` are also present; no job ID is recorded for them. All tasks of 45830542 and 45831015 COMPLETED.
* `collect.sh <JOBIDfile>` refuses unless every task COMPLETED (`sacct`), then sums the `ALPHA` lines per α.
* `../check/sweep_a333.tsv` holds the 512 `ALPHA` lines of job 45831015. The raw output is in
  `archive/cannon_2026-09-26.tgz` under `PartnersRivalsM2f_alpha/`.

## `archive/`

* `cannon_2026-09-26.tgz`, made on 2026-09-26: the kit `~/PartnersRivalsM2f_wfc` with its finished outputs `out/`,
  `out_probe1/`, `out_probe2/` and `out_f2_probe/`, and the whole of `~/PartnersRivalsM2f_alpha` and
  `~/PartnersRivalsM2f_task2`, outputs included.
* `cannon_out_f2_2026-09-26.tgz`: `out_f2/`, the output of the N = 100 evolutionary run (see `wfc/README.md`).
