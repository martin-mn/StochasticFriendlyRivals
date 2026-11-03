# `PRFull/cannon/cw2` — internal review campaign, memory-one cube: error rates E3, E5 and the time series TS

Internal review, round 1, items 7 and 8 (internal review notes, not deposited: report 6 and §1f, 2e).
It asks two things of the memory-one Wright–Fisher run over the cube [0,1]^4 (PartnersRivalsM1f `cannon/cw`, run F1):

* **E3, E5.** Does the run change if the error rate is 1e-3 or 1e-5 instead of 1e-4? These are 10 replicates × 1e7 generations
  at 71 games, in the F1 regime (N = 1000, β = 100, u = 1e-2).
* **TS.** Is the run converged? Time series at 4× production length (4e7 generations) at 24 games × 10 replicates, one replicate
  per task. For every block of 2e4 generations over the whole run, the task records the exact atom counts, the efficiency, the
  number of live lineages and the number of mutants.

**Status (2026-09-27): built and validated locally. NOTHING HAS BEEN COPIED TO OR SUBMITTED ON CANNON.**

## 1. Source: `cwx.f` = production `cw.f` + three switches

`ref/` holds byte copies of the production sources (sha256 in `ref/SHA256SUMS`). `cw.f` is `37ec579b…c729a2ebb`, the F1
production source that mkcw2.py asserts. Also in `ref/`:
- `prod/`: the archived production outputs of isl 3 and 1401 (run1, job 48366968);
- `probe1e6/`: the Cannon probe output of isl 3 (job 48340155);
- `packed/`: `cw1_e4.tsv`, `cw1_e4.rep.tsv`, `kappa_e4.tsv`, the comparator for the analysis.

`compile.sh` refuses to build if a copied file differs from `ref/`. `mt.f mtb.f payf2.f payf3.f m1cls.c m1cls.h games.dat` are
byte copies of production's. Only `cwx.f` and `kappa.c` are edited, and **`diff ref/cw.f cwx.f` is the audit**: 19 hunks, 209
added lines (128 of them code, the rest comments), 3 changed lines. Every hunk is tagged `[cw2 n]`.

| tag | switch | what it changes | default |
|---|---|---|---|
| [cw2 1] | `-DEPS1=<v>` (compile) | The value of eps slot ie = 1 in `data epsv /EPS1,1.d-03,1.d-04/` (changed line 1). The E5 arm builds `-DEPS1=1.0d-5` and runs at isl = 3·ipt+1. Slots 2 (1e-3) and 3 (1e-4) are untouched. `kappa.c` gets the same macro (2 hunks). | `1.d-02`, production's |
| [cw2 2] | `cwx.x <isl> [irep0 irep1 [seedbase]]` (runtime) | Runs only replicates irep0..irep1, with seeds seedbase + (isl−1)·10 + irep (changed lines 2 and 3: `do irep=irep0,irep1`, `iseed=jseedb+…`). With explicit arguments, the h-file gets one extra `# reps:` comment line. A replicate run alone is the same trajectory as inside the full task: the pair cache is keyed by never-reused serials, so it only decides whether a payoff is recomputed. Only the w-file's cumulative `gcount` differs. | 1..10, base 40000000: the production call `cw.x <isl>` |
| [cw2 3] | `-DTSERIES [-DNBLK=k]` (compile) | Classifies every lineage **at birth** (atom only, through the same `m1cls` call, `tscls`). Over the whole run, in NBLK blocks, it accumulates the exact individual-generations of each atom from the live lineages, taken at the point where cw.f samples, plus the summed payoff, the live-lineage count and the mutants. Output `t<isl>`: `TS irun iblk it_last eff meanlive nmut c0..c7`. No random number is drawn and no production variable is written. **Asserted every replicate:** the second-half blocks sum to cw.f's own `atomcnt`, atom by atom, and every block holds exactly lblk·N individuals. TSERIES with BINMODE is an `#error`. | off; NBLK = 2000 |

Builds (`bash compile.sh`, then `bash compile.sh probe`):

| binary | source + defines | used by |
|---|---|---|
| `cw.x` | `cwx.f` (no switch) = the production program | E3 (isl = 3·ipt+2), VAL |
| `cw_e5.x` | `-DEPS1=1.0d-5` | E5 (isl = 3·ipt+1) |
| `cw_ts.x` | `-DTSERIES -DITEND=40000000` (lblk = 2e4) | TS, TSX |
| `kappa.x`, `kappa_e5.x` | `kappa.c`, `-DEPS1=1e-5` | KA (probe) |
| `cwp.x`, `cwp_e5.x`, `cwp_ts.x`, `cw_ts1.x` | itend 1e6, 1e6, 4e6; TSERIES at 1e7 | probe only |

The compilers are production's: `module load intel/25.2.1-fasrc01`, then `ifx -O3 -xHost -fp-model precise -fpp` and
`icx -O2 -xHost -fp-model precise -ffp-contract=off`. Run on sapphire or test only. Locally, gfortran 14 and clang use
`-O2 -ffp-contract=off`.

## 2. Validation (local, 2026-09-27, gfortran 14.1 / Apple clang 15; at most 2 cores, nice 10)

| check | result |
|---|---|
| **(i) default build vs Cannon production, full length.** `cw.x 3 1 1` (replicate 1 alone, production seed) vs `ref/prod/3` line 1 and `ref/prod/w3` line 3 (run1, 48366968, built with ifx on Cannon) | **REP line and w-line bit-identical**; 306.7 s locally |
| default build vs Cannon probe (48340155, itend 1e6): `cwp.x 3 1 1` vs `ref/probe1e6/3` line 1 and the w-line | **bit-identical** (gfortran locally = ifx on Cannon) |
| default build (and the pristine `ref/cw.f`) at itend 4000 vs the archived local chain `PartnersRivalsM1f/cannon/cw/run_local` (gfortran), 8 cells (ipt 0 29 186 246 466 480 500 415), all 10 replicates | `<isl>`, `w<isl>`, `h<isl>`: **24 of 24 files identical** |
| **(ii) TSERIES**, full length: `cw_ts1.x 3 1 1` | REP line and w-line **bit-identical to production**. All three files identical to the default build. In-program block check passed (2000 blocks). Cost +9.1 % (334.6 s vs 306.7 s) |
| TSERIES, itend 4000 (2000 blocks of 2), 6 cells at all three eps slots, 10 replicates | `<isl>`, `w<isl>`, `h<isl>` identical to the default build at 6 of 6 cells. Also identical with the TS arm's call (single replicate, seed base 43000000) at 2 cells |
| replicate range: replicate k alone (k = 2, 7, 10) vs line k of the full task, 3 cells | REP line identical 9 of 9. The w-line differs only in the cumulative `gcount` |
| seed base: `1 10 40000000` given explicitly vs the production call | `<isl>` and w-file identical. The h-file differs only by the `# reps:` line. Base 43000000 gives a different trajectory |
| EPS1: `-DEPS1=1.0d-5` vs the default build at eps slots 2 and 3 (6 cells) | identical at 6 of 6. At slot 1 the h-file says eps 0.100E-04 (default 0.100E-01) |
| kappa: `kappa.x` vs `ref/kappa.c` at 6 isl, and `kappa_e5.x` vs `kappa.x` at slots 2 and 3 | identical |
| kernel masses at eps 1e-3 and 1e-5, 71 games, M = 1e7 (local run of the probe's KA rows) vs production kappa_e4 | max \|z\| = 2.96 over 1136 comparisons: **eps-independent within sampling**, as exact arithmetic says (p = (s−ε)/(1−2ε) has the same law at every ε) |
| **(iii) every arm locally**: `tasks_local.txt` through `task.sh` (test builds at itend 2e4: E3, E5, two TS replicates, KA at 24 games), then `collect.py --local` | all complete, and idempotent (a rerun skips). A wrong-eps binary is refused (exit 3). A truncated h-file is detected and rerun, and the rerun reproduced the file |
| the real arm binaries, briefly | `cwp.x` E3 at ipt 0, 1e6: E = 0.983 (production at eps 1e-4: 0.985). `cwp_e5.x` at ipt 466: E = 0.974. `cwp_ts.x` at ipt 0, 4e6: 135.6 s, 2000 TS lines |
| dry runs with fake `sbatch`/`squeue`/`sacct` (scratch copy, HOME in scratch) | submit.sh: refuses an arm before PROBE_OK, refuses a second probe, refuses TS (240 > 100) without QUEUE_OK. It writes JOBID_\<arm\> and ledger lines of the form `date PRFull cw2 jobid ntasks [arm] note`. verify.sh probe runs the VAL/VALTS comparison, KA and rates |

**Worth knowing before E3 is read.** At ipt 466 (u = −8.75, v = 8.42, deep PD), a 1e6-generation E3 replicate stays at
E = 0.005, against 0.974 at eps 1e-5 and 0.973 in production. This is the effect the binary paper [PR2026] reports at eps = 1e-3:
deep in the Prisoner's Dilemma every error costs the exploited partner R − S. It is physics, not a switch: E3 runs the unchanged production code path at eps slot 2. The full-length E3 replicate
at 466 is in the probe.

## 3. Games

* **`games71.txt`, the E3/E5 subset.** It is the 64 games plus the special games of the internal review notes. The 64 are given as
  **games.dat line numbers**, i.e. wfc task numbers, as `wfc/README.md` says. The special games 466 411 301 237 326 415 500 246
  337 are **ipt**, the 0-based p### of the tags (for example ipt 466 = V_p466_PDW, (u,v) = (−8.754, 8.416)). ipt 301 is
  line 302 and ipt 237 is line 238, and both are already among the 64. **So the union is 71 games, not 72.** Per wedge: W 19,
  S 16, N 18, E 18.
* **`games24.txt`, the TS subset**, 6 per wedge, following report 6 of the internal review notes. It is **the memory-two kit's list**
  (`../wfc2/mk_tasks.py` `TS_IPT`), adopted so that both memories are diagnosed at the same games:
  - W: 466 411 301 plus the typical PDW 186 377 29;
  - S: the split SHS 184 417 451 plus the typical 208 17 171;
  - N: the split PDN 500 390 445 479, plus 337 (special, split SDN) and 256 (typical SDN);
  - E: the special split HAE 237 326 415 plus the split SDE 109 164 117.

  "Split" means the ten replicates of memory-two F1 do not all name the pooled leader. "Typical" means unsplit, among the 64,
  with efficiency nearest its wedge-and-quadrant median. Six of the 24 are outside the 71. **If either kit changes the list, both
  must change it** (`python3 mktasks.py` here).
* **`games_tsx.txt`, OPTIONAL arm TSX**: 246 75 311 185. Of the 24 TS games only 256 is split at memory one, and only barely
  (000:1/100:9). The memory-one splits are Snowdrift games of N and E, 000 against 100. The internal review notes (§1f) say that this patchwork needs replicate diagnostics. TSX
  adds them for ~14 core-h.

## 4. Seeds

Replicate irep of task isl uses seed = base + (isl−1)·10 + irep, with MT19937 `init_genrand` (cw.f).

| run | isl | base | seeds |
|---|---|---|---|
| production F1 | 3·ipt+3 | 40000000 | 40000000 + (3·ipt+2)·10 + irep |
| **E3** | 3·ipt+2 | 40000000 | 40000000 + (3·ipt+1)·10 + irep — other isl, disjoint |
| **E5** | 3·ipt+1 | 40000000 | 40000000 + (3·ipt)·10 + irep — disjoint |
| **TS, TSX** | 3·ipt+3 | 43000000 | 43000000 + (3·ipt+2)·10 + irep — disjoint (TS and TSX share no game) |
| probe VAL, VALTS | 3 | 40000000 | = production isl 3 replicate 1, **on purpose** |

The bases 41M (kappa: 41000000 + isl) and 42M (F2), and 6M–39M (other campaigns), are not touched. The probe rows PE3, PE5 and
PTS run replicate 1 of one task of each arm with that arm's own seed. `verify.sh <arm>` checks that the arm reproduces them.

## 5. Arms, tasks, cost

The rates are measured on sapphire in the production array 48366968, i.e. the same program and the same hardware. At the 71
games a production task (10 × 1e7 generations) took 2403–3094 s, 2836 s on average. At the 24 TS games it took 2576–3094 s.
TSERIES adds 9.1 % (measured locally).

| arm | sbatch | tasks | per task | --time | core-h |
|---|---|---|---|---|---|
| probe | `probe.sbatch` (test) | 13 | 30 s – 23 min | 1:00:00 | ~1.1 |
| E3 | `e3.sbatch` | 71 = 71 games × 10 replicates | ~47 min (≤ 52) | 1:45:00 | **55.9** |
| E5 | `e5.sbatch` | 71 | ~47 min | 1:45:00 | **55.9** |
| TS | `ts.sbatch` | 240 = 24 × 10, one replicate each | ~21 min (≤ 23) | 0:50:00 | **84.6** |
| TSX (optional) | `tsx.sbatch` | 40 = 4 × 10 | ~20 min | 0:50:00 | 13.6 |
| **total** (without TSX) | | 395 | | | **~197** |

- Every arm uses `--partition=sapphire` (probe: test), `--job-name=PRFull.cw2.<arm>`, `--comment=PRFull`, `--cpus-per-task=1`
  and `--mem=400M` (production MaxRSS was 34 MB).
- Logs go to `out/log/<arm>_%A_%a.{out,err}`.
- Each `--time` is 2× the worst production task.
- `verify.sh probe` prints the worst measured probe rate per arm and says `RAISE --time` if 2× that exceeds the script's value.
  eps may change the cost: fewer live lineages in a defecting population are cheaper, a richer population is dearer.
- Disk on Cannon:
  - E3/E5: about 7 KB per task, plus one 72 KB copy of games.dat per arm directory.
  - TS: about 390 KB per task (the 316 KB t-file plus the task directory's games.dat), so **~95 MB** for the arm.
  - Everything is small: raw per-task output only, no dumps.

## 6. Protocol (Cannon; Martin opens the session; the main session does the cluster work)

Copy the kit without Mac binaries or outputs. Suggested remote path: `~/CL1/PRFull/cannon/cw2`.

    rsync -a --exclude build --exclude '*.x' --exclude out --exclude 'results*' --exclude validation \
          /Users/martin/Documents/CL1/PRFull/cannon/cw2/ cannon:CL1/PRFull/cannon/cw2/

On a login node, in the kit:

    bash compile.sh                  # cw.x cw_e5.x cw_ts.x kappa.x kappa_e5.x   (ifx/icx; refuses if ref/ sha256 fails)
    bash compile.sh probe            # cwp.x cwp_e5.x cwp_ts.x cw_ts1.x
    bash submit.sh probe             # = sbatch --parsable probe.sbatch -> JOBID_probe + ~/jobs.tsv line
    bash verify.sh probe             # when done: sacct per task, stderr, products, VAL/VALTS bit-identity, KA, rates -> PROBE_OK
    #   raise --time in e3/e5/ts.sbatch if verify says so, then
    bash submit.sh E3                # 71 tasks
    bash submit.sh E5                # 71 tasks
    squeue -u mnowak -h -t PENDING,RUNNING -o '%A' | wc -l      # fair share: report the number, ask before TS (240 tasks)
    QUEUE_OK=1 bash submit.sh TS     # 240 tasks
    bash verify.sh E3; bash verify.sh E5; bash verify.sh TS     # until each says "clean"

`submit.sh` is the protocol in one step:
- `jid=$(sbatch --parsable <arm>.sbatch)`, `echo $jid > JOBID_<arm>`;
- one ledger line `date PRFull cw2 jobid ntasks [arm] note` appended to `~/jobs.tsv`;
- it refuses a second submission of an arm with a JOBID file (`RESUBMIT=1` keeps `.prev`);
- it refuses the arms before `PROBE_OK`;
- it refuses more than 100 tasks without `QUEUE_OK=1`;
- it refuses binaries that are not cwx.f builds and task files that differ from `mktasks.py`'s output.

Completion is decided **only** by `verify.sh <arm>`:
- `sacct -j $(cat JOBID_<arm>)` per task, where anything not COMPLETED is listed;
- empty stderr;
- `task.sh --check` on every task line, i.e. complete products naming the right eps and itend;
- the core-hours.

Never `squeue -u` as a basis, never `scancel` without the JOBID. A FAILED task is investigated against its files before any
rerun. `task.sh` is idempotent: a complete task is skipped.

Back on the Mac (only when Martin says the outputs may be fetched):

    bash collect.sh fetch CL1/PRFull/cannon/cw2     # rsync out/ and JOBID_* back, then the tables

## 7. Analysis: `collect.py` (local, numpy/scipy), formats of the paper's census (`Figures/fig45.py`)

- **E3, E5.** Per wedge (W S N E all):
  - efficiency min/median/max;
  - the most abundant and the most enriched atom (games per atom);
  - the median friendly-rival (111) share and the median competitive (001+011+111) share, each with a 95 % bootstrap interval
    over replicates within games;
  - the games whose ten replicates all name the pooled leader.

  The same table for production at the same games, then the comparison: games with the same leader and the same most enriched
  atom, ΔE with Welch t per game, and the success criteria. Per-game TSV: `results/E3_games.tsv`, `E5_games.tsv`. Kernel masses
  (enrichment denominators) are production's `kappa_e4` for the arm AND for production: the kernel's law does not depend on
  eps, and one common denominator keeps Monte-Carlo noise in kappa (the 3/M floor at atoms of a few draws) from moving the
  most enriched atom. The probe's KA masses are the eps-independence check (`--kappa-check`); the table also says at how many
  games the arm's most enriched atom would change with them.
- **TS (and TSX).** Per game and replicate:
  - Q3 against Q4 of the production window (5e6, 1e7], with batch-means SE from 20 batches, for efficiency and for the leading
    atom's share (the production pooled leader);
  - Geyer's initial-positive-sequence IAT and ESS of the leading-atom share over (5e6, 1e7] (250 blocks) and (2e7, 4e7] (1000
    blocks);
  - the run-length check: leader and efficiency of (5e6, 1e7] against (2e7, 4e7] across replicates;
  - comparison of the TS window with production (independent chains);
  - 111 above ½: replicates ever, with a Clopper–Pearson interval, and takeovers and escapes per 1e7 generations, with Garwood
    intervals.

  Outputs: `results/TS_games.tsv`, and `results/TS_series.npz` with all block series for the SI figure (time series, cumulative
  means).
- **Pure-python helpers** used by `verify.sh` on Cannon: `--kappa-check` and `--probe-times`.

## 8. Success criteria (from the internal review notes, stated before the runs)

**E3, E5** (as in the binary paper [PR2026]). Each is printed PASS or FAIL:
1. the most abundant atom is unchanged at ≥ 95 % of the 71 games;
2. efficiency is within replicate spread: |Welch t| ≤ 2 at ≥ 95 % of games and none above 3;
3. friendly rivals (111) still lead all 16 S games;
4. the W/N/E medians of the 111 share and of the competitive share lie inside production's bootstrap intervals.

The deep-PD effect above suggests that E3 may fail (2) there and at a few other games. As in the binary paper, that is itself the result
to report: the atoms are limit properties and do not move, while efficiency deep in the PD does.

**TS.** At every game whose production replicates agree on the leader:
- the quarters agree, i.e. at most 2 of 10 replicates beyond 2 SE, for efficiency and for the leading share;
- the ESS of the leading-atom share over the production window is ≥ 50 in every replicate (the internal review notes' "ESS at least 50 per replicate"; the median is printed too);
- leader and efficiency agree between (5e6, 1e7] and (2e7, 4e7] within replicate spread.

At split and takeover games, report whether the 111 state is absorbing (the takeover and escape counts) and restate the affected
counts as finite-horizon statements.

## 9. Files

| file | role |
|---|---|
| `cwx.f` | the simulator: `ref/cw.f` + [cw2 1–3] |
| `kappa.c` | the kernel-mass program + the EPS1 macro |
| `ref/` | byte copies and archived production outputs for the checks, with sha256 |
| `mt.f mtb.f payf2.f payf3.f m1cls.c m1cls.h games.dat` | byte copies of production's |
| `games71.txt games24.txt games_tsx.txt` | the subsets, with the reasons |
| `mktasks.py` | writes every `tasks_*.txt` (`--check` compares) |
| `tasks_e3.txt tasks_e5.txt tasks_ts.txt tasks_probe.txt tasks_tsx.txt tasks_local.txt` | one line per task: `arm exe isl ipt irep0 irep1 seedbase eps itend outdir` |
| `compile.sh` | the builds (run set, probe set, local test set) |
| `task.sh` | one task, idempotent, with a completeness and build check (`--check`) |
| `probe.sbatch e3.sbatch e5.sbatch ts.sbatch tsx.sbatch` | the Slurm scripts |
| `submit.sh verify.sh` | the JOBID/ledger protocol and the sacct-based verification |
| `collect.py collect.sh` | the analysis |
| `validation/` | the local full-length outputs behind section 2 (default and TSERIES, isl 3 replicate 1; 1e6 probe) |
| `out/local/`, `results_local/` | the local smoke test (itend 2e4): **not a result** |

## 10. Adversarial verification (2026-09-27)

Independent re-run of the local validation in a scratch copy (fresh gfortran 14 / clang builds, at most 2 cores, nice 10):
- `cw.x 3 1 1` (1e7): REP line and w-line **bit-identical** to `ref/prod` (run1, 48366968); 340 s.
- `cw_ts1.x 3 1 1` (1e7): `3`, `w3`, `h3` identical to the default build's; its REP/w lines = production; 2000 blocks of
  5000, each exactly lblk·N; the second-half block atom counts reproduce the REP shares exactly, and the block efficiencies
  average to the REP `ef` within its printed precision. 372 s (+9 %).
- `cwp.x 3 1 1` (1e6) = the Cannon probe archive `ref/probe1e6` (REP and w line).
- itend 4000, 13 cells over all three eps slots, 10 replicates: pristine `ref/cw.f` = `cwx.f` default = `-DTSERIES` (all 39
  files); `-DEPS1=1.0d-5` identical at slots 2 and 3, different (eps 0.100E-04) at slot 1; = `PartnersRivalsM1f run_local`
  at 5 cells. Replicate k alone = line k of the full task (12 of 12; w-line apart from `gcount`). Explicit seed base 40000000
  = the production call except the `# reps:` line. `kappa.x` = `ref/kappa.c`; `kappa_e5.x` = `kappa.x` at slots 2 and 3.
- `tasks_local.txt` through `task.sh`, then `collect.py --local`: every arm runs (0–7 s); a truncated t-file is detected.
- `submit.sh` dry run with fake `sbatch`/`squeue`: JOBID file and a ledger line with the job ID in field 4.

Changes made by the verifier:
1. `task.sh`: the games.dat staging used ONE shared temp name (`.games.dat.new`) in a directory that 71 E3 (E5) tasks, and
   pairs of probe rows, enter at the same moment: one task could rename another's half-written copy into place. Now a
   temp name unique per host/job/process, then an atomic rename. `--check` no longer writes games.dat.
2. `task.sh`: the TSERIES detection `*_ts*` missed the local test build `cw_tts.x` (its t-file was never checked); now `*ts*.x`.
3. `collect.py` (E3/E5): the enrichment denominators are now production's `kappa_e4` for the arm and for production alike
   (the kernel's law does not depend on eps). Before, the arm used the probe's independent KA estimate and production used
   kappa_e4, so kappa noise (and the 3/M floor applied to one estimate but not the other) could change "most enriched" with
   no change in any share. The KA masses remain the eps-independence check, and the table reports how many games would
   change with them.
4. `collect.py` (TS) and §8: the ESS criterion is the internal review notes' "at least 50 per replicate" (the minimum over replicates),
   not the median; the median is still printed.
5. `results_local/` regenerated with the fixed scripts; the local KA5 smoke row completed for the current `games24.txt`.
