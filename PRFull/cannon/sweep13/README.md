# `cannon/sweep13`: the rarity sweep at exactly α = 1/3

This kit redoes the rarity sweep of SI Table 3, Figure 3d, Methods "Rarity" and the main-text rarity paragraph
of the paper. That sweep drew 1 633 280 000 memory-two strategies from the three-outcome kernel: each
coordinate is 0 with probability α, 1 with probability α, and uniform on (0,1) otherwise. Each strategy was
classified exactly in the ε → 0 limit:

- efficient below or above the switch line;
- competitive for T > S and for T < S, against all 65536 binary co-players;
- friendly rival in each wedge W, S, N, E.

It ran at **α = 0.333**, not 1/3 (Cannon array 45831015, 512 tasks × 3 190 000 draws, 1935.2 core-hours). The
program read α as a permille (`atol(argv[4])/1000.0`). This kit runs the same campaign at **exactly α = 1/3**, the
kernel of the evolutionary runs. Everything else is the same: the same size, engine, classifier and tolerances.

**Status: run on Cannon.** The probe (array 48928481, `JOBID_probe`, output in `out_probe/`) and the full run
(array 48928780, `JOBID_13`, 512 tasks, output in `out_13/`) are done; `collect.py` collected the run into
`sweep13.json`, and SI Table 3 and the rarity numbers of the paper come from it.

## Files

| file | what |
|---|---|
| `census.c` | `PartnersRivalsM2f/cannon/alpha/census.c` plus the patch below (`census.patch`). |
| `m2f.c`, `witnesses.h` | The ε → 0 engine and the co-player screen. Byte-identical copies of the archived kit's files. |
| `orig/` | The archived kit, untouched: `census.c`, `compile.sh`, `task.sh`, `run_333.sbatch`, `probe.sbatch`, `collect.sh`, from `PartnersRivalsM2f/cannon/alpha/`. Every file except `collect.sh`, which was local-only, is byte-identical to the Cannon copy `~/PartnersRivalsM2f_alpha` inside `PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz`. The JOBID files are left out on purpose: those arrays are not this kit's. |
| `compile.sh` | `gcc -O3 -march=native` (`CC=`/`CFLAGS=` override it). Writes `BUILD.txt`: compiler version and sha256 of the sources and binary. Creates `out/` and `out_probe/`. |
| `task.sh` | One array element. The full run takes `ALPHA`, `NSWEEP` and seed `SEEDBASE + i`. The probe takes one line of `TASKFILE`. |
| `run_13.sbatch` | The full run: sapphire, `--array=0-511`, 1 CPU, 4G, 8 h, `ALPHA=1/3 NSWEEP=3190000 SEEDBASE=3000000`. |
| `probe.sbatch`, `tasks_probe.txt` | The probe: partition `test`, 4 tasks, 30 min. |
| `ref/probe45830542_task4.{tsv,err}` | The archived output of probe 45830542 task 4, the target of probe task 0. |
| `collect.py` | SI Table 3 and the rarity numbers of the text, from the task outputs. It also has the reproduction mode (`--reproduce`), the probe mode (`--probe`) and the sacct gate (`--jobid`). |
| `stage.sh`, `submit.sh` | Copy and build the kit on Cannon. Submit under the JOBID protocol. |
| `validate/` | The local checks below: `validate_local.sh` (7 min on 2 cores; writes only `validate/build/`), plus `replay.c` (stream replay against archived W points), `invert.c` (seed recovery by inverting the xorshift) and `kthr.c` (exact threshold test). The replay of the efficiency test: `rcheck.c` (with its copy of `witnesses.h`), `rcheck.sh`, `rcheck_compare.py` and its results `rc.tsv`, `flags.tsv` (see below). |

## The patch (`census.patch`, psweep only; `strata` untouched)

- **The permille form, `333`, is unchanged.** It is still the threshold `u < am/1000.0` on the double `urnd()`.
  So every archived psweep run reproduces bit for bit (see Validation).
- **The exact-rational form, `num/den` (e.g. `1/3`), uses no floating point to decide a pin.**
  - Let k be the 53-bit integer behind `urnd()` (so `urnd() = k/2^53`), taken by a new `k53()` with the same
    xorshift step, and let q = ⌊2^53/den⌋.
  - k is redrawn while k ≥ den·q. Then **k < num·q gives 0** and **k < 2·num·q gives 1**.
  - Each has probability num·q/(den·q) = num/den **exactly**.
  - At 1/3: q = 3002399751580330, and only the 2 largest of the 2^53 values of k are redrawn (probability
    2.2e-16).
  - At the same α the two forms draw the same strategies from the same stream, except with probability about
    2·10^-13 per coordinate. `333/1000` and `333` differ only on the few hundred values of k between the rounded
    thresholds and on the 992 redrawn values. The probe checks this (task 1 against task 0).
  - The interior value is `urnd()`, redrawn in the probability-2^-53 case that it is exactly 0.0. The interior
    therefore lies strictly in (0,1), and a 0 is always a pin.
  - For contrast, a double threshold `u < 1.0/3.0` would pin 0 with probability 1/3 + 3.7e-17.
- The ALPHA line carries α as it was given (`333` or `1/3`), so the collector can tell the two apart. With `num/den`,
  the log has one more line: `psweep exact alpha=1/3: <k> redraws`.
- `pdraw <n> <seed> <alpha>` is a new mode. It runs the same `kdraw()` with the same seeding and does no
  classification. It is used for the law check.
- This follows `cannon/wfc` in taking α as `num/den`. Unlike `wfc.c`, which divides two `atof`s and compares
  doubles, the comparison here is on the random integer, as asked.

## Seeds

- The seed map is unchanged: `rs = seed·2862933555777941757 + 3037000493` (xorshift64).
- The seeds used before:
  - 45831015: `100000 + 1000·AM + i` = 433000 + i.
  - Probe 45830542: 90000 + i, over its α grid (0, 100, 200, 300, 333, 360, 400, 500). The probe's `task.sh`
    was not archived. The seeds were recovered here by inverting the xorshift from an interior coordinate of
    an archived W point. They were then confirmed by replay on its three non-empty W files (α = 0.36, 0.4,
    0.5) and by reproducing task 4's output exactly.
  - Every permille seed of the old `task.sh` lies in [100000, 600999].
- **The full run here uses seeds 3000000 + i, i = 0..511, and the probe's cost tasks use 3900000 and 3900001.**
  None of these is a seed of any archived psweep run: 45831015, 45830542, or the permille family. So the new
  sample is independent of the 0.333 sample it will be compared with.
- Two streams overlap only if two seeds land within about 7·10^7 steps of each other on the one 2^64 − 1 cycle.
  Over about 10^3 streams that has probability about 10^-5.

## Validation (local, Apple M2 Pro, Apple clang 15, at most 2 cores, nice 10)

`bash validate/validate_local.sh [<extracted PartnersRivalsM2f_alpha/out>]` re-runs checks 1–3 and 5.

1. **Bit-identical to the archived Cannon output, built with gcc 8.5.0 on x86-64.**
   - `census psweep 100000 90004 333` reproduces probe 45830542 task 4 exactly. Both the ALPHA line and the
     stderr summary match, for all 11 counts over 10^5 draws at α = 0.333 (`ref/`, `cmp` identical).
   - The same draws through the rational form, `333/1000`, give identical counts, with 0 redraws.
2. **The draw stream of the archived full run.**
   - For every one of the 512 tasks of 45831015, the patched `kdraw()` was replayed with seed 433000 + i
     (`validate/replay.c`).
   - All **6388 of 6388** archived W friendly rivals (`pt_333_i.txt`, `%.17g`) turn up in order, equal in all
     16 coordinates.
   - The complete patched program, classification included, then reproduces the start of two tasks' W files
     exactly:
     - task 232, first 3000 draws: W at draw 126, line 1;
     - task 55, first 9000 draws: W at draws 631 and 8327, lines 1–2.
   - The program writes only totals, so a shorter run's output is not a prefix of a full task's output.
     The W file is, and the probe check in 1 covers all counts over 10^5 draws.
3. **The law at exactly 1/3.**
   - `pdraw 1000000 3000000 1/3` gives these fractions over 1.6·10^7 coordinates:

     | 0 | 1 | interior | χ²(2) | p | per coordinate, χ²(32) | p |
     |---|---|---|---|---|---|---|
     | 0.333353 | 0.333426 | 0.333221 | 1.04 | 0.60 | 25.5 | 0.79 |

   - At 10^8 draws the fractions are 0.333341 / 0.333333 / 0.333326. Against the α = 1/3 law, χ² = 0.59
     (p = 0.74). Against the α = 0.333 law, χ² = 3271.
   - There were 0 redraws.
   - The pin thresholds are exact integer comparisons. Setting the state so that `k53()` returns a chosen k:
     - k = q − 1 gives 0, and k = q gives 1;
     - k = 2q − 1 gives 1, and k = 2q gives an interior value;
     - k = 3q − 1 gives an interior value, and k = 3q and 2^53 − 1 are redrawn.
4. **The tasks run.**
   - `task.sh` was rehearsed under simulated Slurm variables, in a copy of the kit built by `CC=cc bash
     compile.sh`. It ran the full-run form (1/3, 3000 draws, seeds 3000000 and 3000001) and all four
     `tasks_probe.txt` lines at their real sizes.
   - Every task exited 0, the logs were as specified, and the stand-in Slurm stderr was empty.
   - `collect.py --probe out_probe` printed IDENTICAL for tasks 0 and 1 and exited 0.
   - `collect.py out --tasks 2 --n 3000 --alpha 1/3 --seedbase 3000000` passed every check.
   - It refuses, as it should, on a wrong alpha, a wrong seed base or a non-empty stderr file.
   - Local cost (M2 Pro, nice 10): 2.86 ms per draw at 0.333 (10^5 draws) and 3.08 at 1/3 (2 × 2·10^4).
5. **`collect.py --reproduce` on `PartnersRivalsM2f/check/sweep_a333.tsv` reproduces the published numbers:
   41 of 41 agree.**
   - The four SI Table 3 rows, character for character.
   - The legend marginals (0.1442, 0.0648, 0.01827).
   - The main-text rates, factors and intervals (673 [657, 690], 12 [11.8, 12.0], 5.1 [5.095, 5.097]).
   - The Methods N − E z = 0.77 and the dispersion range 0.95–1.12.
   - The old Table 4's 0.023 %, 255,700 and 74.
   - Every internal review (round 1) rate, Clopper–Pearson interval, delta-method ratio interval and task-bootstrap
     interval, the bootstrap bit-identical to `scratchpad/rarity_ci.py`.
   - On the archived per-task directory it also gives 41 of 41, with all 512 logs, W files and
     `j_333_45831015_*.err` consistent:

         python3 collect.py <archive>/PartnersRivalsM2f_alpha/out --glob 'a_333_*.tsv' \
                 --errglob 'j_333_45831015_*.err' --alpha 333 --reproduce

     - The directory holds the same 512 lines as the tsv.
     - The tsv is in Cannon's en_US glob order (0, 100–109, 10, 110, …), and the bootstrap resamples by
       position. So in the directory's index order the bootstrap agrees only to Monte Carlo precision
       (≤ 3·10^-4 relative).
     - Every other number is identical.

## Replay of the efficiency test (`validate/rcheck.c`)

The census reads a draw's self-play class from its limit weights with thresholds 1e-12 (`selfclass()` in `census.c`).
SI Box 1 and the Computational Supplement ("The efficiency test") report how often that reading differs from the exact
class. `validate/rcheck.c` measures it on the run at exactly α = 1/3 (array 48928780, `out_13/`):

- It regenerates each task's 3 190 000 draws from its seed 3000000 + i, with `kdraw()` and the seeding copied
  verbatim from `census.c`.
- It computes each draw's self-play limit once and reads its class twice: from the weights, exactly as
  `selfclass()`, and exactly from the support Σ, the states that carry the least exponent (no tolerance).
- A draw is flagged when the two readings differ. With `classify`, a flagged draw also gets the psweep
  competitiveness test (the screen and the rival margin against all 65536 binary co-players, tolerance −1e-9), to
  see which friendly-rival count it entered.

Build and run, from `validate/` (C99, no other dependency; about 4.5 s per task on an M2 Pro as an arm64 binary):

    cc -O2 -ffp-contract=off -fno-fast-math -o rcheck rcheck.c -lm
    ./rcheck <task> 3190000 classify > replay_<task>.out

`JOBS=3 bash validate/rcheck.sh [TASK ...]` does this for the listed tasks (default all 512: 47 core-minutes), in
`validate/build/`. All 512 take about 17 min on 3 cores in the step `rcheck` of `reproduce_certificates.sh`, whose
`cc`, called from the author's x86_64 Python, builds an x86_64 binary that runs under Rosetta: 5.8 s of CPU per task
(3 tasks timed on 2026-11-03: 6.8 s wall, 17.3 s of CPU, with the build). It then runs
`validate/rcheck_compare.py`, which compares each replay with the deposit: the RC line with `rc.tsv`, the FLAG lines
with `flags.tsv`, and the read counts with the task's ALPHA line `out_13/a_1o3_<task>.tsv`. A complete replay prints

    replayed 512 tasks: RC lines differ 0, FLAG lines differ 0, ALPHA counts differ 0

The step `rcheck` of `PRFull/release/reproduce_certificates.sh` runs it on tasks 1 and 3 (`--quick`) or on all 512.

- `rc.tsv` (no header) has one RC line per task, as `rcheck` prints it. The columns are: `RC`, task, n (draws);
  the classes read from the weights, C (efficient below the switch line, W_CC = 1), A (above it,
  W_CD + W_DC = 1) and Z (on it, W_DD = 0); the same three read exactly from Σ; the number of flagged draws; and the
  flagged draws that entered the friendly-rival counts of W, S, N and E (read as C or A and competitive for T > S or
  T < S).
- `flags.tsv` has a header line and one FLAG line per flagged draw: task, draw index, the class read and the exact
  class (`-` for none), the least weight in Σ, the four weights, the results of the two competitiveness tests, and
  the 16 coordinates (`%.17g`).
- Totals: 512 tasks, 1 633 280 000 draws, 137 flagged draws in 120 tasks; contaminated W/S/N/E 0/15/0/0. The read
  counts of every task equal the ALPHA line of the Cannon run.
- 18 of the 137 are competitive for T < S: 15 of them were read as C, and so entered the S count, and 3 as Z.
  None is competitive for T > S.

## Sizing (from 45831015, `sacct`)

- 45831015: 512 tasks, 3:13–3:55 per task, 3:58 from first start to last end, 1935.2 core-hours, 4.27 ms per
  draw on sapphire, 17 MiB peak RSS.
- At 1/3, 0.1 % more coordinates are pinned than at 0.333. Locally, 1/3 cost 3.08 ms per draw (2 × 2·10^4
  draws) against 2.86 at 0.333 (10^5 draws), 8 % more. That may be real or scatter: the cost is dominated by the
  full scans of screen survivors, about 750 per 2·10^4 draws. The probe measures it on the cluster.
- The new binary is built by gcc 16.2 rather than 8.5, and its speed is measured by the probe, not guessed.
- `--time=8:00:00` is twice the longest old task. It holds unless the new build is more than about 1.9× slower
  per draw.
- Expect about 1950–2100 core-hours and about 4 h wall-clock if all 512 tasks start together.
- Output is about 5 MB in all. Every file is tiny: W files hold about 12 lines.

## Protocol (`~/.claude/CLAUDE.md`)

Martin opens `ssh cannon`, and the session keeps it warm (`ssh -o BatchMode=yes cannon true` every ~25 min)
while any array is outstanding.

    bash stage.sh              # mkdir + rsync to ~/CL1/PRFull/cannon/sweep13, then compile.sh there (BUILD.txt)
    bash submit.sh probe       # = jid=$(ssh cannon "cd ~/CL1/PRFull/cannon/sweep13 && sbatch --parsable probe.sbatch")
                               #   -> JOBID_probe here and there + ~/jobs.tsv line (date PRFull sweep13.probe jid 4 note)
    # when squeue -j $(cat JOBID_probe) is empty:
    python3 collect.py --probe out_probe --jobid JOBID_probe --fetch CL1/PRFull/cannon/sweep13/out_probe

The probe mode checks three things:

- **sacct per task**: every task COMPLETED, and `squeue -j` is empty.
- **Empty Slurm stderr.**
- **The two bit-identity checks.** `p333_s90004` and `p333o1000_s90004` must both print IDENTICAL.
  - That shows the gcc 16.2 build reproduces the gcc 8.5 archived probe on the cluster.
  - If they differ, the most likely cause is a floating-point contraction difference between the two
    compilers. Rebuild with `CC=/usr/bin/gcc bash compile.sh` (the system gcc, which built the archive), or
    compare with `-ffp-contract=off`, before running the full array.

The probe also prints the measured ms per draw at 1/3 and the projected task time. `--time` holds if the
longest projected task times 1.5 is under 8 h.

    bash submit.sh 13          # 512 tasks > 100: prints mnowak's pending/running count and stops; report it
    bash submit.sh 13 --yes    # once agreed: sbatch --parsable run_13.sbatch -> JOBID_13, ledger line
                               #   (date PRFull sweep13 jid 512 'rarity sweep at exactly alpha=1/3, ...')
    # when squeue -j $(cat JOBID_13) is empty:
    python3 collect.py out_13 --alpha 1/3 --seedbase 3000000 --jobid JOBID_13 \
            --fetch CL1/PRFull/cannon/sweep13/out \
            --compare /Users/martin/Documents/CL1/PartnersRivalsM2f/check/sweep_a333.tsv --json sweep13.json

`collect.py` refuses to report unless all of the following hold:

- `sacct -j <JOBID> -X` shows every task COMPLETED (a FAILED task is investigated against its files, not rerun
  blindly);
- there are exactly 512 ALPHA lines of alpha `1/3` and n 3 190 000;
- task indices run 0..511, with seed 3000000 + i and no seed repeated;
- each census summary line equals its ALPHA line;
- each W file has FR_W lines;
- all 512 `j_13_<JOBID>_*.err` are empty.

`--compare` gives each rate at 1/3 against 0.333 with a 95 % interval. It prints the census-predicted shifts
(W +0.88 %, S +0.12 %, N and E +0.69 %) beside the measured ones.

Without the helper scripts, the submission is exactly:

    jid=$(ssh cannon "cd ~/CL1/PRFull/cannon/sweep13 && sbatch --parsable probe.sbatch") && echo "$jid" > JOBID_probe
    ssh cannon "cd ~/CL1/PRFull/cannon/sweep13 && echo $jid > JOBID_probe && printf '%s\tPRFull\tsweep13.probe\t%s\t4\t%s\n' \$(date +%F) $jid 'probe on test: repro of 45830542 task 4 + cost at alpha=1/3' >> ~/jobs.tsv"
    ssh cannon "squeue -u mnowak -h -t PENDING,RUNNING -o '%A' | wc -l"        # courtesy check, report before 512
    jid=$(ssh cannon "cd ~/CL1/PRFull/cannon/sweep13 && sbatch --parsable run_13.sbatch") && echo "$jid" > JOBID_13
    ssh cannon "cd ~/CL1/PRFull/cannon/sweep13 && echo $jid > JOBID_13 && printf '%s\tPRFull\tsweep13\t%s\t512\t%s\n' \$(date +%F) $jid 'rarity sweep at exactly alpha=1/3, 512x3.19e6 draws, seeds 3000000+i' >> ~/jobs.tsv"
    ssh cannon "sacct -j $(cat JOBID_13) --format=JobID,State,ExitCode -X --noheader | awk '\$2 != \"COMPLETED\"'"   # must be empty

## After the run (not part of this kit)

- `Figures/fig3.py`, and `PartnersRivalsM2f/Figures/fig1.py` if it is reused, read `sweep_a333.tsv`. `fig3.py`
  asserts `f[1] == "333"` and `fig1.py` parses `int(f[1])`, so both need a reader that accepts `1/3`, pointed at
  `cat out_13/a_1o3_{0..511}.tsv`. `fig3.py` also asserts that the recomputed rates and independence products
  equal its hard-coded 0.333 values (`TAB4`) to 6e-4 relative, and plots `TAB4` rather than the recomputed rates:
  at 1/3 the expected shift in W (+0.88 %) breaks that assert, so `TAB4` must be replaced by the new numbers too.
- The text sentences that carry 0.333 and the carried-to-1/3 shifts can then state the measured rates at 1/3
  directly:
  - Methods "Rarity", including the `% TODO` in the paper;
  - SI Table 3 and its legend;
  - SI Box 1;
  - the Fig. 3d legend;
  - the main-text paragraph;
  - SI Table 6's cost row (now CS Table 2);
  - the Methods sentence on the friendly-rival kernel mass.
- Storage: keep the kit on Cannon, archive `out/` locally and verify it, and delete only with Martin's OK.
