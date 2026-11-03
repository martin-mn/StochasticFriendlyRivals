# `cannon/stratcert`: one exact friendly rival in every retained support stratum (internal review, round 2, item 1 / P1)

Remote directory: `~/CL1/PRFull/cannon/stratcert`. Everything is in the eps -> 0 limit.

## What it proves

The exclusion certificates of `check/ubound` prove that at most 364 (W), 5022 (N), 5022 (E) and 7726 (ray T>S)
open support strata can hold a friendly rival. The M2f census found friendly rivals, in double precision with a
tolerance of 1e-9 ("validated"), in 362, 5020, 5022 and 7721 of them; the re-censuses of W and N (and the C/D
mirror) found the other 2 + 2 + 1 (see Points). This kit replaces "found" by "proved".
For every retained stratum it certifies one exact rational point p of the open stratum that is:

* **efficient.** The exact efficiency class is checked three ways, and all three must agree:
  * `cls.bin` of `ub.c classes`, through the class column of `units.tsv` that `mkstrata.py` read from it (`cls.bin`
    itself is not staged, so the direct `cls.bin` lookup in `stratcert.py` is inactive on Cannon and in the kit root);
  * `en.efficiency_class(pattern)` (lemma L1);
  * the exact rational self-play support at p (`exact.py`).
* **competitive.** `ce.certify_point(p, sign)` decides the sign of sign·(w_DC − w_CD) against each of the 65 536
  binary co-players. Each decision is one of three: an outward-rounded interval enclosure of one strict sign; a
  structural zero (w_CD = w_DC = 0 by exponents); or the exact rational value, when the enclosure contains 0.
  Not one co-player may be negative. By the SI Proposition, this covers co-players of any memory.

With the stratum counts' proved upper bounds, the counts 364 / 5022 / 5022 / 7726 then become **exact**. So do the
126 + 238 strata of the two W families, p6 = p8 = 0 and p6 = p8 = 1.

## Units

`units.tsv` has one line per unit, where a unit is a (stratum, sign) pair. It holds 12 748 units:
* 7726 units with sign + (the ray T>S strata). Of these, class C are also W (364), class A are also N (5022),
  and class Z are the ray only (2340).
* 5022 units with sign − and class A (E).

Competitiveness depends only on the sign of T − S. Efficiency on the switch line is w_DD = 0 (class C, A or Z).
So one certificate per sign-+ unit serves both its wedge (W or N) and the ray.

How the list was built (`prep/mkstrata.py`, log `prep/mkstrata.out`):
1. `ub classes` computed the classes of all 3^16 strata.
2. `ub cover` ran with the 1681 co-players of `check/ubound/results/pool_c2_q.txt` (= `pool_q.txt`). A stratum is
   retained when none of them excludes it.
3. The result reproduces `check/ubound/results/status_c2_alldims.out` exactly: W 364 (126 at p6 = p8 = 0, 238 at
   p6 = p8 = 1, none elsewhere), N 5022, E 5022, ray+ 7726, with the same dimension histograms. It also confirms
   mirror(N) = E, and every census friendly-rival stratum of `M2f/check/st_eff.txt` is retained.

The rebuilt `cls.bin` and `c2.*` (180 MB) are not staged or kept. To regenerate them:
`cd prep && ./ub classes cls.bin 3 && ./ub cover cls.bin ../pool_q.txt c2 3`, about 2.5 min.

## Points

`points.tsv` holds the census points. They come from the M2f census, Cannon array 45751662:
`~/CL1/PartnersRivalsM2f/ca1/out/pt_*.txt`, which was read-only here.
* That census stored 5 784 521 confirmed points: up to 4 per stratum and sign, as `%.17g` doubles, so each is an
  exact binary rational. Of these, 50 897 belong here: sign P (T>S) of every class, and sign M of class A.
* Every one of them belongs to a unit.
* The md5 of the 512 source files is in `prep/pt_md5_cannon.txt`. Each point's source is recorded as `file:line`.
* 12 743 of the 12 748 units carry 1 to 4 census points. Five retained strata were missed by the census and have no
  census point: `010***01*1*1*010` and `*10***01*1*1*010` (N), `1*1100111*1**0*0` and `1*1*00111*1**0*0` (W),
  and `11****01*111*010` (ray, class Z). The re-censuses that found them (W: 512 points per candidate stratum,
  local; N: Cannon 45785345) did not store their points, so the kit finds its own witnesses for these five
  (the ladder below).
* Re-checked 2026-09-27: the 512 source files still have the md5 values of `prep/pt_md5_cannon.txt`; they hold
  5 784 521 lines (class A 40 154, C 4 981 110, Z 763 257); 30 random lines of `points.tsv` equal their
  `file:line` source as doubles.

## The ladder for each unit (`stratcert.py`)

Candidates are tried in order. The first certified candidate ends the unit.
1. `census:k`: the stored census points of the unit, in their stored order.
2. `mirror:k`: class A only. These are the census points of the C/D-mirror unit, p_s -> 1 − p_{15−s}, as exact
   rationals (they are not doubles).
3. `round10:k`, `round16:k`: the census points rounded to the grid 2^-10, and then to 2^-16.
4. `random:k`: up to 4096 dyadic draws n/2^20 in [0.02, 0.98], seeded by the uid.
   * A float screen proposes candidates: first the 1681 pool co-players, then all 65 536. The screen never counts
     as proof.
   * At most 8 screened draws go on to `certify_point`.
5. If no candidate is certified, the unit is marked `certified 0`, with the reason for the last failure.

Output: `out/t_TASK.tsv`, one line per unit. The columns are:
* `uid sign pattern dim class loci certified route source n_tried`;
* `n_pos_interval n_pos_exact n_struct_zero n_exact_zero n_neg n_unres_interval`. These are the co-players by
  decision route; the first four sum to 65 536;
* `min_margin` (exact), `selfplay` (exact w_CC, w_CD, w_DC, w_DD);
* `witness`: all 16 coordinates as exact fractions a/b;
* `seconds reason`.

A task is idempotent. Finished uids are skipped on restart, and a line cut short by a killed task is dropped and
redone. `t_TASK.done` marks a complete task.

## Validation (Mac)

* 45 units, covering all four loci and dimensions 0–8, were all certified by `census:0`. The five census-less units
  were certified too: the two N strata by `mirror:0`, and the two W strata and the Z stratum by `random`.
* `prep/validate.py` (`prep/validate.out`) checked 20 direct calls of the **original** library (`PRFull/check/cert`,
  `ce.certify_point(witness, wedge=locus)` for each locus of 12 units). All 20 agree with the kit's lines field by
  field: certified, counts per route, exact minimum margin, exact self-play and class.
* The ladder was tested with a non-competitive point put first. It is rejected with a negative exact margin, and
  the next candidate is certified. With no other candidate, the unit ends `certified 0` with the reason. With the
  random search enabled, it is certified by `random:k`.
* The Cannon login node gave byte-identical lines for 4 units, apart from `seconds`. Python there was 3.10.9, with
  numpy 1.26.4 and gmpy2 2.3.1; the Mac has gmpy2 2.1.2.
* After the probe, 8 probe units (the lowest- and highest-dimensional unit of each of the four kinds +C, +A, +Z,
  −A) were rerun on the Mac: all 8 lines equal the probe's compute-node lines in every column except `seconds`.
* The five census-less units, rerun on the Mac: certified in 1.8–8.6 s each, by `mirror:0` (the two N strata) and
  `random:2`, `random:71`, `random:16` (the Z stratum and the two W strata), one `certify_point` call each.

## Build and run

    ./compile.sh                  # sha256 of the library copy (SHA256SUMS), cert/libcert.so, and on Linux the
                                  # user-level install of numpy<2 and gmpy2 into ./pylib (pip --target; no Cannon
                                  # python module ships gmpy2 -- checked python/3.10.9, 3.12.x, intelpython/3.9.16)
    sbatch --parsable probe.sbatch    # tasks 0-3 of the 64-task partition (uid % 64), test partition -- DONE
    sbatch --parsable full.sbatch     # tasks 4-63, sapphire (0-3 are the probe's, already done) -- DONE (48981870)
    python3 collect.py out 64         # completeness (every uid exactly once, no malformed line, every .done), an audit
                                      # of every certified line (below), per-locus table, certified/retained by
                                      # dimension, the verdict line; writes out/stratcert_all.tsv

`collect.py` audits every certified line independently of the task: the witness has 16 coordinates in [0,1] and lies
in the open stratum (0/1 exactly where pinned, strictly inside (0,1) where free); the four decision counts add up to
65 536, `n_neg` = 0, `n_unres_interval` = `n_pos_exact` + `n_exact_zero`, and the minimum margin is >= 0; the recorded
exact self-play sums to 1, its class equals the unit's class and meets every locus, and the sign matches the loci;
and the exact self-play is recomputed from the witness with `cert/exact.py` (gmpy2, or `fractions.Fraction` under the
login node's python 3.6) and must equal the recorded one (`--no-recheck` skips this, and then the verdict is False).
The verdict is `True` only if nothing is missing, duplicated or malformed, all 64 `.done` markers exist, the audit
finds nothing, and every retained stratum of every locus is certified. (Before the review of 2026-09-27 the verdict
ignored duplicates, malformed lines, the `.done` markers and the line contents.)

`task.sh TASK NTASKS OUTDIR` runs `stratcert.py` on `$SLURM_CPUS_PER_TASK` single-threaded worker processes. Python is
`/n/sw/Mambaforge-22.11.1-4/bin/python3`, which is what `module load python/3.10.9-fasrc01` puts on PATH, with
`PYTHONPATH=./pylib` (numpy 1.26.4, gmpy2 2.3.1, installed by `compile.sh`; present on Cannon).

The library in `cert/` is a verbatim copy of `PRFull/check/cert`: `cmp` passes (re-checked 2026-09-27, after the
probe), and the sha256 values are in `SHA256SUMS`; `sha256sum -c SHA256SUMS` passes on Cannon. Local timing is
8.6 core-s per unit on average (0.9–19 s). The cost is almost entirely the exact rational evaluation of the
co-players whose enclosure contains 0. For the 45 units above, every one of these was an exact tie (margin exactly
0). None was a positive margin too close to call.

## Probe 48975432 (test partition, 4 tasks x 8 cpus, Sapphire Rapids)

`sacct -j 48975432 -X`: all 4 tasks COMPLETED, exit 0:0, elapsed 3:28, 3:28, 3:36, 3:33 (task 0 is the worst,
3:35.9 wall, 1651 user-s); all four `.err` files empty; each `t_N.done` present; MaxRSS 0.78–0.79 GB per task
(130 MB per worker). Outputs fetched to `out/` here.

* 800 units (tasks 0–3 = uid % 64 in {0,1,2,3}), **800 certified, 0 failed**, every one by `census:0` (the first
  stored census point, exactly as stored). Loci: W 17 (W00 1, W11 16), N 315, E 316, ray+ 484; dimensions 0–8.
* At every witness the exact minimum margin over the 65 536 co-players is exactly 0. Every co-player whose
  enclosure contained 0 was an exact tie: `n_pos_exact` = 0 in all 800 lines.
* Per-unit cost, single core (compute node), mean [max] seconds:

| kind (loci) | n | mean | median | max | by dimension, mean |
|---|---|---|---|---|---|
| +C (W, ray+) | 17 | 3.05 | 2.28 | 11.98 | d1 5.7, d2 1.8, d3 3.2, d4 2.4, d5 2.5 |
| +A (N, ray+) | 315 | 7.64 | 6.99 | 15.31 | d0 8.0, d1 7.7, d2 7.6, d3 7.8, d4 7.5, d5 6.9, d6 8.8, d7 8.4, d8 9.6 |
| +Z (ray+ only) | 152 | 9.86 | 9.16 | 16.06 | d0 11.4, d1 9.6, d2 10.5, d3 10.2, d4 10.0, d5 8.9, d6 7.5, d7 9.3, d8 10.7 |
| −A (E) | 316 | 7.88 | 7.25 | 16.30 | d0 6.2, d1 6.9, d2 7.8, d3 8.0, d4 8.1, d5 8.2, d6 7.2, d7 8.7 |

  Overall 8.06 core-s per unit on the compute node; the Mac is about 1.05x slower on the same units. The cost
  does not grow with the dimension: it is set by the number of exact ties, not by the stratum.

## Projection of the full run

Mean probe rate per (kind, dimension) times the unit counts (+A 5022, −A 5022, +Z 2340, +C 364):
* all 12 748 units: 1.02e5 core-s = **28.3 core-h**; the 60 remaining tasks (4–63): 9.55e4 core-s = **26.5 core-h**
  of compute, about 60 x 8 x 3.6 min = 29 core-h allocated;
* projected wall per task 195–203 s (probe tasks: projected 201–202 s, measured 202–212 s of compute);
* `full.sbatch`: `--array=4-63`, 8 cpus, `--mem=4G` (5x MaxRSS), `--time=00:15:00`. 2 x the worst probe task is 7:12;
  the limit is raised to 15 min because the probe never used a fallback route, and a unit whose four census points
  all fail could take up to 24 `certify_point` calls of at most 16 s (about 6.5 min). A task that still hits the
  limit resumes where it stopped when resubmitted (idempotent).
* `sbatch --test-only full.sbatch` (no job created; `sacct` has no row for the ID it prints) accepts the script;
  on 2026-09-27 21:10 it estimated the start at 2026-09-28 12:58 on sapphire (backfill of 15-min jobs usually beats
  such estimates) and 22:23 the same evening on `shared`. `shared` is Cascade Lake: its rate was not measured,
  so it would need its own probe first.

## Adversarial review (2026-09-27, before the full run)

Independent re-checks, all on the Mac with the original library `PRFull/check/cert` (scripts and logs in the review
session's scratchpad, not part of the kit):
* All 800 probe lines: uid % 64 = task, fields equal `units.tsv`, the witness lies in the open stratum and equals
  (exactly) the census point its `source` names, counts add up to 65 536 with `n_neg` = 0, the exact self-play
  recomputed from the witness equals the recorded one and its class is the unit's class and meets every locus.
  No problem.
* Eight witnesses (probe units 1475 W00 d1, 7680 W11 d5, 5056 N d8, 12675 E d7, 6976 Z d8, and the census-less
  7712 W `random:71`, 6259 Z `random:2`, 6897 N `mirror:0`, rerun here with `stratcert.py --only`):
  (a) `ce.certify_point(witness, wedge=locus)` for every locus agrees with the kit line field by field;
  (b) the exact rational margin against all 65 536 co-players, with no interval arithmetic, gives the same
  positive / zero / structural-zero counts and no negative margin;
  (c) an exact finite-eps chain (matrix-tree determinants, no leading-order arithmetic) at eps = 10^-60 for 32
  co-players per witness (structural zeros, exact ties, the 6 smallest and 6 random positive margins) differs from
  the limit margin by at most 7e-58, while the smallest positive limit margin is 4.5e-5 or more; the exact
  finite-eps self-play at eps = 10^-60 is within 4e-59 of the recorded limit, so each witness's class is confirmed
  without the leading-order code.
* The five census-less units rerun through the kit locally: certified in 1.9-9.0 s by the same routes as before.
* Fixed in the kit: `collect.py` (verdict now also requires no duplicate or malformed line, all 64 `.done`
  markers, and a clean audit of every certified line; NTASKS defaults to 64); README submission snippet (absolute
  path for `JOBID_full`, guard against a failed `sbatch`) and resubmission rule. `stratcert.py`, the library copy,
  `units.tsv`, `points.tsv`, `task.sh` and the batch scripts are unchanged, so `SHA256SUMS` and the probe stand.

## Job IDs

| file | array | tasks | state |
|---|---|---|---|
| `JOBID_probe` | 48975432 | 4 (test, 8 cpus each; ledger line in `~/jobs.tsv`) | COMPLETED 4/4, 800/800 certified |
| `JOBID_full` | 48981870 | 60 (tasks 4-63, sapphire, 8 cpus each; ledger line in `~/jobs.tsv`) | COMPLETED 60/60, with the probe 12 748/12 748 certified, verdict True |

## Full run 48981870 (sapphire, 60 tasks x 8 cpus, submitted 2026-09-27)

`sacct -j 48981870 -X`: all 60 tasks COMPLETED, exit 0:0, elapsed 0:03-0:04 per task (longest 3:21), first start
21:53:13, last end 21:57:23; 28.9 core-hours (CPUTimeRAW), TotalCPU 27.4 h, MaxRSS at most 775 MiB per task
(`PRFull/release/accounting`). All 60 `logs/full_*.err` are empty, and all 64 `out/t_N.done` markers exist. Outputs
fetched to `out/`.

`python3 collect.py out 64` on tasks 0-63 (printed output deposited as `collect.out`; about 6 s, rerun by
`PRFull/release/reproduce_certificates.sh`, step `stratcert`):
* units 12 748, lines 12 748, duplicated 0, malformed 0, missing 0; tasks without `.done`: 0; the audit of all 12 748
  certified lines (exact self-play recomputed for every witness) finds no problem;
* every retained stratum certified: W 364 (W00 126, W11 238), N 5022, E 5022, ray+ 7726, in every dimension;
  **verdict True**;
* routes: the first stored census point (`census:0`) in 12 743 units; `mirror:0` for the two census-less N strata,
  `random` for the two census-less W strata and the Z stratum;
* the exact minimum margin over the 65 536 co-players is 0 at all 12 748 witnesses, and no co-player needed the
  exact step for a positive margin (`n_pos_exact` = 0 throughout): every enclosure that contained 0 was an exact tie;
* 105 939 core-s recorded by the units (8.3 s per unit).

The full run was submitted with (it wrote `JOBID_full` here and on Cannon, and one ledger line, before anything else):

    K=/Users/martin/Documents/CL1/PRFull/cannon/stratcert
    jid=$(ssh -o BatchMode=yes cannon "cd ~/CL1/PRFull/cannon/stratcert && sbatch --parsable full.sbatch")
    if [[ "$jid" =~ ^[0-9]+$ ]]; then
      echo "$jid" > $K/JOBID_full; ssh -o BatchMode=yes cannon "echo $jid > ~/CL1/PRFull/cannon/stratcert/JOBID_full"
      ssh -o BatchMode=yes cannon "printf '%s\tPRFull\tstratcert\t%s\t%s\t%s\n' \"\$(date +%F)\" $jid 60 'full: one exact friendly rival per retained stratum, tasks 4-63 of 64 (sapphire)' >> ~/jobs.tsv"
    else echo "sbatch failed: '$jid'"; fi

A task that ends TIMEOUT, FAILED or OUT_OF_MEMORY is rerun by array index, which resumes it (finished uids are
skipped): `sbatch --parsable --array=<indices> full.sbatch` (the command-line `--array` overrides the script's). That
is a new array job ID: record it as `JOBID_full2` (here and on Cannon) with its own ledger line, exactly as above, and
check both IDs with `sacct`. Look at the task's `.err` first: a Python exception in a worker ends the task and will
recur at the same unit.
Resubmit an index only once `sacct` shows its previous run finished: two runs of one index at the same time append to
the same `t_N.tsv` (`collect.py` would then report duplicates and the verdict False).

Completion only from `sacct -j $(cat JOBID_full) --format=JobID,State,ExitCode -X --noheader | awk '$2 != "COMPLETED"'`
(empty), empty `logs/full_*.err`, then `python3 collect.py out 64` must print `missing 0`, `tasks without .done: 0` and
the verdict `True`.
