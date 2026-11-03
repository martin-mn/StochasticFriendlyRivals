# `cannon/evoaudit` — direct interval/exact audit of the memory-two evolutionary atoms (internal review, round 2, item 2 / P2)

Production (`PartnersRivalsM2f/cannon/wfc/wfc.c` with `m2f.c`) gave every distinct sampled strategy of the memory-two
Wright–Fisher runs its atom (efficient, stable, competitive; `atom = 4E + 2S + C`) in the limit ε → 0 in double
precision: efficiency from the self-play limit (`W_CC > 1 − 1e-12` below the switch line, `W_CD + W_DC > 1 − 1e-12`
above it), competitiveness from the 65 536 binary co-players with a margin counted non-negative from −1e-9, stability
from the same co-players with a deviation counted as a gain if π(q,p) > π(p,p) + 1e-9, and an efficient competitive
strategy set stable by the theorem of PR2026 without a scan. This kit re-decides every bit for a stratified random
sample of the stored strategies with no tolerance at all.

## What is decided, and how (`src/audit.py`)

For a stored strategy p (the 16 doubles of the dump, taken as exact rationals) at its game (R,S,T,P) = (1, cs, ct, 0)
(the exact doubles of `games.dat`; no game lies on the switch line S+T = 2 or on T = S):

| bit | rule | arithmetic |
|---|---|---|
| E | support class of the exact self-play limit: `C` (w_CC = 1) required below the switch line, `A` (w_CD + w_DC = 1) above it | exact rationals (`ex.lead_exact(p,p)`), cross-checked with `en.efficiency_class` of p's stratum (lemma L1) |
| C | sign(T−S)·(w_DC − w_CD) ≥ 0 against every binary q | one outward-rounded interval evaluation per q (`en.i_scanW`), margin enclosure formed as `certeng.c`'s `imargin`; structural zero when no support state has recent outcome CD or DC; exact rationals where the enclosure contains 0 |
| S | π(q,p) ≤ π(p,p) against every binary q, π(p,p) = R w_CC + S w_CD + T w_DC of the exact self-play limit, π(q,p) = R w_CC + T w_CD + S w_DC + P w_DD in the focal-p frequencies of the (p,q) chain | the same enclosures, one outward step after every operation; structural when the limit support holds a single outcome (then π(q,p) is that outcome's payoff exactly); exact rationals where the enclosure contains 0 |

Weak inequalities: a tie is not a gain and satisfies competitiveness. An enclosure containing 0 is sent to exact
arithmetic only when the bit depends on it (C when no margin is certified negative, S when no gain is certified
positive); otherwise the bit is already decided. Stability is scanned for every strategy, 111 included, so the audit
also tests the theorem on the sample. The co-player quantifier (binary q suffice) is the paper's box-vertex lemma and
is not re-examined here; the audit checks the arithmetic of the classification, not the lemma.

The production classifier is re-run on the same doubles (`src/prodcls.c`, which links production's `m2f.c`
unchanged and writes wfc.c's expressions verbatim; built with production's compiler and flags, GCC 8.5.0
`-O3 -march=native`) for two purposes: to check that the dump's atom is reproduced, and to name the co-players that
decided each production bit and their double values (`pr` in the output; `n_neg_tol` / `n_gain_tol` count the q at
which production's tolerance was engaged: a double margin in [−1e-9, 0), a double gain in (0, 1e-9]).

## The sample (`src/sample.py` → `sample/sample.npy`, `sample/design.json`)

Population: all 4 772 732 records of the production dumps (F1, N = 1000: 1 453 629; F2, N = 100: 3 319 103); a record
is one distinct sampled strategy of one replicate at one game, `count` its samples (40 000 samples per game). Strata:
regime × wedge × production atom (54 non-empty). Poisson sampling, seed 20260927: record i is included with the
exactly known probability π_i = 1 − (1−π^G_i)(1−π^A_h)(1−π^P_i), the union of

* G — random individuals of each wedge's population: π^G_i = min(1, b_w c_i / C_w), b_w = 3000 in W, N, E and 1000 in S
  per regime (C_w = 40 000 × games of the wedge);
* A — equal probability per strategy within the stratum: min(1, a_h / N_h);
* P — probability proportional to samples within the stratum: min(1, a_h c_i / C_h);

with a_h = 150 for the production atoms 110, 111, 011, 001 (oversampled) and 50 otherwise. Drawn: 27 261 strategies
(F1 12 023, F2 15 238; by production atom 110: 14 855, 111: 2 901, 011: 1 864, 001: 2 615, 000: 2 697, 100: 1 459,
010: 870). Every stratum with at most its target is taken whole (e.g. F1 W 111: all 50).

Estimates (`src/analyze.py`): the population share with property y (changed atom, a changed bit, unresolved) is the
Horvitz–Thompson sum Σ y_i (c_i / C)/π_i with its Poisson-sampling variance; when no audited strategy of a group has y,
the 95 % upper bound −ln 0.05 / κ, κ = min over the group's records with π_i < 1 of −ln(1−π_i)·C/c_i, holds whichever
records carry y. With the design above these bounds are 0.06 % per regime (F1 0.0603 %, F2 0.0582 %), at most 0.092 % in each of W, N, E
and 0.25 % in S (F1 0.2450 %, F2 0.2366 %) (`analyze.py`, `bound95_if_none`; recomputed independently in the review of
2026-09-27).

## Files

| | |
|---|---|
| `lib/cert/` | copy of `PRFull/check/cert` (engine.py, exact.py, certeng.c, compile.sh, `__init__.py`) |
| `src/m2f.c` | production's `m2f.c`, byte-identical (sha256 511bf5b1…e401) |
| `src/prodcls.c` | wfc.c's `classify` with full scans (reproduction and deciding co-players) |
| `src/load.py`, `src/sample.py`, `src/audit.py`, `src/analyze.py` | loaders, design, audit, estimates |
| `src/project.py` | per-strategy cost from the probe and the projection of the full run |
| `src/test_local.py`, `src/test_exact_full.py` | local validation (below); logs in `logs/test_*.txt` |
| `src/probe_version/` | `audit.py` and `analyze.py` exactly as they ran in the probe (see Status for what changed since) |
| `check.sh <arm>` | completion check of an array from its `JOBID_<arm>`: sacct per task, empty stderr, complete chunks, no errors |
| `games.dat` | production's (= `cannon/wfc2/games.dat`) |
| `data/F1`, `data/F2` | local only (not staged, and not in the release): read-only extraction of the dumps `S_1o3_<task>_1.bin` and `r_1o3_*` from `PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz` (`PartnersRivalsM2f_wfc/out/`) and `cannon_out_f2_2026-09-26.tgz` (`out_f2/`); `PRFull/release/tools/evoaudit_dumps.py extract` writes the `S_` files again from the archives (byte-identical, checked 2026-09-27) |
| `compile.sh`, `task.sh`, `probe.sbatch`, `run_full.sbatch`, `stage.sh`, `collect.sh` | the kit |
| `JOBID_probe`, `JOBID_full` | the probe array (48975385) and the full run (48981982) |
| `out/`, `results/` | the full run's chunks 4-127 (fetched) and the analysis of all 128 chunks (`summary.json`, `confusion.txt`, `records.tsv.gz`, and `analyze.out`, the printed output of `analyze.py`) |
| `out_probe/`, `results_probe/` | the probe's chunks 0-3 (fetched, md5 identical to Cannon) and their (partial) analysis |

Output: `<outdir>/chunk_<k>.jsonl`, one JSON line per strategy (production atom, reproduced production atom, validated
atom, per bit the method and the counts, the exact extreme margin and gain with their co-players, the production
diagnostics, and for a changed bit the margin or deviation responsible), and `chunk_<k>.done`. A task is idempotent:
records already written are skipped, a line cut by a kill is dropped and redone (tested locally). A strategy whose audit
raises an exception is written with `val: null` and the error text and counted as UNRESOLVED by the analysis (it is
never dropped; the `.done` line counts the errors); the task still succeeds.

Chunk k holds the sample rows i with i mod 128 = k (213 strategies for k <= 124, 212 for k = 125..127).

## Local validation (Mac, 3 threads, each run < 2 min)

* `src/test_local.py 24 1` (`logs/test_local_1.txt`) and `src/test_local.py 10 5` (`logs/test_local_2.txt`, after the
  post-probe changes): one strategy from each of 24 and 10 random strata. For each: (1) the margin classes formed
  from `en.i_scanW` equal the library's own `en.i_scan`, co-player by co-player; (2) the mirror identity used for
  π(q,p) holds exactly at 60 random binary q (R w_CC + T w_CD + S w_DC of the (p,q) chain = q's payoff
  R w'_CC + S w'_CD + T w'_DC in the (q,p) chain, `ex.lead_exact` both ways); (3) every sign decided by an interval
  or a structural rule agrees with the exact sign at those q; (4) production atom = reproduced atom = validated atom.
  All passed; nothing flagged.
* `src/test_exact_full.py 4 3` (`logs/test_exact_full.txt`): brute force, exact rationals against all 65 536 binary q
  and no interval or structural rule, for one production 110, 111, 011 and 001 each: the four atoms agree with the
  audit's, and so do the counts of negative / zero margins and of gains / ties wherever the audit resolved every
  enclosure of the bit (it skips the exact step for a bit already decided by a certified sign).
* restart: a chunk cut in the middle of a line and rerun completes, dropping the cut line and redoing it; a forced
  exception is recorded as an unresolved strategy and the analysis counts it as such.
* Mac against Cannon: the 27 probe strategies of rows i mod 1024 = 0, re-audited locally, give output identical to
  Cannon's field by field.  This is not general for production's double diagnostics (`pr`): the Mac build (clang,
  no FMA) can differ from Cannon's (GCC 8.5.0 -O3 -march=native) in the last bit, e.g. sid 3712 has pr.pp and pr.gmax
  one ulp apart (same atom, same counts).  The reproduction check is judged on Cannon only.

## Environment on Cannon

`python/3.12.8-fasrc01` has no numpy and no gmpy2 (nor do the other python modules; `intelpython/3.9.16` has numpy but no
gmpy2). The kit therefore has its own venv, made once on the login node:

    cd ~/CL1/PRFull/cannon/evoaudit && module load python/3.12.8-fasrc01 && python3 -m venv venv && ./venv/bin/pip install numpy gmpy2

(numpy 2.5.3, gmpy2 2.3.1; `exact.py` then uses gmpy2.mpq; without gmpy2 it falls back to fractions.Fraction, about
4x slower, same results). `task.sh` calls `./venv/bin/python` directly, no module needed on the compute nodes.
C: `/usr/bin/cc` = GCC 8.5.0.

## Protocol

    bash stage.sh                                   # rsync (no data, outputs, venv, Mac binaries) + compile on Cannon
                                                    # NOBUILD=1 bash stage.sh: copy only, keep the Cannon binaries
    jid=$(ssh -o BatchMode=yes cannon "cd ~/CL1/PRFull/cannon/evoaudit && sbatch --parsable probe.sbatch")
    echo $jid > JOBID_probe   # and on Cannon, and one line in ~/jobs.tsv
    bash check.sh probe                             # sacct per task, empty stderr, complete chunks: all sections empty
    bash collect.sh out_probe && python3 src/analyze.py out_probe results_probe && python3 src/project.py out_probe

Full run (submitted 2026-09-27 as 48981982, see Status): chunks 4-127 (`run_full.sbatch`, 124 tasks, sapphire, 1 core,
1G, `--time=00:11:00`).
Before submitting, the courtesy check of the global protocol (124 > 100 tasks):
`ssh cannon "squeue -u mnowak -h -t PENDING,RUNNING -o '%A' | wc -l"`, and ask if a large campaign is queued. Then

    NOBUILD=1 bash stage.sh
    jid=$(ssh -o BatchMode=yes cannon "cd ~/CL1/PRFull/cannon/evoaudit && sbatch --parsable run_full.sbatch")
    echo "$jid" > JOBID_full; ssh -o BatchMode=yes cannon "cd ~/CL1/PRFull/cannon/evoaudit && echo $jid > JOBID_full && \
      printf '%s\tPRFull\tevoaudit\t%s\t%s\t%s\n' \"\$(date +%F)\" $jid 124 'referee-2 P2 audit of memory-two atoms: full run, chunks 4-127 of 128' >> ~/jobs.tsv"
    bash check.sh full                              # after sacct shows every task finished
    bash collect.sh out && python3 src/analyze.py out_probe,out results

A task that timed out or failed is completed by resubmitting only its index (`sbatch --parsable --array=<k> run_full.sbatch`,
with its own JOBID file and ledger line); nothing already written is recomputed.

## Status

**Probe: JOBID_probe = 48975385** (sapphire, 4 tasks = chunks 0-3, 852 strategies, one core each), submitted 2026-09-27,
ledger line written. `sacct`: all 4 tasks COMPLETED, exit 0:0, Elapsed 5:22, 4:40, 4:51, 5:09 (TotalCPU 5:16, 4:35, 4:47,
5:02), MaxRSS 45-49 MB; all four `.err` empty; every `.done` reports 213 of 213 records. Outputs fetched, md5 identical.

Probe results (852 strategies, both regimes, all four wedges, every production atom; `results_probe/`):

* validated atom = production atom for **all 852**; production's classifier reproduced the dump's atom for all 852;
  no validated 101 (the theorem efficient + competitive => stable held on all 90 production-111 strategies, whose
  stability production had not scanned); the exact support class equals `en.efficiency_class` of the stratum (L1) for all.
* how the bits were decided: E by the exact support (852); C by interval alone 627, interval + structural zeros 150,
  exact rationals needed 75; S by interval 206, interval + structural 307, structural gain 1, exact rationals 338.
* production's 1e-9 tolerance decided a bit at 18 strategies (C at 16: no double margin below -1e-9 but some in
  [-1e-9, 0); S at 2: no double gain above 1e-9 but some in (0, 1e-9]). At all 18 production's bit is right: no exact
  margin is negative (no exact gain positive), and every co-player sent to exact arithmetic for the bit in question
  (its enclosure contained 0; 5 792 to 43 832 per strategy) is an exact tie, none on either side of 0.
* the deciding values: over the probe, the double values production tolerated were rounding noise of at most 1.1e-16
  (margins) and 2.2e-16 (gains), while every non-competitive strategy has an exact worst margin <= -1.7e-4 and every
  unstable strategy an exact best gain >= 5.6e-6. The tolerance 1e-9 lies in an empty gap from 2.2e-16 to 5.6e-6,
  more than six orders of magnitude above the largest rounding error and more than three below the smallest
  deciding violation (`tolerance_margins` in `results_probe/summary.json`).
* a first confusion table (per regime and wedge): `results_probe/confusion.txt` -- diagonal throughout. The shares
  there are indicative only (4 of 128 systematic chunks, inclusion probabilities scaled by 4/128).
* Mac (clang, 3 threads) and Cannon (GCC 8.5.0, 1 thread) outputs are identical, field by field, on the 27 probe
  strategies re-run locally (rows i mod 1024 = 0), including production's double diagnostics (not general: sid 3712
  differs by one ulp in pr.pp and pr.gmax, a compiler difference; every audit field agrees).

Cost (`src/project.py out_probe`): 1.39 s per strategy on one sapphire core (sum of `sec` 1185 s for 852); by production
atom: 000 0.37 s, 100 0.37, 011 0.38, 010 1.20, 110 1.51, 001 1.81, 111 2.75 (max 11.6 s; the cost is the number of
co-players sent to exact rationals, up to 49 000). Projection for the full sample of 27 261: 38 400 s = 10.7 core-hours,
of which the 124 remaining chunks 10.3; per chunk mean 300 s, sd 24 s; slowest of 128 chunks 369 s on average, 422 s at
the 99.9th percentile of 1000 bootstrap replicates. `--time=00:11:00` = 2 x the worst probe task (5:22), which is 1.6 x
the bootstrap 99.9% slowest chunk. Reserved at most 124 x 11 min = 22.7 core-hours.

Changes after the probe (the audit logic, and therefore every bit, is unchanged; the probe's four chunks are final):
`audit.py` records an exception as an unresolved strategy instead of ending the task, counts errors in the `.done` line,
and, when a production bit changes, makes sure the co-player named as responsible is one whose exact margin is negative
(gain positive) -- a fallback that never fired in the probe (no bit changed); `analyze.py` merges several output
directories, handles a partial set of chunks and error records, and adds ratio (Hajek) shares, per-bit flip
directions, the production shares of each group and the tolerance gaps above. The probe-era files are in
`src/probe_version/`.

Deviation from the plan: the sample (27 261 strategies) is larger than the "few thousand" suggested, because the probe
showed the whole sample costs about 11 core-hours; the larger sample makes the 95% bound on a changed share, when no
change is found, 0.06% of each regime's population (0.09% in W, N, E, 0.25% in S) instead of about 1%.

**Full run: JOBID_full = 48981982** (sapphire, 124 tasks = chunks 4-127, one core each), submitted 2026-09-27 after the
review below, ledger line written. `sacct`: all 124 tasks COMPLETED, exit 0:0, elapsed 0:04-0:06 per task (longest
6:15), first start 21:53:57, last end 22:00:52; 10.3 core-hours (CPUTimeRAW), TotalCPU 10.0 h, MaxRSS at most 49 MiB
(`PRFull/release/accounting`). All 124 `logs/full_*.err` are empty; every `chunk_<k>.done` of the 128 chunks reports all
its records and 0 errors. Outputs fetched to `out/`.

Results of the whole sample (`python3 src/analyze.py out_probe,out results`, about 8 s; `results/`, printed output
`results/analyze.out`; rerun by `PRFull/release/reproduce_certificates.sh`, step `evoaudit`):

* **27 261 of 27 261 audited**, 0 errors, 0 missing, 0 unresolved; production's classifier reproduced the dumped atom
  of all 27 261; the exact support class equals `en.efficiency_class` of the stratum (L1) for all; no validated 101.
* **2 changed atoms**, both `010` -> `000`, both at the stability bit: sid 10397 (F1, N = 1000, W, task 436) and sid
  15610 (F2, N = 100, S, task 146). Both are stability near-ties: their self-play limit puts a weight of 4.9e-11 and
  6.5e-10 on DD, and binary deviations (co-player 499 named) gain exactly that much, 4.874925678887806e-11 and
  6.459012207986248e-10 exactly, below production's tolerance of 1e-9; production's double values agree to 1.1e-16
  and it tolerated them. **No bit of efficiency or competitiveness changed** (95% upper bound on a changed-C share:
  0.060% at F1, 0.058% at F2).
* all 136 sampled W friendly rivals (50 F1, 86 F2) confirmed, and no other sampled W strategy is one.
* how the bits were decided: E by the exact support (27 261); C by interval 19 881, interval + structural zeros 4443,
  exact rationals 2937; S by interval 6637, interval + structural 10 725, structural 136, exact rationals 9763.
* production's 1e-9 tolerance decided C at 855 strategies and S at 61; exact arithmetic agrees with production at all
  of them except the 2 above. For C the tolerance lies in an empty gap (accepted double margins >= -1.1e-16, rejected
  strategies have exact minimum margins of -3.06e-7 or less); for S it does not (the two gains above are below 1e-9).
* the population shares (Horvitz-Thompson) and the 95% bounds per regime and wedge are in `results/confusion.txt` and
  `summary.json`; 37 434 core-s recorded by the strategies (1.37 s per strategy).

## Adversarial review (2026-09-27, before the full run)

Checked against wfc.c/m2f.c, the cert library and the dumps; nothing submitted to Cannon, nothing deleted there.

* Conventions: pi(q,p) = R w_CC + T w_CD + S w_DC (+ P w_DD, P = 0) in the focal-p frequencies (CD = p cooperates,
  q defects: q earns T); pi(p,p) = R w_CC + S w_CD + T w_DC of the exact self-play limit; margin sign(T-S)(w_DC - w_CD)
  (p's payoff minus q's is (T-S)(w_DC - w_CD)); wedges W (+1,C), S (-1,C), N (+1,A), E (-1,A), equal to the tag letter
  of games.dat at all 512 games; atom = 4E + 2S + C: the dumped atoms summed by replicate reproduce the a000..a111
  columns of every REP line of the paper's runs (PartnersRivalsM2f/check/wfc_a333 and wfc_f2, 5120 + 5120 replicates,
  0 mismatches; the kit's r_ files are byte-identical to those), and the REP gid/ipt/u/v match games.dat at every task.
* Design: pi_i recomputed independently from the dumps (max difference 0), the draw reproduced (27 261, same rows), the
  stored p equal to the dumps (spot check every 997th); the 95 % bounds are 0.0603 % / 0.0582 % per regime, at most
  0.092 % in W, N, E and 0.2450 % / 0.2366 % in S (README said 0.24 % for S; corrected to 0.25 %).
* Independent recomputation of 12 probe strategies with the library, not audit.py (E from the exact self-play
  frequencies; C by `ce.certify_point`; S in the q-FOCAL chain, `en.i_eval(q, q, p)` enclosures combined in exact
  rationals, `ex.lead_exact(q, p)` wherever they do not exclude 0, so the mirror identity is not used): sids 3712,
  19075, 21888, 3584, 16256, 13184, 2817, 12800, 8961, and in full exact rationals against all 65 536 q (margin from
  the (p,q) chain, pi(q,p) from the (q,p) chain) sids 4096, 12547, 18048.  All 12 atoms equal the audit's and
  production's; every count the audit resolved completely (negatives, zeros, gains, ties) is equal; the mirror identity
  holds exactly at all 3 x 65 536 q.  The current `src/audit.py` re-run locally on 5 probe records reproduces Cannon's
  records in every audit field.
* Wedge W: the paper's share of 111 there is production's exact count ex111 (wfc.c settle111), which tests only
  strategies with p0 = 1 and p15 = 0 at birth and uses the same double tests (W_CC > 1 - 1e-12, margin >= -1e-9) as the
  dumped atoms.  All 136 sampled W friendly rivals (50 F1, 86 F2) are in the sample with pi = 1, all have p0 = 1 and
  p15 = 0.  `analyze.py` now also reports `W_friendly_rivals` (production and validated E = C = 1 among the sampled W
  strategies, the flips either way, and any validated one outside the ex111 filter); probe: 6 and 6, no flips, none
  outside.  The ex111 count itself (strategies never sampled) is not audited.
