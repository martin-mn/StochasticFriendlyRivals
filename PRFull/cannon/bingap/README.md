# `bingap`: exhaustive census of the binary-pair margins (run on 2026-09-27; result at the end)

## Why this kit exists

The Methods (ms.tex:326), `notes/facts_m2.md:114` and M2f `Text/methods.tex:76` all say that the smallest non-zero |w_DC − w_CD| over binary memory-two pairs is 1/37. That is false.

- The pair p = 15204, q = 21477 gives exactly 1/963. Both `PartnersRivalsM2f/check/exact.py` and `check/cert/exact.py` agree.
- In a sample of 10⁵ random pairs, 637 fall below 1/37 (see `check/cert/tests.py`, section d).

This kit computes the true minimum over all 65536² pairs, with no sampling and no tolerance. It also computes the largest D\*, where D\* is the number of minimum-order spanning trees.

## Method

This is the integrality argument of `check/cert`.

1. For a binary pair every edge's leading coefficient is 1.
2. So A\*_k is the number of minimum-order spanning trees rooted at k, D\* = Σ A\*_k, and N = D\*·(w_DC − w_CD) is an integer.
3. The outward-rounded interval engine (`certeng.c`, a copy of `check/cert/certeng.c`) encloses D\* and the margin. Each integer is accepted only when its enclosure isolates it. Any pair where that fails is counted as `flagged` and printed to stderr. None was flagged locally.
4. The mirror symmetry gives m(q,p) = −m(p,q) and D\*(q,p) = D\*(p,q), so only q ≥ p is computed. That is 65536·65537/2 = 2 147 516 416 pairs. `--checksym` verifies the symmetry on every 97th q < p.
5. Task t of 64 handles every focal genome with p ≡ t (mod 64), which keeps the load balanced.

## Measured rate (local, native arm64, 1 core, `nice 10`)

1.70 µs per pair. Two slices were run: tasks 0 and 3 of 4096, 557 056 pairs each, 0.9 s each. Neither slice flagged a pair or broke the symmetry. Their smallest non-zero margin was 1/49 and their largest D\* was 624.

At that rate the full census is about 1.0 core-hour. The Cannon rate must come from the probe, not from rescaling this number.

## Sizing

| array | tasks | pairs per task | local-rate time per task |
|---|---|---|---|
| `bingap_probe.sbatch` (`test` partition, `--pmax 8192`) | 8 | ≈ 7.9 M | ≈ 13 s |
| `bingap.sbatch` (`shared`) | 64 | ≈ 33.6 M (see `tasks.txt`) | ≈ 57 s |

The full array is `--time=00:30:00` and 500 MB per task. Re-size it from the probe's `us_per_pair`.

## Protocol (Martin's CLAUDE.md)

1. Stage the kit:
   ```bash
   rsync -a PRFull/cannon/bingap/ cannon:PRFull/cannon/bingap/
   ssh cannon "cd PRFull/cannon/bingap && ./compile.sh"
   ```
2. Submit the probe and record its job ID in the same step:
   ```bash
   jid=$(ssh cannon "cd PRFull/cannon/bingap && sbatch --parsable bingap_probe.sbatch")
   echo "$jid" > JOBID_PROBE
   ssh cannon "printf '%s\tPRFull\tbingap\t%s\t8\tprobe, binary-pair margin census\n' $(date +%F) $jid >> ~/jobs.tsv"
   ```
3. Check the probe:
   ```bash
   sacct -j $jid --format=JobID,State,ExitCode,Elapsed -X
   ```
   Every task must be COMPLETED and `out/*.err` must be empty. Read `us_per_pair` from `out/task_*.txt`.
4. Submit the full array with `bingap.sbatch`, recording `JOBID` and a `jobs.tsv` line in the same step.
5. Verify the full run, then collect:
   ```bash
   sacct -j $(cat JOBID) --format=JobID,State,ExitCode -X --noheader | awk '$2 != "COMPLETED"'
   ```
   This must print nothing, and every `out/*.err` must be empty. Then run `python3 collect.py`.

## Output fields (one line per task)

- `pairs`, `zero`: exact zeros.
- `flagged`: enclosure did not isolate an integer.
- `symmetry_violations`
- `below_1/37`, `below_1e-9`: counts of non-zero margins under each bound.
- `min_nonzero N/D at p, q`
- `max_Dstar at p, q`
- `seconds`, `us_per_pair`

## What the result decides

- If the true minimum is at least 10⁻⁷, the double-precision binary census decisions at tolerance 10⁻⁹ are rigorous, and ms.tex:326 can quote the minimum.
- A priori the only bound is 2⁻³⁰ ≈ 9.3·10⁻¹⁰. The de Bruijn graph B(4,2) has 4¹³ arborescences per root, 2³⁰ in all, and that bound is *below* the 10⁻⁹ tolerance.

## Result (2026-09-27)

Probe 48916871 (4 tasks, test partition, 2.4 µs per pair on Cannon), then the full array 48917106 (64 tasks, all
COMPLETED, empty stderr, 2.1 core-hours). Over all 65536² binary memory-two pairs (2 147 516 416 computed with
q ≥ p, the rest by the mirror symmetry, verified on every 97th q < p):

- exact zeros 502 294 197; enclosures that failed to isolate an integer: 0; symmetry violations: 0;
- smallest non-zero |w_DC − w_CD| = **1/4140** ≈ 2.4 × 10⁻⁴ (p = 14002, q = 29490); 13 785 666 pairs lie in (0, 1/37);
- no non-zero margin below 10⁻⁹;
- largest number of minimum-order spanning trees D* = 36 864 (p = q = 13107).

So the binary census decisions at tolerance 10⁻⁹ are rigorous, and ms.tex's "1/37" is replaced by 1/4140. Outputs:
`out/`, `probe_out/`, `results.txt`.
