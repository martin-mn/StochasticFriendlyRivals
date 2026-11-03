# `cannon/cw`: the Wright–Fisher run over the continuous memory-one strategy space

    cw.x     the production program: sf1.f's process with s in [eps,1-eps]^4, the exact-alpha kernel,
             each sampled individual credited to the exact eps->0 atom (E,N,R) of its strategy
    cwp.x    the same at itend = 1e6, for the Cannon rate probe only
    kappa.x  the kernel's own mass on each atom and stratum per game: the reference for enrichment
    cwb.x    BINMODE, the audit build: sf1's 16 corners through cw's machinery, bit-identical to sf1.x

    N = 1000, m = 1, beta = 100, u = 1e-2, eps = 1e-4, itend = 1e7 (second half sampled), 10 replicates,
    the 512 sunflower games of DiskM1WF (Voronoi cells), cell ipt = task 3*ipt+3, array 3-1536:3

    run F2 (section 11): cw2.x / cwp2.x / cwb2.x, the same programs built from cw2.f = mkcw2.py(cw.f):
    N = 100, m = 1, beta = 3, u = 1e-4, eps = 1e-4, itend = 1e8, 10 replicates, the same 512 games

Figs 4 and 5 already have a Memory-1 row, but that row is the **binary** run: 16 strategies, DiskM1WF
run F1 (`DiskM1WF/Cannon/mw/sf1.f`). This kit runs the same process over **all** of memory one. Two changes make
that possible:

1. A strategy is a lineage that carries its own real 4-vector.
2. Every lineage is classified exactly when it dies, so each cell ends with 8 exact atom counts.

Fig 4(b,c) and Fig 5 are drawn from those counts. Fig 4(c) also needs the kernel's own atom masses, which
`kappa.x` computes. The figure scripts are `../../Figures/fig4c.py` and `fig5c.py`. They draw three rows:
Memory-1 binary, **Memory-1 continuous**, and Memory-2 binary.

**Status (2026-09-24): Cannon-ready. NOTHING HAS BEEN SUBMITTED.**
- Every validation below passes with the final sources.
- The complete chain has been run locally, end to end, at a short itend: run, pack, both figures.
- The Cannon rate is not yet known. The probe measures it; see Protocol.

**Run F2 (2026-09-25): the N = 100 run for SI Figs 4 and 5 -- Cannon-ready, NOTHING OF F2 SUBMITTED.** The F1
production array **48366968** is live while F2 is prepared: **nobody relinks cw.x (or cwb.x, cwp.x, kappa.x)
on Cannon while it runs** -- build only `bash compile.sh f2` and `bash compile.sh probe2` there. Section 11.

---

## 1. The process: sf1.f's, with these deviations and no others

`cw.f` is `ref_sf1.f` (a byte copy of `DiskM1WF/Cannon/mw/sf1.f`, sha256
`7d4649e312080253d1d8acd9019da54a5d5c1ef7ef22595796bd3a821063da7a`) with tagged edits. `diff ref_sf1.f cw.f` is
the audit, and `CHANGES.md` lists every hunk. `mt.f`, `mtb.f`, `payf2.f`, `payf3.f` and `games.dat` are byte
copies (`cmp`-identical to `DiskM1WF/Cannon/mw/`).

**Unchanged from sf1:**
- N = 1000, m = 1, beta = 100, u = 1e-2, itend = 1e7, init = itend/2, nrep = nper = 10, eps = 1e-4 (the ie = 3
  slot), and games.dat.
- The task decode.
- The interaction loop: each i draws j != i uniformly, then `f(i) += vl`, `f(j) += vh`.
- The synchronous update:
  - draw r;
  - if r < u, i becomes a mutant;
  - otherwise i draws a role model j uniformly from the whole population, including itself;
  - Fermi with `be = beta/(2m)`, `thr = 25` and the sign-split formula, then one draw against prob.
- The RNG: `mtfill` into `rb(8192)`, and `init_genrand` per replicate.
- pay and ef accounting, and `Emax = max(R,(S+T)/2)`.

| tag | deviation | why |
|---|---|---|
| [cw 1] | **Strategy storage.** Each individual holds a lineage slot. A slot carries `strat(4,ic)` (the effective vector), a serial `gid` that increases for the whole task, and a count. Dead slots are recycled after each update. Mutants are staged and placed after the recount. This is m1m2.F's class table, with `maxc = 2n+2` | a real-valued strategy cannot be an index into a 16-entry table |
| [cw 2] | **Payoffs.** An ordered-pair cache is keyed by `(gid(ia), gid(ib))` (m1m2.F's `pkey/pc1/pc2`). On a miss, `pay4` solves the 4-state pair chain by GTH state reduction (subtraction-free) and forms sf1's own two lines, `vl = rr*w0+ss*w1+tt*w2` and `vh = rr*w0+tt*w1+ss*w2`. The co-player acts on the swapped state | the payoff table becomes a cache. `pay4` agrees with `payf3` on all 256 corner pairs to 4.4e-16 |
| [cw 3] | **Mutation and initial population: the exact-alpha biased kernel** (section 2). The initial population is n independent kernel draws | the chosen kernel, which puts mass on every face of the cube |
| [cw 4] | **Sampling.** In the second half, lineage abundance `iab(ic) += icnt(ic)` (int64) per generation. The per-individual `ef` loop is sf1's | integer counts, so the invariant below is exact |
| [cw 5] | **Classification** of a lineage's abundance when it dies, and at the end of each replicate for every live lineage: atom, stratum, corner. The 16 corners are classified once per task and cached. **Asserted every replicate**: `sum(atomcnt) == sum(strcnt) == (itend/2)*n`, otherwise the run stops | the atoms are exact predicates of p, so they are applied once per lineage and never per individual per generation |
| [cw 6] | **Outputs** (section 6). sf1's file names, sf1's flush discipline | the h-file carries atoms, strata and corners in place of 16 strategies |
| [cw 7] | **`-DBINMODE`**: the 16 corners, sf1's single integer mutation draw, sf1's 4 fair coins, and sf1's `payf3` table read on a cache miss, all through the lineage, cache and abundance machinery of 1, 2 and 4 | the audit: BINMODE is **bit-identical to sf1** (V3) |
| [cw 8] | `-DITEND=<n>` overrides itend (tests and the probe build) | a short run with the production source |
| [cw 9] | `-DCLSREP=<k>`: the timing build, which repeats each `m1cls` call k times | measures the classifier's share of the time (section 9) |

**Seeds.** Continuum: `40000000 + (isl-1)*10 + irep`. This is disjoint from every other campaign (6M, 8M, 9M, 10M,
11M, 12M, 13M, 20–25M, 30–32M and 39M are taken). BINMODE uses sf1's `10000000`, so that it can be compared with `sf1.x`.
kappa: `41000000 + isl`, in its own process with its own MT19937.

## 2. The kernel

This is m1m2.F's `umut4` on the **exact-alpha** path (`-DALPHA`), ported to sf1's `rb/irb/mtfill` buffer. It
draws one uniform U per coordinate, with alpha = 1/3:

    uw = (1 - 2 eps)/(1 - 2 alpha),   ua = eps - alpha*uw,   s_k = clip(ua + uw*U, eps, 1-eps)

So each coordinate is independently:
- exactly `pmin = eps` (p_k = 0) with probability 1/3;
- exactly `pmax = 1-eps` (p_k = 1) with probability 1/3;
- uniform on the interior otherwise.

Each of the 81 strata therefore has mass exactly 1/81. The kappa check confirms this at every cell (section 8),
and 1-(2/3)^4 = 0.8025 of draws need an `m1cls` call (measured 0.8024).

**This deliberately does NOT reproduce the original `s.f` alpha quirk.** That code used the literal
`ua = -1, uw = 3`, which gives alpha = (1+eps)/3.

## 3. The classifier (`m1cls.c`, called once per lineage through `cwcls1`, `bind(C)`)

`m1cls(s, eps, u, v, atom, stratum, neartie)`. The game is (R,S,T,P) = (1, u, 1+v, 0), with u = cS and v = cT-1
taken exactly from games.dat.

**Preimage.** `s_k == eps` gives p_k = 0 (code 0). `s_k == 1.0-eps` gives p_k = 1 (code 1). Otherwise p_k is
(s_k-eps)/(1-2eps), interior (code 2). The stratum is `sum code_k 3^k`. Classification is the **exact eps->0
atom of p (convention A): limit first, weak inequality after**.

**Atom** = 4E + 2N + R, where E, N and R are:
- **E**: the closed forms of `Text/friendly.md` §1 = `check/m1lib.py`: `in_eff_cc` for u+v < 1, `in_eff_alt` for
  u+v > 1, `in_eff_line` on u+v = 1. They depend only on the stratum, so E is exact. The sign of u+v-1 is exact
  (TwoSum).
- **R**: `in_riv_plus` for T > S, `in_riv_minus` for T < S, and every strategy on T = S. The one non-stratum
  clause, p2+p3 <= 1, is decided exactly from s2+s3-1 by TwoSum. R is exact.
- **N**: Nash iff max over the 16 deterministic memory-one deviators sigma of [pi_lim(sigma,p) - pi_lim(p,p)] <= 0.
  This is the quantifier collapse of `Text/partners.md` §1 and `Text/nash.md`.
  - The limits come from the leading-order (coefficient, exponent) GTH of `../m1f_pairs.c`: `gth` (lines 36–80),
    `branch` (85–90), `swapo` (92), `freqs_cont` (96–114). `m1cls.c` documents the port.
  - Two deliberate changes to the port: the branch terms are computed once per strategy, and 1-p of an interior
    coordinate is formed from s directly. Before that fix, 1-p lost all relative accuracy within 1e-12 of p = 1;
    the worst margin error is now 1.8e-15.
  - Decision: tau = 1e-12·scale, where scale = max(1,|u|,|1+v|).
  - The deviators run in the order ALLD, ALLC, 1..14, with early exit only on a gap > 1e-9·scale. The early-exit
    result is therefore **identical** to the full test by construction, and that was checked at every test point.
- **neartie**: 2 if the margin is exactly 0.0 (a structural tie, decided by the weak inequality); 1 if
  0 < |margin| <= 1e-9·scale (rounding-sensitive); 0 otherwise. `cw.f` counts 1 and 2 separately (w-file).

**Validation** (`python3 test_cls.py`, 41 s; rerun with the final sources on 2026-09-24):

| test | what | result |
|---|---|---|
| V1 (C) | 16 corners x 512 games against `DiskM1WF/Data/m1_atoms_512.npz` (exact rationals) | **0 mismatches in 8192**. Atom totals [5072, 1310, 409, 486, 291, 0, 372, 252] (index order 0..7) = npz column sums. Early exit = full test at all 8192. neartie at the corners: 6673 of type 0, 0 of type 1, 1519 of type 2 |
| **V1 (linked binary)** | `cwc.x` (cw.f at itend = 2) run at **all 512 cells**. Its h-file line `# corner atoms`, which is the corner cache built by cw.f's own `cwcls1 -> m1cls` call, compared with the npz | **0 mismatches in 8192** |
| V2 | 3280 points (3120 non-corner and 160 corners), all 81 strata, 29 games (681 W, 678 S, 678 N, 678 E, 565 on boundary lines), against an exact-Fraction oracle from `m1lib` at eps = 1e-100 | E 0, R 0, stratum 0 mismatches; early exit vs full: 0. **N: 30 mismatches, all explained.** Each is a C "Nash" with an exact margin in (0, 1e-12], every one neartie-flagged, all among the 1318 adversarial draws within 1e-6..1e-14 of a face. The kernel, dyadic and p2+p3=1 draws have 0 mismatches. Max \|C margin - exact\| = 1.78e-15. Closed forms vs definitions at 1e-100: E 0, R 0 |
| V2 constructed ties | 14 points on the ZD sheet or on the diagonal x*: exact ties that the double preimage can only approximate | all 14 flagged; they differ from the oracle by at most 3.4e-16 |
| Fortran link | gfortran with {clang, gcc} objects, implicit external (`m1cls_`) and `bind(C)` (`m1cls`) | all 4 builds agree with C. ifx was tested only through cw.f's `bind(C)` call on a Cannon login node (V3 ifx) |

**What the N tolerance costs in production.** A kernel draw lands within 1e-12 of a face with probability about
1e-11. That is roughly 0.005 misclassified lineages per task, out of 4e8 classifications. In the 512-cell local
run, **0 of 82,230,872** `m1cls` calls returned neartie = 1.

## 4. Enrichment, and the grey rule

- **Share** s(a): the time-averaged fraction of the population in atom a over the sampled half of the 10
  replicates. It is computed as `atomcnt(a) / total` from the h-file's exact counts.
- **Kernel mass** kappa(a): the fraction of M = 1e7 kernel draws at that game lying in atom a (`kappa.x`, seed
  41000000+isl). It is the continuum's "share of the strategy space". In the binary rows, that role is played by
  |a|/16 and |a|/65536.
- **Enrichment** e(a) = s(a)/kappa(a). Fig 4(c) shows argmax_a e(a) in the 7 atom colours, the same categorical
  key as the original, so no enrichment colour range enters the figure.
  - Atoms with kappa = 0 and s = 0 are excluded, as the original excludes |a| = 0.
  - If kappa_hat = 0 but s > 0, the true kappa is only below the resolution 1/M. It is set to 3/M, the 95 % upper
    bound, so that e is a conservative lower bound.
  - The ceiling is 1/kappa, not 16. `fig4c.py` prints the winners' enrichment range, their smallest kappa (in
    kernel draws), and how many winners are not resolved from the runner-up by kappa_hat's own noise.
- **Fig 5** shows the share holding each property: efficient 1** (100, 110, 111), Nash \*1\* (010, 110, 011, 111),
  competitive \*\*1 (001, 011, 111). A cell is grey only where the kernel's mass on the property is 0 **and** the
  population holds none of it.
- **Atom 101 is empty** at every game in theory (efficient + rival implies Nash). `pack.py` refuses a cell whose
  101 count is non-zero, and `cwdata.py` refuses a pack with non-zero s101 or k101.

## 5. Protocol (Cannon)

Copy these files to `<remote-kit>` on Cannon: `cw.f ref_sf1.f mt.f mtb.f payf2.f payf3.f games.dat m1cls.c m1cls.h
kappa.c compile.sh task.sh probe.slurm prod.slurm kappa.slurm submit.sh verify.sh pack.py check_v4.py`. Do not
copy the local outputs `run_local/ kappa_local/ packed_local/` or the Mac binaries. Then, on a login node:

```bash
bash compile.sh                  # cw.x cwb.x kappa.x, with ifx -O3 -xHost -fp-model precise / icx
bash compile.sh probe            # cwp.x (itend 1e6)
ITS=2000 bash compile.sh test    # optional: re-audit BINMODE vs sf1 under ifx (V3); run in a scratch dir
bash submit.sh probe             # 5 tasks, partition test  -> probe/, JOBID_probe, IDX_probe, ~/jobs.tsv
bash submit.sh kappa             # 1 task (loops over the 512 cells), shared -> kappa/, JOBID_kappa
bash verify.sh probe             # sacct per task, stderr, and the RATE -> "set prod.slurm -t to about ..."
#   edit prod.slurm's -t (and --mem from maxRSS) from the probe, THEN
bash submit.sh prod              # 512 tasks, --array=3-1536:3, partition sapphire -> run1/
bash verify.sh prod              # repeat until every task is COMPLETED and 512 of 512 cells are complete
bash verify.sh kappa
python3 pack.py run1 kappa packed    # -> packed/cw1_e4.tsv, cw1_e4.rep.tsv, kappa_e4.tsv (no numpy needed)
```

Then copy `<remote-kit>/packed/` to this directory's `packed/` on the Mac and run
`python3 ../../Figures/fig4c.py` and `fig5c.py`. They write `Figures/Fig4c.{pdf,png}` and `Fig5c.{pdf,png}`.

- **`submit.sh`** does the following in one step:
  - stages games.dat by temp file and rename;
  - refuses a binary not linked with the real classifier (it checks the `m1cls_detail` symbol);
  - refuses a second submission while that tag's array is still queued;
  - cross-checks the ledger's task count against the script's own `--array=` line;
  - reports the user's queued/running count before the 512-task array, and refuses above 200 without
    `QUEUE_OK=1`;
  - submits with `sbatch --parsable` into `JOBID_<tag>`, writes `IDX_<tag>`, and appends
    `date  PartnersRivalsM1f  cw  jobid  ntasks  [tag] note` to `~/jobs.tsv`.

  Job names are `PartnersRivalsM1f.cw{probe,prod,kappa}`, with `--comment=PartnersRivalsM1f`.
- **`verify.sh <tag>`** decides completion **per task from sacct**, never from `squeue -u`. It lists every task not
  COMPLETED, any non-empty stderr log, the core-hours, and the cells complete within `IDX_<tag>`. For the probe it
  also prints the measured s per 1e4 generations, the implied production wall (x10), a suggested `-t` (3x the
  worst), and maxRSS.
- **`task.sh <isl> <rundir> <exe>`** runs one cell idempotently.
  - A cell is **complete** when the summary has 10 lines **and** the h-file has its `# total` line and exactly 105
    count lines. cw.f writes the h-file tail in one burst at the end, so 10 summary lines alone do not mean done
    (prod1.slurm's lesson, kept).
  - It skips complete cells, reruns incomplete ones from scratch, and exits non-zero if the products are
    incomplete afterwards, so sacct shows FAILED.
- **Why the probe is a separate build (`cwp.x`, `-DITEND=1000000`) in its own directory `probe/`.**
  - `cw.x` stays exactly the production program, and no argument can shorten a production task.
  - A probe cell can never pass for a production cell: `pack.py` refuses itend != 1e7.
  - The cost per generation does not depend on itend, because init = itend/2 in both and the lineage population
    equilibrates within ~1e3 generations. So production wall is ~10x the probe wall.
  - The five cells are mw probe1's (isl 3, 1527, 1530, 1533, 1536: centre PD, far PD, HA, SD, SH), so the rate can
    be set beside the binary run's measured 42–69 min at the same games.
- **Why partition sapphire for production** (mw used shared): the probe runs on `test`, which is Sapphire Rapids.
  Only production on Sapphire Rapids makes the probe's rate a measurement of the production hardware; shared is
  Cascade Lake. Code built with `-xHost` on a login node also runs safely on the newer SPR cores.
- **kappa is one task, not 512**: a cell takes ~2 s. The task writes each `kappa/k<isl>` by temp file and rename,
  skips finished cells on a rerun, and was tested locally (a removed cell is recomputed: "1 computed, 511 already
  present").

**`pack.py`** checks every cell, rather than assuming:
- **the task decode**: the h-file's `# cell:` isl/gid/tag against games.dat row (isl-1)/3+1, and kappa's isl/gid
  likewise. A deliberately swapped cell was caught: "h-file says isl 6 gid 2001 V_p001_HAW, games.dat row 3 is gid
  2002 V_p002_SDW -- TASK DECODE";
- (u,v) = (cS, cT-1) = the display columns;
- the run parameters (ie 3, eps, umut, beta; nper 10 and itend 1e7 unless `--dev`);
- 8/81/16 count lines, and atoms = strata = total = nper·(itend/2)·N;
- corners = corner strata, and 101 = 0;
- 10 replicate lines, each naming its gid, with shares summing to 1 and Emax matching;
- ef = mean pay / Emax against the mean of the replicate efs (1e-9);
- kappa: 98 fields, both groups summing to 1, k101 = 0, and M = 1e7 unless `--dev`.

The formats are exactly those that `Figures/cwdata.py` defines and reads.

## 6. Outputs per task (sf1's file names)

- **`<isl>`**: one line per replicate, flushed: `gid ie u v Emax pay ef meanlive solves_per_gen share0..share7`.
  The atom index is `4E+2N+R` (6 = `110`, 7 = `111`, 5 = `101` is empty).
- **`w<isl>`**: two header lines, then one line per replicate:
  `# isl irun ncls ncall nneartie neartie_abund nhw cache_hit_rate meanlive solves_per_gen gcount nexacttie exacttie_abund`.
- **`h<isl>`**:
  - sf1's header (`# cell`, `# game`, `# run`, `# ic`) and `# corner atoms (code 0..15)`;
  - at the end, in one burst: `# atoms 8` with 8 `index count` lines, `# strata 81` with 81, `# corners 16` with
    16 (sf1's bit code), then `# total N`, where N = 10·(itend/2)·n = 5e10 in production.
- **`kappa/k<isl>`**: one line of 98 fields: `isl gid ie eps u v M`, 8 atom fractions (index 4E+2N+R),
  81 stratum fractions, and the neartie-1 and neartie-2 fractions.
- BINMODE writes sf1's three files exactly, with diagnostics on stdout.

## 7. Build

`compile.sh` targets are `cw cwb cwp kappa` (the Cannon set; with no argument it builds `cw cwb kappa`), `probe`
(= cwp), `test` (= `sf1t cwt cwbt` at ITS = 20000), `cwc` (itend 2, the V1 re-check) and `cwr` (`CLSREP` timing).
- **On Cannon:** `module load intel/25.2.1-fasrc01`, then `ifx -O3 -xHost -fp-model precise -fpp` and
  `icx -O2 -xHost -fp-model precise`.
- **Locally:** `gfortran -O2 -ffp-contract=off -fno-fast-math -cpp` and `gcc -O2 -ffp-contract=off`. m1cls.c must
  never be built with fast-math, because its exact sums depend on strict IEEE arithmetic.
- The classifier is `m1cls.c`. The old test stub is quarantined in `quarantine/m1cls_stub.c`. `CLS=<file>`
  overrides the classifier for tests, and **a stub is refused for the cw, cwb, cwp and kappa targets**.
- Do not relink a binary whose array is live: the compiler writes it in place.

## 8. Validation of the simulator, with the final sources (2026-09-24, gfortran 14 on the Mac)

**V3: BINMODE is bit-identical to sf1**, now with the real classifier linked. `sf1t.x` (ref_sf1.f with only its
`itv` line changed to 20000; compile.sh asserts that one-line diff) and `cwbt.x` (`-DBINMODE -DITEND=20000`), both
built with the same flags, were run at isl 3, 768 and 1536. `cmp` of `<isl>`, `h<isl>` and `w<isl>` gives **9 of 9
identical** (1570 B, 806 B and 14331 B per cell). The simulator builder's earlier ifx run on a Cannon login node
(itend 2e3, isl 3 and 1536) gave 6 of 6.

**V4: continuum, `cwt.x` (itend 2e4, 10 replicates), real classifier**. The second group of cells is the PD, SH
and SD set used for timing:

| isl | game | h total = atoms = strata | max\|sum shares - 1\| | s101 | high water / maxc | classified | m1cls calls | neartie 1 | neartie 2 | kappa (M = 1e6): sum-1, k101 |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | V_p000_PDW | 1e8 = 1e8 = 1e8 | 0 | 0 | 1000 / 2002 | 1,001,649 | 804,038 | 0 | 193,561 | -1.1e-16, 0 |
| 45 | V_p014_SHW | 1e8 | 0 | 0 | 1000 / 2002 | 999,239 | 801,188 | 0 | 349,374 | 0, 0 |
| 9 | V_p002_SDW | 1e8 | 0 | 0 | 1000 / 2002 | 1,000,642 | 802,597 | 0 | 140,302 | 0, 0 |
| 768 | | 1e8 | 0 | 0 | 1000 / 2002 | 999,280 | 801,797 | 0 | 209,968 | |
| 1536 | | 1e8 | 0 | 0 | 1000 / 2002 | 999,947 | 802,629 | 0 | 480,388 | |

- At every cell the per-replicate invariant held (it is asserted), and the corner counts equal the corner strata.
  `python3 check_v4.py <dir> <isl>...` reproduces the table.
- At the PD centre (isl 3) the population sits 96.4–97.3 % in atom 110 (efficient Nash non-rivals, the
  partners), replicate by replicate.
- The high-water mark of 1000 is the initial population. It is never exceeded, so `maxc` has a 2x margin.

**The whole chain, locally, at all 512 cells** (a SHORT LOCAL RUN, not a result):
- `cwt.x` built at ITS = 4000, 10 replicates per cell, through `task.sh`, 6 processes at a time: 607 CPU-s,
  1 min 44 s wall.
- `kappa.x` at M = 1e6 for all 512 cells: 128 CPU-s.
- `check_v4.py`: 512 of 512 PASS.
- `pack.py run_local kappa_local packed_local --dev`: "all 512 cells complete and consistent". In production mode
  the same files are refused ("itend 4000, production is 10 / 10000000").
- Over the 5120 replicates:
  - 102,491,136 lineages classified and 82,230,872 m1cls calls;
  - **0 neartie 1**; 19,466,049 exact ties (neartie 2) carrying 71 % of the sampled abundance;
  - high water 1000, cache hit rate 0.968, 29.6 live lineages on average, 32.3 pair solves per generation.
- In the pack: max\|sum s - 1\| = 2.2e-16, s101 = 0, k101 = 0, and the strata are within 4.2 sigma of 1/81 at
  every cell.
- `fig4c.py --short` and `fig5c.py --short` then render `Figures/Fig4c_localshort.{pdf,png}` and
  `Fig5c_localshort.{pdf,png}`. They are watermarked SHORT LOCAL RUN, and their footer states itend = 4000 and
  M = 1e6.

## 9. Cost: LOCAL measurements only (the Cannon probe measures the real rate)

Single thread, gfortran -O2, Martin's M2 Pro (load average 4–7 from other jobs). Each figure is the wall time of
`cwt.x` (10 x 2e4 = 2e5 generations), measured twice:

| cell | wall (2e5 gen) | s per 1e4 generations | cwr.x (CLSREP = 11) | m1cls share | ns per m1cls call |
|---|---|---|---|---|---|
| isl 3, PD | 5.68, 5.68 s | 0.284 | 7.60, 7.37 s | 3.2 % | 225 |
| isl 45, SH | 5.81, 5.80 s | 0.290 | 8.20, 8.19 s | 4.1 % | 298 |
| isl 9, SD | 5.67, 5.84 s | 0.288 | 7.23, 7.21 s | 2.5 % | 183 |

The m1cls share is computed as (t_cwr - t_cwt)/10 / t_cwt, and `cwr.x` output equals `cwt.x` output byte for byte.

- **Extrapolated LOCAL estimate for one production task** (10 x 1e7 = 1e8 generations): about **2850–2900 s, or
  48 min**, of which m1cls is 2.5–4 % (about 4e8 calls per task).
- For comparison, sf1.x (the binary run) measured 0.228 s per 1e4 generations locally, and 42–69 min per task on
  Cannon.
- Peak RSS locally is 66 MB (the pair cache is at most 96 MB), so `--mem=400`. kappa.x uses 1.4 MB and takes 2.07 s
  at M = 1e7.
- `prod.slurm` has `-t 360` as a **provisional** value, 7.5x the local estimate. Set it from `verify.sh probe`
  before submitting production.

## 10. Files

| file | role |
|---|---|
| `cw.f` | the simulator (continuum, and BINMODE for the audit) |
| `ref_sf1.f` `mt.f` `mtb.f` `payf2.f` `payf3.f` `games.dat` | sf1's, byte copies |
| `m1cls.c` `m1cls.h` | the classifier |
| `kappa.c` | the kernel-measure program |
| `compile.sh` | builds everything |
| `task.sh` | one cell, idempotent (local and Slurm) |
| `probe.slurm` `prod.slurm` `kappa.slurm` | the three Slurm scripts |
| `submit.sh` `verify.sh` | submission (JOBID, IDX, ledger) and per-task completion from sacct |
| `pack.py` | run1/ + kappa/ -> packed/ (`--dev` for a short local run) |
| `check_v4.py` | the V4 invariant checks on any output directory |
| `test_cls.py` `test_cls.c` `test_fcall.f` `test_fbind.f90` | classifier validation (V1, V2, Fortran link, kappa) |
| `CHANGES.md` | every hunk of `diff ref_sf1.f cw.f`, and the integration revision |
| `quarantine/m1cls_stub.c` | the old test stub (atom = 0); refused for production targets |
| `run_local/` `kappa_local/` `packed_local/` | the short local run (itend 4000, M 1e6): a test of the chain, not data |
| `mkcw2.py` `cw2.f` | run F2: cw2.f is GENERATED from cw.f by mkcw2.py (asserted substitutions); never edit cw2.f |
| `ref_sf2.f` | byte copy of `DiskM1WF/Cannon/mw/sf2.f` (sha256 `e227b4c8...2061edf`), for the F2 BINMODE audit |
| `probe2.slurm` `prod2.slurm` | the F2 Slurm scripts (test / sapphire), calling the unchanged task.sh |
| `run_local2/` `packed_local2/` | the short local F2 chain (itend 1e4, kappa_local reused): a test, not data |

---

## 11. Run F2: N = 100, for SI Figs 4 and 5

`PartnersRivals1/Revised2/figures/SI/SIFig4.pdf`
and `SIFig5.pdf` are the N = 100 versions of Figures 4 and 5. Their Memory-1 binary row is DiskM1WF's run F2
(`DiskM1WF/Cannon/mw/sf2.f`, pack `DiskM1WF/Data/mw2_e4.tsv`). F2 here is that process over all of memory one,
exactly as F1 is sf1's.

### What differs from F1, and nothing else

`cw2.f` is **generated** from `cw.f` by `mkcw2.py`: an asserted list of 13 substitutions, each of which must match
exactly once (mw's `mksf.py` discipline). mkcw2.py first asserts cw.f's sha256 (`37ec579b...c729a2ebb`, the F1
production source). It then asserts the diff: 14 changed lines, 7 of them code and 7 comment, and the line count
unchanged. **`diff cw.f cw2.f` is the audit.** The code lines are:

| line | cw.f (F1) | cw2.f (F2) | = sf2.f vs sf1.f |
|---|---|---|---|
| 47 | `#define ITEND 10000000` | `#define ITEND 100000000` | `itv` 1e7 -> 1e8 (cw.f's itv default is the ITEND macro) |
| 56 | `parameter (n=1000,m=1)` | `parameter (n=100,m=1)` | n |
| 66 | `parameter (iseedb=10000000)` (BINMODE) | `parameter (iseedb=11000000)` | iseedb: sf2's own, so BINMODE-F2 compares with sf2 |
| 68 | `parameter (iseedb=40000000)` (continuum) | `parameter (iseedb=42000000)` | a new disjoint base |
| 94 | `data uval /1.d-02/` | `data uval /1.d-04/` | uval |
| 95 | `data bval /100.d+00/` | `data bval /3.d+00/` | bval |
| 240 | w-file header string `(maxc 2002)` | `(maxc 202)` | the class-table size stated in the output (maxc = 2n+2) |

The 7 comment lines are those that name F1's parameters or F1's comparison program (sf1.x -> sf2.x, sf1's 1e7 ->
sf2's 1e8, the seed bases). eps, nrep = nper = 10, games.dat, the kernel, the classifier, the task decode, the
generation loop, the RNG discipline and every output format are F1's. So is the file layout: the h-file names
`umut 0.100E-03  beta 3.0` and `itend 100000000`, and pack.py checks them.

- `bash compile.sh` refuses to build any F2 target unless `python3 mkcw2.py --check` passes, i.e. unless cw2.f is
  byte for byte mkcw2.py's output of the present cw.f.
- `submit.sh` refuses probe2/prod2 unless the same check passes **and** the binary's strings name `maxc 202)`.
  cw.x names `maxc 2002)`, so an F1 binary cannot be submitted as F2.

**Seeds.**
- Continuum F2: `42000000 + (isl-1)*10 + irep`. This is disjoint from F1's 40M, kappa's 41M, and the 6M, 8M, 9M,
  10M, 11M, 12M, 13M, 20–25M, 30–32M and 39M bases of the other campaigns.
- BINMODE-F2: sf2's `11000000`.
- **kappa is shared with F1.** The kernel's atom masses depend on the game and eps, not on N, beta or u. F2
  therefore uses the same `kappa/` directory and the same `kappa_e4.tsv`, and there is no F2 kappa job.

**Numerical range.** itend = 1e8 fits int32. Every product that could overflow is already int64 in cw.f:
- `nlook` = n·m·itend = 1e10;
- the invariant `(itend/2)·n` = 5e9 per replicate;
- the h total = 5e10, the same as F1.

gcount reaches ~1e7 per task, far below the asserted 2^31.

### Validation (2026-09-25, gfortran 14 on the Mac, real classifier)

**Build.**
- `ref_sf2.f` is `cmp`-identical to `DiskM1WF/Cannon/mw/sf2.f`, sha256
  `e227b4c8472ace1d7727481632a96a025c3d5d045bef402b6c849e0cf2061edf`.
- `diff ref_sf1.f ref_sf2.f` is sf2's 5 parameter lines and 2 comment lines.
- compile.sh gained the targets `cw2 cwb2 cwp2 cwt2 cwbt2 sf2t` and the groups `f2` (= cw2 cwb2), `probe2`
  (= cwp2, ITP2 = 1e7) and `test2` (= sf2t cwt2 cwbt2). The no-argument set is still `cw cwb kappa`.
- `diff` of the old and new compile.sh removes exactly one line: the stub-refusal `case`, which now also lists
  cw2, cwb2 and cwp2. Every existing build line is unchanged.
- All 9 existing targets (cw cwb kappa cwp cwt cwbt cwc cwr sf1t) were rebuilt before and after the edit. **All 9
  binaries are byte-identical** (sha256 compared). cw.x (`7fd8c9e0...eb16af`) and cwb.x (`dace4ae1...5c148534`)
  are also identical to the binaries that were in the kit before this session.

**V3-F2: BINMODE is bit-identical to sf2.**
- Test builds: `sf2t.x` is ref_sf2.f with only its `itv` line changed to 200000; compile.sh asserts that one-line
  diff. `cwbt2.x` is cw2.f with `-DBINMODE -DITEND=200000`. Both use the same flags.
- They were run at isl 3 (PD centre), 9 (SD), 45 (SH) and 1536 (SH, the farthest). `cmp` of `<isl>`, `h<isl>` and
  `w<isl>`: **12 of 12 identical** (1570, 806 and 14331 B per cell).
- Each run takes 4.8–6.5 s.
- WSLS (code 9, `CDDC`) is the rank-1 strategy at the PD centre in all 10 replicates (share 0.9992–0.9997), as
  in mw's F2.
- BINMODE's stdout gives `sum(atomcnt) = 10000000 = (itend/2)·n` for every replicate.

**V4-F2: the continuum.** `cwt2.x` (itend 2e5, 10 replicates) was run at 3 cells; `python3 check_v4.py --n 100 . 3 45 9` gives PASS at all three:

| isl | game | h total = atoms = strata | max\|sum shares - 1\| | s101 | nhw / maxc | classified | m1cls calls | neartie 1 | exact ties | meanlive | solves / gen | wall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | V_p000_PDW | 1e8 = 1e8 = 1e8 | 0 | 0 | 100 / 202 | 10,064 | 8,085 | 0 | 2,018 | 1.046 | 0.0214 | 5.73 s |
| 45 | V_p014_SHW | 1e8 | 1.1e-16 | 0 | 100 / 202 | 10,037 | 8,101 | 0 | 3,432 | 1.038 | 0.0209 | 5.41 s |
| 9 | V_p002_SDW | 1e8 | 0 | 0 | 100 / 202 | 9,868 | 7,932 | 0 | 1,337 | 1.044 | 0.0212 | 5.47 s |

- The per-replicate invariant `sum(atomcnt) == (itend/2)·n` is asserted by the program. It held in every
  replicate, and the corner counts equal the corner strata.
- **Mean live lineages ≈ 1.04.** At N·u = 0.01 the population is monomorphic almost always.
- Rates per generation:
  - pair solves: 0.021, with a cache hit rate of 0.99979;
  - lineages classified: 0.0050 per generation (0.010 per sampled generation, about one per mutant);
  - m1cls calls: 0.0040 per generation.
- The class table never goes above the initial population (100 of 202).
- At the PD centre, 84.9 % of the sampled abundance (10 replicates pooled) sits in atom 110, the partners. The
  per-replicate share ranges from 9.8 % to 99.98 %: at itend 2e5 a replicate sees only ~2000 mutants.

**The whole F2 chain, locally, at all 512 cells** (a SHORT LOCAL RUN, not a result):
- `cwt2.x` at ITS = 1e4, run through the unchanged `task.sh` into `run_local2/`, 6 processes at a time: 153 CPU-s,
  27 s wall.
- `check_v4.py --n 100`: **512 of 512 PASS**.
- Over the 5120 replicates:
  - 261,593 lineages classified and 210,105 m1cls calls;
  - **0 neartie 1**; exact ties carry 71.5 % of the abundance;
  - high water 100, meanlive 1.16 (1.01–3.10), cache hit rate 0.9992.
- `pack.py run_local2 kappa_local packed_local2 --dev --run=F2`: "all 512 cells complete and consistent". It writes
  `cw2_e4.tsv`, `cw2_e4.rep.tsv`, `cw2_e4.diag.tsv` and `kappa_e4.tsv`. max\|sum s - 1\| = 2.2e-16, s101 = 0,
  k101 = 0. The run line is `# N = 100, beta = 3, u = 1e-4, eps = 1e-4, itend = 10000, 10 replicates,
  iseedb = 42000000`.
- Refusals tested:
  - in production mode (no `--dev`) the same run is refused: "itend 10000, production is 10 / 100000000";
  - F1 cells under `--run=F2`, and F2 cells under F1, are refused on umut and beta;
  - a truncated h-file, a swapped pair of cells (TASK DECODE, both directions), a 101 count of 1 and a near-tie
    share of 2e-5 are each refused;
  - a `kappa_e4.tsv` already in `<out>` that comes from a different kappa directory is refused. The shared F1
    file, with identical data lines, is accepted and left untouched.

**pack.py F1 is unchanged.** `pack.py run_local kappa_local <dir> --dev` (and with `--run=F1`) was compared with the
pre-F2 pack.py's output: all four files are `cmp`-identical (cw1_e4.tsv 97114 B, cw1_e4.rep.tsv 859854 B,
kappa_e4.tsv 397188 B, cw1_e4.diag.tsv 19879 B). The first three are also identical to `packed_local/`.

### Cost: LOCAL only (probe2 measures the Cannon rate)

- cwt2.x runs at **0.27–0.29 s per 1e5 generations** (isl 3, 45, 9). sf2t.x, the binary process, runs at
  0.24–0.26 s at the same three cells (0.31 s at isl 1536).
- At N = 100 with ~1 live lineage, the continuum machinery costs almost nothing. The per-individual loop sets the
  cost, as it does in sf2.
- **LOCAL estimate for one full F2 task** (10 × 1e8 = 1e9 generations): **~2700–2900 s = 45–48 min.**
- mw's binary F2 took 56.6–71.7 min per task on Cannon (shared, Cascade Lake; 1.43x its local rate) and 526 core-h
  for the 512 cells.
- Expect roughly 400–650 core-h for F2 here. **This is not a measurement: take it from `verify.sh probe2`.**
- Memory: maxc = 202, so the pair cache is 3 × 202² × 8 B = 1 MB. `--mem=400` as prod.slurm.

### Protocol (Cannon; the main session does all cluster work)

**THE RULE WHILE 48366968 (F1 production) IS LIVE: nobody relinks cw.x on Cannon.**
- The compiler rewrites the executable in place under the running tasks.
- Never run plain `bash compile.sh`, nor the targets `cw`, `cwb`, `cwp`, `kappa`, `probe` or `test`, in the remote
  kit until `verify.sh prod` shows every F1 task finished.
- The F2 targets write only `cw2.x`, `cwb2.x`, `cwp2.x`, `cwt2.x`, `cwbt2.x` and `sf2t.x`, plus `build/m1cls.o`,
  which a running binary never reads.

Copy these files to the same `<remote-kit>`:

    cw2.f mkcw2.py ref_sf2.f compile.sh probe2.slurm prod2.slurm submit.sh verify.sh pack.py check_v4.py

Do not copy cw.f, task.sh, prod.slurm, probe.slurm, kappa.slurm, games.dat or any of the other untouched files:
they are byte-identical to the remote copies and are in use by the live array. Leave out the local outputs too.
submit.sh, verify.sh and pack.py are not read by running tasks, so replacing them is safe. Then, on a login node:

```bash
python3 mkcw2.py --check            # cw2.f = mkcw2.py(cw.f), against the REMOTE cw.f (sha asserted)
bash compile.sh f2                  # cw2.x cwb2.x   (ifx -O3 -xHost -fp-model precise)
bash compile.sh probe2              # cwp2.x (itend 1e7)
# optional, in a scratch copy: ITS=20000 bash compile.sh test2 -> V3-F2 under ifx (cmp sf2t vs cwbt2 outputs)
bash submit.sh probe2               # 5 tasks, test -> probe2/, JOBID_probe2, IDX_probe2, ~/jobs.tsv
bash verify.sh probe2               # sacct per task, stderr, and the RATE -> "set prod2.slurm -t to about ..."
#   edit prod2.slurm's -t (provisional 240) and --mem from the probe, THEN
QUEUE_OK=1 bash submit.sh prod2     # 512 tasks, sapphire -> run2/, JOBID_prod2 (NOT JOBID), IDX_prod2
bash verify.sh prod2                # until every task is COMPLETED and 512 of 512 cells are complete
python3 pack.py run2 kappa packed --run=F2   # -> packed/cw2_e4.tsv, .rep.tsv, .diag.tsv; kappa_e4.tsv shared
```

Before `QUEUE_OK=1`, report the queued count that submit.sh prints, with the ledger breakdown (the global
protocol). The F1 array alone is up to 512 tasks.

- **JOBID files.** `JOBID_probe2` and `JOBID_prod2` are F2's own. The plain `JOBID` stays F1's campaign id,
  48366968: submit.sh writes `JOBID` only for tag `prod`. The RESUBMIT and QUEUE_OK guards apply unchanged.
- **Ledger notes** say `[probe2]` / `[prod2] ... run F2 ... N=100 beta=3 u=1e-4 itend=1e8 eps=1e-4 ... seeds
  42000000+`.
- **Logs.** `log/{o,e}_f2probe_%a.txt` and `log/{o,e}_f2prod_%a.txt`. The names are chosen so that verify.sh
  probe's glob `e_probe*` and prod's `e_prod*` never match F2 logs.
- **Partitions.** probe2 runs on **test** (Sapphire Rapids, the five cells of probe.slurm and mw probe2). prod2
  runs on **sapphire**, as F1 does: the login node builds `-xHost` Sapphire Rapids code, so **never run it on
  shared**.
- **`verify.sh probe2`.** Its rate line uses F2's generation counts: 1e8 generations per probe task, so s/1e4 per
  1e4 generations. The production task is 10x the probe task, and the suggested `-t` is for prod2.slurm.
- **Dry run of submit.sh and verify.sh**, in a scratch copy of the kit with fake sbatch, squeue and sacct, and
  HOME in scratch:
  - probe2 submits and writes JOBID_probe2, IDX_probe2 (5) and a ledger line;
  - a second probe2 is refused, both when not queued (RESUBMIT guard) and when live;
  - prod2 with 300 queued is refused without QUEUE_OK;
  - prod2 with QUEUE_OK submits and writes JOBID_prod2 and IDX_prod2 (512); RESUBMIT=1 keeps JOBID_prod2.prev;
  - the plain JOBID stayed 48366968;
  - cw.x copied over cw2.x is refused;
  - F1's verify.sh output (probe, prod, kappa) and the F1 submit refusals are byte-identical to the pre-F2
    scripts'.

The F2 figures (the N = 100 counterparts of fig4c.py and fig5c.py, reading `cw2_e4.tsv` and the shared
`kappa_e4.tsv`) are not part of this kit.

## Harvested and removed from Cannon, 2026-09-26

`run1/ run2/ kappa/ probe/ probe2/ log/ vcheck/ vcheck2/` are here, md5-identical to Cannon's (7254 files), and were
removed from `~/CL1/PartnersRivalsM1f/cannon/cw` together with the test binaries and `build/`. The Cannon kit keeps its
sources, scripts, JOBID/IDX files and the production binaries; see its `HARVESTED.md`.
