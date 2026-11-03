# Code and data: "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity"

Martin A. Nowak, Harvard University. Release 1.0.0: repository https://github.com/martin-mn/StochasticFriendlyRivals,
tag `v1.0.0`; archived on Zenodo at https://doi.org/10.5281/zenodo.23923779. MIT licence (`LICENSE`).

This release holds every program behind the paper, except those of a few second checks (see "Left out of the
release"), with its build script, compiler and flags; the Slurm scripts and task files that fix the seed of every
cluster run (Computational Supplement, Table 1); the output of the census, the sweeps and the evolutionary runs; SI
Data 1 to 3; the exact audit of a sample of the evolutionary atoms (SI Table 6); the exact certificates and the 25
exact witnesses, with the programs that check them, and the exactly certified points of the stratum counts; the exact
proofs of the memory-one statements of the SI and of the offer range of main text Eq. (4) (`PRFull/check/ref3`); the
Slurm accounting of every array (Computational Supplement, Table 2); the figure scripts; and a SHA-256 manifest of
every file. Three scripts redraw every figure from the deposited code and data, recompute every cell of SI Table 1,
write or recompute SI Tables 3 to 6 and CS Table 2, and rerun the exact certificates and proofs (and check the
deposited output of the certified points and of the audit). None of them needs a cluster. SI Table 1 is written as a
labelled plain-text table (`PRFull/Figures/sitable1.out`), not in LaTeX. No script writes SI Table 2 or SI Boxes 1 and
2; the certificate steps recompute the numbers of SI Table 2 and SI Box 2.

**This release does not contain the paper.** The paper's LaTeX sources and PDF, and those of its Computational
Supplement, are added in release 1.1.0, when the article is published; `reproduce_all.sh` then also rebuilds both. Until then the scripts skip that step and say so. References below to sections, tables and figures of the paper,
its SI and its Computational Supplement (CS) are to those documents.

The Computational Supplement (in release 1.1.0: `PRFull/Text/cs.pdf`, source `cs.tex`), posted with the paper,
describes the implementation behind every result: the limit computation and the full counts of its checks, the campaigns of the
census and the rarity sweep, the implementation and running times of the certificates, the protocols of the
evolutionary runs and of the audit, the programs (its Table 1, formerly SI Table 7) and the cost of every array
(its Table 2, formerly SI Table 8). The paper and the SI refer to it by section and table.

## Layout

The release keeps the directory layout in which the programs ran. The release root holds `PRFull/` (this paper) and
the parts of the projects it builds on that its programs read:

| directory | what is in it |
|---|---|
| `PRFull/check/` | the certificates (SI section "Certificates", SI Box 2, SI Table 2; Computational Supplement, Section 3) and their programs: `cert/` (interval engine `certeng.c`, exact engine `exact.py`, `certify.py`, `tests.py`; `results/` the certified boxes; `witnesses/` the 25 exact witnesses with their exact margins against all 65 536 binary co-players), `ubound/` (the exclusions and upper bounds), `wsep/` (the separation of the two families in W), `families/` (the segment p(t), the (a,b) family), `ref3/` (memory one: the partner conditions on the 13 boundary strata, proved twice, and SI Box 3 written from the proof; the other stability statements and the dimension of the stable set; the family (0,p_CD,p_DC,1); and the offer range of the point of main text Eq. (4); README there) |
| `PRFull/cannon/` | the Cannon kits of the revisions: `bingap` (all 65 536² binary pairs), `cw2` and `wfc2` (sensitivity runs and 4x run length, memory one and two), `sweep13` (rarity sweep at alpha = 1/3; in `validate/` the replay of its efficiency test, `rcheck.c`), `stratcert` (one exactly certified friendly rival in every retained stratum of W, N, E and the ray T>S) and `evoaudit` (the exact audit of 27 261 sampled memory-two evolutionary strategies, SI Table 6); each with its sources, `compile.sh`, `task.sh`, Slurm scripts, task files, `JOBID_*` files, README and output (the last two: see "The kits of the stratum points and of the audit") |
| `PRFull/Figures/` | the figure scripts, the figures (PDF and PNG) and their statistics (`fig45_stats.txt`, `sifig3_stats.txt`, `sitable4.out`), and SI Table 1 as `sitable1.py` recomputes it (`sitable1.out`) |
| `PRFull/Text/` | SI Data 1-3 (`SIData1_games512.csv`, `SIData2_sensitivity.csv`, `SIData3_runlength.csv`) and SI Table 4 as written by `Figures/sitable4.py`, without its legend (`si_table4_body.tex`); the paper (`ms.tex`, `refs.tex`, `reorder_bib.py`, `ms.pdf`, with `si_table4.tex`, and SI Box 3 with its title and legend, `si_box3_strata.tex` and `si_box3_legend.tex`) and the Computational Supplement (`cs.tex`, `cs_bib.py`, `cs.pdf`) follow in release 1.1.0 |
| `PRFull/release/` | this file, the scripts below, `accounting/` (the Slurm accounting), `census_class/` (the support-based efficiency class of every stratum), `figure_inputs.json` (every file each figure script reads, traced), `tools/` |
| `PartnersRivalsM1f/` | the memory-one programs (`cannon/m1f_pairs.c`, the evolutionary program `cannon/cw/cw.f` with `m1cls.c`, `kappa.c`), their Slurm scripts, per-replicate output (`cannon/cw/run1`, `run2`, `packed`, `kappa`) and the checks of SI section "Memory one" (`check/`) |
| `PartnersRivalsM2f/` | the memory-two programs (`cannon/m2f.c`, `census.c`, `dimprobe.c`, `wfc/wfc.c`, ...), the census output (`check/st_eff.txt`, and in `cannon/ca1_out/` what was only on Cannon: the confirmed points at full double precision `pt_*.txt`, the raw output of the dimension probe `dim_*` and of the predicate-free sweep `sw_*`, the census probe `probe_45749304_*`; `cannon/ca1_out.md5` is the md5 of each of these 2440 files taken on Cannon, and every copied file matches it), the evolutionary output (`check/wfc_a333`, `check/wfc_f2`, and the raw output with the sampled strategies in `cannon/archive/*.tgz`), the alpha = 0.333 sweep (in `cannon/archive/cannon_2026-09-26.tgz`) and the exact checks (`check/`); and, from the scratch of the derivation pass, which the release otherwise leaves out, the programs and output of five checks (`check/working/`, README there): `t2_fw.c` (the efficient sets of all 3^16 strata by minimum arborescences), `effcc_task1/effcc.c` with its logs and `t2_altform2.c` (the closed forms of efficiency below and above the switch line, each checked as a set equality on all 3^16 strata), `verify_task1/` (every stratum of W of dimension 8 to 15: `dim8W.py`, `dim8b.py`, `dimNW.py` and their logs) and `task2/kit/` (the kit of the re-censuses of W and N: `t2.c`, `t2b.c`, their build and Slurm scripts and `*_JOBID` files), with the collected lines of its two arrays, `task2/rec_W.txt` (364 strata) and `task2/rec_cannon.txt` (5022) |
| `PartnersRivals1/`, `PartnersRivals/`, `DiskM1WF/`, `DiskM2WF/`, `Donation/` | only the files the figure scripts and the kit analyses read (`figure_inputs.json`): the drawing code of the binary paper (PR2026) and the data of its binary runs, which Figures 4 and 5 and SI Table 4 compare with (0.73 GB). The binary paper's own release is PRFR v1.1.0 (https://github.com/martin-mn/PRFR, https://doi.org/10.5281/zenodo.23079603). None of the 17 files from `PartnersRivals1/` is byte-identical to a file of that release (checked on 2026-10-22 against the tag v1.1.0, commit e6ca720): 12 are not in it (the run dumps `Data/f1_e4.*` and `f2_e4.*`, `Cannon/atoms512.bin`, and the helper scripts `Figures/disklib.py`, `FinalFigures/atomcounts.py`, `atomshares.py`, `wfrep.py`); `Figures/m2_masks.npz` holds the same arrays as its `data/census/m2_masks.npz`, with a different provenance note; and `FinalFigures/fig4.py`, `fig5.py` and `wfdata.py` are the working versions from which its published scripts were made. They are deposited here as the figure scripts of this paper read them |

`CONTENTS.tsv` lists every file with its size and the promise of the statement it belongs to; `MANIFEST.sha256` gives
its SHA-256.

## Requirements

Tested on macOS 26.6 (Apple M2 Pro, 12 cores, 16 GB), with Anaconda Python 3.10.9 (x86_64, under Rosetta).

| | version used | needed for |
|---|---|---|
| Python | 3.10.9 (the only version tested) | everything |
| numpy | 1.23.5 | everything |
| scipy | 1.10.0 | figures (fig45, sifig3), the kit analyses, the defensibility LP |
| matplotlib | 3.7.0 | figures |
| gmpy2 | 2.1.2 | the exact engine of the certificates (`check/cert/exact.py`) |
| sympy | 1.11.1 | the exact checks of memory one (`PRFull/check/ref3`, `PartnersRivalsM1f/check/symbolic_rivals.py`); `reproduce_certificates.sh` stops at its start without it |
| C compiler | Apple clang 15.0 locally; GCC 8.5.0 on Cannon | the certificate engines (`-O2 -ffp-contract=off -fno-fast-math`: the outward rounding needs it), the census, the memory-two evolution, the replay of the efficiency test (`cannon/sweep13/validate/rcheck.c`) |
| Fortran compiler | Intel oneAPI ifx 2025.2.1 on Cannon (`-O3 -xHost -fp-model precise`); GNU Fortran 14.1 locally (`-O2 -ffp-contract=off -fno-fast-math`) | only to rebuild the memory-one evolutionary program `cw.f` (not needed for the verification) |
| pdflatex | TeX Live 2024 (pdfTeX 1.40.26) | rebuilding the paper and the Computational Supplement (release 1.1.0) |
| bash, rsync, shasum or sha256sum | bash 3.2 or later | the scripts |

Hardware: the quick verification needs 3 cores, 3 GB of memory and 0.9 GB of disk for its working copies (0.2 GB,
and while the step `evoaudit` runs the 0.67 GB of strategy dumps it extracts from the archives; the release itself is
3.1 GB). The full rerun of the certificates writes about 1.5 GB of intermediate files (`check/ubound`); its
peak memory was not measured as one run (each of its steps was run before, on the same machine, within 16 GB). No
compiled code is deposited; the scripts build what they need.

## Minimal verification (no cluster; about 13 minutes)

From the release root:

    shasum -a 256 -c MANIFEST.sha256 | grep -v ': OK$'       # nothing printed: every file is the deposited one
    python3 PRFull/release/relocate.py                       # point the programs' absolute paths at this tree
    bash PRFull/release/reproduce_all.sh --quick             # figures, tables, a subset of the certificates (and, in 1.1.0, the paper and the CS)

`relocate.py` is needed because many programs name their inputs by absolute paths on the author's machine
(`/Users/martin/Documents/CL1/...`); it rewrites that prefix in the program files only and keeps a record in
`PRFull/release/relocate.json` while the tree is relocated. `relocate.py --undo` restores the files and removes the
record, after which the manifest checks again and the tree holds no file that the manifest does not list (apart from
the scripts' working directory, if it was left in the tree; see below).

### From the GitHub repository alone

The repository carries everything except the large data, which the Zenodo record holds in archives
(`PRFull/release/hosting/hosting_plan.txt` lists what is where; `hosting/zenodo_archive_of_file.tsv` gives the archive
of every file). In a clone of the repository, `shasum -a 256 -c --ignore-missing MANIFEST.sha256` checks the files it
has, and `reproduce_all.sh --quick` runs every step but three, which need data the repository does not carry: `fig45`
and `sitable4` (Figures 4 and 5, SI Figures 1 and 2, SI Table 4: the binary runs' strategy dumps) and `evoaudit` (the
audit of the atoms: the production archives of the memory-two runs). Those steps are skipped, each with a line
`SKIPPED: missing ... ; they are in the archive PRFull-1.0.0_<name>.tar.gz of the Zenodo record`, and the summaries
count them. (The raw output of the dimension probe, `PartnersRivalsM2f/cannon/ca1_out/dim_*` (15 MB), which the step
`sitable1` reads, is in the repository; the rest of the census output that was only on Cannon is in the archive
`census`.) To run them,

    bash PRFull/release/fetch_zenodo.sh          # fetches the two archives (about 0.7 GB compressed), unpacks them
                                                 # into the clone and checks every file against MANIFEST.sha256

The full scripts skip two more steps in the same way, whose archives the default fetch does not bring: `table5` of
`reproduce_figures.sh` (SI Table 5: the raw output of the sensitivity runs, `PRFull/cannon/cw2/out` and `wfc2/out`,
in the archive `evolution-runs`) and `separation` of `reproduce_certificates.sh` (the deposited runs of the further
classes of the separation, `PRFull/check/wsep/out/extra`, which its comparison reads, in the archive `certificates`).
`fetch_zenodo.sh evolution-runs certificates` fetches these two, and `--all` all five data archives, everything the
full `reproduce_all.sh` reads.

Apart from the relocation, nothing in the tree is changed by the scripts. Each one first runs `relocate.py` on the
tree (a no-op in the author's tree and in a tree already relocated; `relocate.py --undo` reverses it), then works on
copies in `PRFull/release/work` (or `--out DIR`; a DIR outside the tree adds nothing to it), runs every program
under `nice -n 10`, and compares every file a program writes with the deposited file at the same place. The C code is
built in the copies too: the library `libm2f`, which the self-tests of the certificate library compare with, is built
from copies of its sources in `DIR/certificates/PartnersRivalsM2f`, not in the tree. "The same" allows only what
cannot be reproduced: the creation date inside a PDF, run times, dates and paths in text, lines that several threads
write in another order, and, where a program picks the first of several exactly tied co-players, which one it names
(this is said in the step's verdict). The summary of a run is `work/logs/summary.txt`; the exit status is 0 when
every step passed and every rewritten file is the deposited one. `reproduce_certificates.sh --only A,B,...` reruns
single steps (the names `--list` prints) in an existing working copy, and `RCHECK_TASKS='0 1 2'` chooses the tasks of
the rarity sweep that the step `rcheck` replays. `PYTHON=...` chooses the interpreter of both scripts and `CC=...` the
C compiler of `reproduce_certificates.sh`, also in `check/ubound/run_all.sh` and `check/wsep/run_all.sh`, which its
steps `exclusions` and `separation` run.

Measured on the machine above; the quick run below is the one of 2026-11-03, 09:04-09:18, on a tree assembled from
these files (`reproduce_figures.sh --quick`, then `reproduce_certificates.sh --quick`: the two scripts that
`reproduce_all.sh --quick` runs), and the times of the full runs are sums of the recorded times of their steps:

| script | what it does | wall-clock | peak memory |
|---|---|---|---|
| `reproduce_figures.sh --quick` | Figures 1-5, SI Figures 1-3, SI Tables 1, 3 and 4, CS Table 2, SI Data 1-3, from the deposited data | 2 min (127 s) | 1.4 GB |
| `reproduce_certificates.sh --quick` | the subset of the certificates and proofs listed below, and the checks of the deposited certified points and audit | 11 min (664 s) on 3 cores | 2.5 GB |
| `reproduce_all.sh --quick` | both, and in release 1.1.0 the LaTeX builds (the paper and the Computational Supplement, 5 s on 2026-09-28) | 13 min (791 s) | 2.5 GB |
| `reproduce_figures.sh` | the quick set, the receipts of Figures 2 and 3 (3 min on 3 processes) and the SI Table 5 analyses (31 s) | about 6 min (the sum of its steps, each run separately) | 1.3 GB |
| `reproduce_certificates.sh` | every certificate and proof, and the same two checks | about 2 h 20 min on 3 cores (the sum of the recorded timings of its steps, the replay of all 512 tasks of the rarity sweep about 17 min of it, at 5.8 s of CPU per task for the x86_64 binary the step builds under Rosetta, measured on 2026-11-03; not run as one job here) | not measured |
| `reproduce_all.sh` | everything | about 2 h 30 min | not measured |

So `reproduce_figures.sh` alone takes about 6 minutes, and the full `reproduce_certificates.sh` about 2 hours and 20
minutes, by the
recorded times of their steps; neither has been run as one job, and the full rerun's peak memory has not been
measured.

The certificate steps of the quick run, in seconds: build 8, library self-tests 119, W box 86, 5-cube 40, p(t) 51,
six witnesses 60, exclusions 130 (2.5 GB: the classes of all 3^16 strata), separation 80, binary pairs 20, the
partner conditions of memory one 5 (`m1_boundary`), the two proofs compared 7 (`m1_compare`), SI Box 3 6 (`m1_box`),
the family (0,p_CD,p_DC,1) 5 (`m1_cycle`), the offer range 3 (`offer`), the replay of two tasks of the rarity sweep 7
(`rcheck`), the certified points of the stratum counts 5 (`stratcert`), the audit of the atoms 10 (`evoaudit`: about 3
to extract the dumps from the archives, the rest for `analyze.py`, which redraws the design from all 4.8 M records;
1.2 GB). Of the figure steps, `sitable1` (SI Table 1) took 8 s and 0.6 GB. (The peak memory of a step is that of the
largest process it started.)

What the quick run checks, and what it must print (the scripts assert these lines, or compare the output that holds
them with the deposited file):

* **Figures and tables** (`reproduce_figures.sh --quick`): every figure PDF identical to the deposited one apart
  from its creation date, every PNG pixel for pixel, `fig45_stats.txt`, `sifig3_stats.txt`, `sitable4.out`,
  `si_table4_body.tex`, `sitable1.out` (SI Table 1), SI Data 1-3 and `sweep13.json` (SI Table 3) identical (SI Data 2
  and 3 byte for byte, rebuilt by `tools/mk_sidata23.py`, written for the release: the tree had no program for them);
  SI Table 1 recomputed by `Figures/sitable1.py` from the census (`PartnersRivalsM2f/check/st_eff.txt`), the raw output
  of the dimension probe and the support-based classes, with all its checks passed (`checks failed: 0 of 29`): each of
  the 384 stripes of the probe agrees with its task log, the class that the census and the probe read from the weights
  is the support-based class at every stratum they list, the strata of dimension 0 are, as sets, the binary friendly
  rivals of W, S, N and E (8, 1519, 80 and 80) and the binary strategies efficient on the line T = S (7639 below the
  switch line, 3072 above it, 14 757 at (1,0)), and the C/D mirror maps N onto E; CS Table 2 (the cost, in the
  Computational Supplement; in earlier drafts it was SI Table 8 and SI Table 7, and the files
  `accounting/check_table7.py`, `table7_jobs.tsv`, `sacct_table7.psv` and the step `table7` keep that number)
  recomputed from `accounting/sacct_table7.psv`: every row and the total of 8569 core-hours over 53 arrays agree
  (`total  53 arrays  core-h  8568.81  CS Table 2  8569  ok`,
  `rows or arrays that differ: 0`).
* **Certificates** (`reproduce_certificates.sh --quick`):
  - the certificate library's self-tests (b) engine against exact rationals, (c) the corner lemma L3, (e) the
    certificate primitives: `==== 23 checks, 0 failed, 1 skipped`. The skipped check (in (e), the exclusion of 300
    random efficient strata with p8 free) read a sample file that was kept only in a working directory of the author
    and no longer exists; it says so (`SKIPPED e.exclude_stratum ...`), and the exclusion certificates themselves are
    rerun and checked by the step `exclusions`;
  - the certified box around the W witness of SI Table 2 (7-dimensional, all 2^7 x 65 536 vertex margins settled) and
    the whole open 5-cube of wedge W: `W_box.json`, `cube5_stratum.json` identical to the deposited files;
  - the segment p(t) of the main text, exactly: `CERTIFIED: p(t) is a friendly rival of wedge N for every t in (0,1)`;
  - six of the 25 witnesses regenerated in exact arithmetic: each exact-margins file has the SHA-256 recorded in the
    deposited witness (`6 regenerated witnesses, 0 differ from the deposit`);
  - the exclusions: the exact efficiency class of all 3^16 strata (`classes: none 31580692 C 6206077 A 2788386 Z
    2471566`), the 12-co-player certificate of the upper bounds, its 786 432 rows recomputed in exact rational
    arithmetic, and the census consistency (`class mismatches with the exact classes here: 0`);
  - the separation in W: the negative control (all 364 friendly-rival strata of W survive: `NOT_EXCLUDED 364`), the
    class c9a (18 108 strata, each excluded; `rows_failing 0`), the lemma check;
  - the census of binary pairs on the slice of Cannon's probe (7.9 M pairs): `min_nonzero 1/902`, 0 flagged, 0
    symmetry violations, identical to the probe's output;
  - the certified points of the stratum counts (`cannon/stratcert`, as deposited: the cluster run is not repeated):
    `collect.py out 64` finds every one of the 12 748 units once (`units 12748, lines 12748, duplicated 0, malformed
    0, missing 0 []`, `tasks without .done (of 64): 0 []`), audits every line and recomputes the exact self-play limit
    of every witness (`audit of 12748 certified lines (exact self-play recomputed for 12748): no problem`), and prints
    `verdict: every retained stratum of W, N, E and ray+ has a certified exact friendly rival: True`; its output is
    identical to `cannon/stratcert/collect.out`, and `out/stratcert_all.tsv` to the deposited one;
  - the exact audit of the evolutionary atoms (`cannon/evoaudit`, as deposited): the strategy dumps are extracted
    from the two archives (`dumps: 1024 files, 4772732 records (F1 1453629, F2 3319103)`, `sampled strategies equal
    their records in the dumps: 27261 of 27261`), `src/analyze.py` redraws the design, asserts that it gives the
    deposited sample and analyses the 128 chunks: `audit: 27261 of 27261 sampled strategies audited, 0 errors, 0
    missing, 2 changed atoms, 0 not reproduced` and `audited 27261, errors 0, missing 0, unresolved 0; changed atoms
    2 (010->000 2); changed bits: E 0, C 0, S 2; validated 101: 0`; `results/` (`summary.json`, `confusion.txt`,
    `records.tsv.gz`, `analyze.out`) identical to the deposited files;
  - the partner conditions of memory one on the 13 boundary strata (`check/ref3/eq3_boundary/P1/prove_strata.py`; SI
    Box 3): `OVERALL PROVED: True`, `prove_strata.json` and `prove_strata.log` identical to the deposited files; the
    two proofs compared (`referee/referee_check.py`): `problems: 0`, its output identical to `referee_check.log`;
  - SI Box 3 written from the proof and read back (`box/mk_strata_table.py`, `check_box_vs_exact_partners.py`):
    `wrote si_box3_strata_body.tex: 13 strata, 34 gains that differ from the display, 65 Bernstein factors` and
    `typeset box agrees with exact_partners.limitW: 204 listed gains, 1044 display gains, 78 L>0 checks; 0
    mismatches`; in this release the box is written without its title and legend (`si_box3_strata_body.tex`), which
    are part of the paper;
  - the family (0,p_CD,p_DC,1) on its nine strata (`check/ref3/cycle_family/prove.py`): `positive-dimensional strata
    all match the claim: True`, its output identical to `prove.log`;
  - the offer range of the point of main text Eq. (4), from the deposited exact margins of all 65 536 binary
    co-players (`check/ref3/eq4_example/offer/interval.py`, `check_eight.py`): no numerator negative as x -> 0
    (`numerator negative as x -> 0+: 0 []`), the smallest root 5/13, shared by 400 co-players, so the point is a
    friendly rival exactly for 0 < x <= 5/13, and the margins of the eight binding co-players `001*00011011*0*1`
    equal (5-13x)/(2(25-13x));
  - two tasks (1 and 3) of the rarity sweep replayed (`cannon/sweep13/validate/rcheck.sh`): `replayed 2 tasks: RC
    lines differ 0, FLAG lines differ 0, ALPHA counts differ 0`.
* **The paper and the Computational Supplement** (release 1.1.0; in this release the step is skipped, with the
  message `LaTeX: skipped, the paper's sources are not in this release (PRFull/Text/ms.tex, cs.tex; they are added in release 1.1.0, when the article is published)`): `reorder_bib.py` leaves `ms.tex` unchanged and `cs_bib.py` leaves
  `cs.tex` unchanged, and pdflatex runs on both without errors or undefined references (on 2026-10-26, from the
  paper's files of release 1.1.0 assembled alone: `ms.pdf` 70 pages and `cs.pdf` 14 pages, no LaTeX warnings and no
  overfull boxes in the last pass of either).

In the run of 2026-11-03 every step passed and no rewritten file differs from the deposited one; so did the run of
2026-10-26 (786 s), before the step `sitable1` was added. Earlier quick runs, before the steps of `check/ref3` and the
replay were added, passed in the same way: 766 s on 2026-09-27, and 763 s on 2026-09-28 (figures 115 s, certificates
642 s, the two LaTeX builds 5 s; largest peak memory 2.5 GB, the exclusions).

The full `reproduce_certificates.sh` adds the other certified boxes (N, E, S, the 15-dimensional S box, the large W
box), the (a,b) family and the boundary of the S witness, all 25 witnesses, the N witness in exact Fractions, the
exclusions in every dimension with their independent re-check (14.2 M exclusions), the separation of all classes with
the greedy cover, and all 64 tasks of the census of binary pairs (smallest non-zero margin 1/4140). It also reruns
the rest of `check/ref3`: the first proof's limits against `exact_partners.py`, the second proof (245 exact limits and
the root search), the open facets and the zero-determinant plane, the other memory-one stability statements
(`certify.py`, which leaves one case open and prints FAIL for it, and `d1_switchline.py`, which proves that case), the
exact sweep over all 81 strata, the 405 769 exact tests of the family (0,p_CD,p_DC,1), the exact margins of all 65 536
binary co-players at four values of the offer probability, and the replay of all 512 tasks of the rarity sweep (about
17 minutes on 3 cores). The steps
`stratcert` and `evoaudit` are the same in both modes: they check the deposited cluster output, which neither mode
reruns.

## What needs the cluster

The cluster computations (the census, the sweeps, the evolutionary runs, the census of binary pairs, the certified
points of the stratum counts, the audit of the atoms) took 8569 core-hours on Harvard's Cannon cluster (Computational Supplement, Table 2)
and, apart from the census of binary pairs, which the full `reproduce_certificates.sh` repeats, are not rerun by these
scripts; their output is deposited, the scripts redraw every figure and recompute SI Table 1 and SI Tables 3 to 6
from it, and they check the certified points and the audit as deposited.
Each kit (`PRFull/cannon/*`, `PartnersRivalsM1f/cannon`, `PartnersRivalsM1f/cannon/cw`, `PartnersRivalsM2f/cannon`,
`PartnersRivalsM2f/cannon/wfc`, `.../alpha`) has a README with its build, its probe, its Slurm scripts and task files
(seeds set by the task index or listed in the task files) and its collection scripts; the kit of the re-censuses of W
and N is `PartnersRivalsM2f/check/working/task2/kit`, described in `check/working/README.md`. Every task of every
production array COMPLETED (`sacct`, per task); tasks were cancelled, timed out or failed only in six probe and
re-census arrays (45630088, 45658311, 48329275; 45785196, 45785348, 45786256), which `accounting/check_table7.py`
lists and CS Table 2 counts. `accounting/sacct_table7.psv` is the `sacct` record of every task and job step of
the 53 arrays (`accounting/dump_sacct.sh` wrote it, last on 2026-09-28 02:46 UTC, when it added the four arrays of
`stratcert` and `evoaudit`, 48975432 and 48981870 of `stratcert`, 48975385 and 48981982 of `evoaudit`, all 192 of their
tasks COMPLETED; the records of the other 49 arrays are byte-identical to the dump before;
`accounting/check_table7.py` recomputes CS Table 2 from it), and `accounting/jobs_tsv_lines.tsv` their lines in the
submission ledger (52: the probe 45656466 has none).

## The census and the efficiency class of every stratum

`PartnersRivalsM2f/check/st_eff.txt` is the census output as first archived: one line `g cls nP nM K` for each of the
1 524 791 efficient strata in which a draw was confirmed competitive, with the class read from the weights at the
stratum's first draw with the tolerance 1e-12. `PRFull/release/census_class/` adds:

* `census_st_45751662.tgz`: the stratum output of all 512 census tasks as the cluster wrote it (3 254 063 lines,
  also the 1 729 272 strata that the census read as not efficient), with the task logs; `census_st_45751662.md5`
  holds the md5 of each file taken on Cannon, and every file matches it;
* `census_support_class.u8.gz`: the efficiency class decided by the support of the self-play limit alone (lemma L1,
  exact integer exponents, no tolerance) for all 43 046 721 strata, one byte per stratum in the census's own code
  (`ub.c`, mode `classes`; `mk_census_class.py`);
* `census_all_support_class.txt.gz`: every census line with its support-based class as a sixth column
  (`join_census.py`), and `st_eff_support_class.txt.gz`, the same for `st_eff.txt`;
* `census_class.json`: the counts and checks. The support-based classes of all strata are those of the Methods (C
  6 206 077, A 2 788 386, Z 2 471 566), and on every one of the 3 254 063 census lines the class the census read
  from the weights equals the support-based class (0 mismatches).

The census output that was only on Cannon (`~/CL1/PartnersRivalsM2f/ca1/out`) was copied read-only on 2026-09-27
(`assemble_release.py --harvest-cannon`: 2440 files, 842 MB; `pt_*.txt`, the confirmed points at full double
precision, up to 4 per stratum and sign, 827 MB; the raw output of the dimension probe `dim_*`, 15 MB, and of the
predicate-free sweep `sw_*`, 0.3 MB, whose summaries are `check/dim.out` and `check/sweep.out`; the probe
`probe_45749304_*`). Every file matches the md5 taken on Cannon (`PartnersRivalsM2f/cannon/ca1_out.md5`), and the
release puts them in `PartnersRivalsM2f/cannon/ca1_out/`.

`check/census.py` and its output `census.out` give the counts of the census alone (W 362, N 5020). SI Table 1 takes
in the dimension probe as well: `PRFull/Figures/sitable1.py` (step `sitable1`) recomputes every cell of its panels a
and b from `st_eff.txt`, the raw output of the probe (`ca1_out/dim_*`) and `census_support_class.u8.gz`, and writes
them, with the strata the probe adds to the census in each locus (two in W and two in N, which give 364 and 5022) and
the dimension histograms, to `PRFull/Figures/sitable1.out`; its counts and histograms of the six blocks of the probe
are those of `check/dim.out`.

## The kits of the stratum points and of the audit

Both ran on Cannon on 2026-09-27; every task of their four arrays COMPLETED (`sacct`, per task), every Slurm stderr
file is empty, and their output is deposited as fetched. Each README gives the build, the probe, the sizing, the
submission and the checks.

`PRFull/cannon/stratcert`: one exactly certified friendly rival in every retained stratum of W, N, E and the ray
T>S (SI section "Certificates", "The stratum counts", and SI Box 2; Computational Supplement, Section 3 and Table 2:
probe 48975432, 4 tasks; full run 48981870, 60 tasks; 30.8 core-hours).

* `units.tsv`: the 12 748 units (stratum, sign): 7726 of sign + (the ray T>S, W's 364 and N's 5022 among them) and
  5022 of sign - (E); `prep/` how they were chosen (`mkstrata.py` with `ub.c`, and the pool of co-players
  `pool_q.txt`; log `mkstrata.out`) and the first checks (`validate.py`); `points.tsv`: the 50 897 census points they
  start from, taken from `pt_*.txt` (md5 in `prep/pt_md5_cannon.txt`);
* `stratcert.py` (the ladder of candidates for each unit), `cert/` (a verbatim copy of `PRFull/check/cert`, checked by
  `SHA256SUMS`), `compile.sh`, `task.sh`, `probe.sbatch`, `full.sbatch`, `JOBID_probe`, `JOBID_full`, `logs/`;
* `out/t_0.tsv` ... `t_63.tsv`: one line per unit, with the witness in exact fractions, its exact self-play limit and
  class, the co-players by decision route and the exact minimum margin; `t_N.done`, `time_N.txt`; `stratcert_all.tsv`,
  all 12 748 lines sorted;
* `collect.py`, the checker (completeness, an audit of every line, the exact self-play of every witness recomputed,
  per-locus tables, the verdict), and `collect.out`, its output on the deposited outputs (verdict True).

`PRFull/cannon/evoaudit`: the exact audit of 27 261 sampled memory-two evolutionary strategies (SI Table 6; Computational
Supplement, Section 4 and Table 2: probe 48975385, 4 tasks; full run 48981982, 124 tasks; 10.7 core-hours).

* `src/`: the loaders, the design of the sample (`sample.py`), the audit (`audit.py`), the analysis (`analyze.py`),
  the projection from the probe, the local tests, production's classifier rebuilt from its source (`prodcls.c` with
  production's `m2f.c`), `probe_version/` (the audit and the analysis as they ran in the probe) and `review/` (the
  scripts of the adversarial review); `lib/cert/`, a copy of the certificate library; `games.dat`; `compile.sh`,
  `task.sh`, `probe.sbatch`, `run_full.sbatch`, `stage.sh`, `collect.sh`, `check.sh`, `JOBID_probe`, `JOBID_full`,
  `logs/`;
* `sample/sample.npy` (the 27 261 strategies with their inclusion probabilities) and `sample/design.json`;
* `out_probe/` (chunks 0-3) and `out/` (chunks 4-127): one JSON line per strategy, and a `.done` line per chunk;
* `results/`: `summary.json`, `confusion.txt`, `records.tsv.gz` (one line per strategy) and `analyze.out` (the output
  of `src/analyze.py out_probe,out results`): 27 261 audited, 0 errors, 0 unresolved, 2 changed atoms (both `010` ->
  `000`, near-ties of stability below production's tolerance), no bit of efficiency or competitiveness changed;
  `results_probe/`, the probe's partial analysis.
* **The strategy dumps it samples are not duplicated.** The kit's `data/F1`, `data/F2` on the author's machine
  (643 MB, not in the release) is a local extraction of the production dumps `S_1o3_<task>_1.bin` (and `r_1o3_*`)
  from the deposited archives `PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz` (`PartnersRivalsM2f_wfc/out/`,
  N = 1000) and `cannon_out_f2_2026-09-26.tgz` (`out_f2/`, N = 100); the 1024 `S_` files there are byte-identical to
  it (checked 2026-09-27). `analyze.py` needs them to redraw the design; `PRFull/release/tools/evoaudit_dumps.py
  extract PartnersRivalsM2f/cannon/archive PRFull/cannon/evoaudit` writes them from the archives (3 s) and checks every
  sampled strategy against its record, which is what the step `evoaudit` of `reproduce_certificates.sh` does, on a
  copy.

Run by hand, `collect.py out 64` and `src/analyze.py out_probe,out results` rewrite `out/stratcert_all.tsv` and
`results/` in place (with the same content); the scripts run them on copies.

## Left out of the release

* The paper and its Computational Supplement, their sources and PDFs (`PRFull/Text/ms.tex`, `refs.tex`,
  `reorder_bib.py`, `ms.pdf`, `cs.tex`, `cs_bib.py`, `cs.pdf`): in release 1.1.0, when the article is published
  (`assemble_release.py --with-paper`).
* Every other manuscript and the author's notes: the companion paper's manuscript (its programs and data that the
  figure scripts read are included), the drafts in `PartnersRivalsM1f/Text` and `PartnersRivalsM2f/Text`, the
  derivation notes `PartnersRivalsM2f/check/derivations`, the project notes and `PRFull/README.md` (this file is the
  release's README), and a third-party article kept as a reference text. Six project READMEs of `PartnersRivalsM1f`
  and `PartnersRivalsM2f` are shipped as release-only copies (`PRFull/release/overrides`): the same technical
  documentation, without the drafts' text and the authors' working notes. SI Table 4 is shipped without its legend
  (`si_table4_body.tex`); `si_table4.tex` and the legend file `Figures/sitable4_legend.tex` are part of the paper. SI
  Box 3 is shipped without its title and legend (`check/ref3/eq3_boundary/box/si_box3_strata_body.tex`);
  `box/si_box3_strata.tex`, `Text/si_box3_strata.tex` and the legend file `Text/si_box3_legend.tex` are part of the
  paper (`mk_strata_table.py` writes the whole box when the legend file is present). Figure 1 takes the values it once read from
  the drafts from `PRFull/Figures/fig1_reference_values.json`. `tools/manuscript_guard.py`, which
  `assemble_release.py` runs on every selected file, checks that no manuscript text is in the release, and
  `assemble_release.py --dest` refuses to assemble while it flags a file.
* The scratch of the memory-two derivation pass, `PartnersRivalsM2f/check/working` (196 MB of one-off drivers, dumps
  and masks), except the programs and output of the five checks that its README lists (`t2_fw.c`, `effcc_task1/effcc.c`,
  `t2_altform2.c`, `verify_task1`, `task2`; see "Layout").
* The programs of some second checks that the SI or the Computational Supplement describes (the Computational
  Supplement says at each that its program is not deposited). Not on the author's disk: the exact determinant program
  that confirmed the rows of ALLD and ALLC and found negative margins at 4611 interior points of 1537 excluded strata;
  the program that decided again every efficiency and ALLD decision of the separation in W, with orders found by
  minimum arborescences; and the re-census of the 60 165 candidate strata of W at 512 points each, with the trial of
  the looser screening rule on 385 088 draws in a tenth of them. On the author's disk but not deposited, for want of
  their logs: the numerical check of the reduction to binary co-players against 1.35 x 10^8 stochastic co-players
  (`claimC.c`, `claimC2.c` and `t3_claimC_cont.c` in `PartnersRivalsM2f/check/working`). Kept only in the scratch of
  the review session of 2026-09-27: the scripts that recomputed every margin of eight certified points of the stratum
  counts in exact rational arithmetic (`PRFull/cannon/stratcert/README.md`, "Adversarial review").
* Compiled code (the sources and build scripts are deposited), `__pycache__`, virtual environments and user-level
  package installs (`venv/`, `pylib/`), LaTeX build files, the working copies of the scripts.
* `PRFull/cannon/evoaudit/data` (see above: the same files are in the archives).
* The rebuilt intermediates of `PRFull/cannon/stratcert/prep` (`cls.bin`, `c2.idx`, `c2.qm`, `c2.qp`, 173 MB): `cd
  prep && ./ub classes cls.bin 3 && ./ub cover cls.bin ../pool_q.txt c2 3` writes them in about 2.5 minutes
  (`c2.log`, the log of the second, is deposited).

## Assembling the release

    python3 PRFull/release/assemble_release.py --dry-run [--files]    # what would be copied, with sizes; what is missing
    python3 PRFull/release/assemble_release.py --dest DIR              # the release tree, MANIFEST.sha256, CONTENTS.tsv
    python3 PRFull/release/assemble_release.py --with-paper --dest DIR # release 1.1.0: the same with the paper and the CS

On 2026-11-03 the dry run selects 27 911 files, 3.08 GB (with `--with-paper` 27 923; `README.md`, `LICENSE`,
`reproduce_all.sh`, `MANIFEST.sha256` and `CONTENTS.tsv` at the top level come on top), and flags nothing, neither
manuscript text nor a promised file that is missing. `--dest` refuses to assemble while the manuscript guard flags a
file. In an assembled release, which does not hold the author's `tools/manuscript_sources.json`, the guard can check
only names (`[guard-off]`), and `--dest` proceeds. This tree was assembled with `--dest` on 2026-11-03.

## Hosting and packing

`hosting/hosting_plan.py TREE` splits an assembled release between the GitHub repository (code, scripts, documentation,
figures, SI Data and the smaller results: about 7000 files, 191 MB, no file over 50 MB) and the Zenodo record (the
whole release: README.md, LICENSE, MANIFEST.sha256, CONTENTS.tsv and six archives: code, binary-paper-dumps,
m2-archives-and-audit, census, evolution-runs, certificates; 3.08 GB). The deposited `hosting/` files were computed
from the assembled tree before they were themselves added to it. `pack_release.sh TREE OUTDIR` writes the archives:
files in sorted order under the prefix `PRFull-1.0.0/`, owner and group root, no extended attributes, no macOS
metadata, no ACLs, `gzip -n`; and `ARCHIVES.sha256`.

## Licence

MIT licence, copyright (c) 2026 Martin A. Nowak; see `LICENSE`. It covers the code and the data of this release,
the files from the companion projects included (they are by the same author).
