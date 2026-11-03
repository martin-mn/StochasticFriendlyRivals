# `check/working/`: the scratch of the derivation pass, and the part of it that is deposited

`working/` holds the scratch of the memory-two derivation pass: about 1100 files (540 at its top level) and 196 MB of
one-off drivers, sweep dumps and masks (`../README.md`, section "`working/` and `derivations/`"). The code and data release leaves it out,
except the programs and outputs of five checks that the paper or its Computational Supplement (CS) describes and that
ran here. They are listed below and deposited as they ran; none of them was changed for the release. This file was
added on 2026-11-03, for the release.

| deposited | check |
|---|---|
| `t2_fw.c` | the efficient sets of all 3^16 strata by minimum arborescences (CS Table 1) |
| `effcc_task1/effcc.c` and its logs `effcc_task1/ver_0.log` ... `ver_11.log` | the closed form of efficiency below the switch line, compared with the limit on every one of the 3^16 strata (SI section on memory two; CS Table 1) |
| `t2_altform2.c` | the characterisation of efficiency above the switch line, compared in the same way (the same) |
| `verify_task1/dim8W.py`, `dim8b.py`, `dimNW.py` and their logs `dim8.log`, `dim8b.log`, `dimN.log` | the second of the two further campaigns in W (CS Section 2): every stratum of dimension 8 to 15 |
| `task2/kit/` (sources, build and Slurm scripts, `*_JOBID`) and `task2/rec_W.txt`, `task2/rec_cannon.txt` | the re-censuses of W and N of the Methods (`t2.c`, `t2b.c`; CS Table 1) |

The programs of the first further campaign in W (the re-census of the 60 165 candidate strata at 512 points each) are
not on the author's disk and are not deposited.

## `t2_fw.c`: the efficient sets by minimum arborescences

For each of the 3^16 = 43 046 721 strata, `t2_fw.c` builds the self-play graph of the pattern at ε = 0, its closed
classes, and for every pair of classes the least number of errors that leads from the one to the other. The support of
the ε → 0 limit is then the union of the classes whose cheapest spanning in-tree, an integer cost found by Edmonds'
minimum arborescence, is least. The program compares this support, state by state, with the support of the limit that
`freqs16` of `../../cannon/m2f.c` computes at one interior draw of the stratum (weights above 1e-13), and counts the
strata accepted by three targets, from the arborescences (`p`) and from the limit computation (`t`), and the strata on
which the two disagree (`m`): ALT, support inside the states 4 to 11 (efficient above the switch line); LINE, support
outside the states 12 to 15 (efficient on it); CC, support inside the states 0 to 3 (efficient below it).

    cc -O3 -o t2_fw t2_fw.c ../../cannon/m2f.c -lm
    for s in 0 1 2; do ./t2_fw $s 3 2> fw_$s.err & done; wait      # stripe s of 3; 50 s on 3 cores

The output of the original run (12 stripes) was not kept. The run above, made on 2026-11-03, wrote

    stripe 0: suppMis 0 maxNclasses 9 | ALT t 1122706 p 1122706 m 0 | LINE t 1449522 p 1449522 m 0 | CC t 0 p 0 m 0
    stripe 1: suppMis 0 maxNclasses 9 | ALT t 1439858 p 1439858 m 0 | LINE t 1880298 p 1880298 m 0 | CC t 0 p 0 m 0
    stripe 2: suppMis 0 maxNclasses 10 | ALT t 225822 p 225822 m 0 | LINE t 8136209 p 8136209 m 0 | CC t 6206077 p 6206077 m 0

that is, no stratum whose support differs, and 2 788 386, 11 466 029 and 6 206 077 strata efficient above, on and
below the switch line, with no disagreement. `fw.c` and `effcc_task1/fw.c` are byte-identical copies of `t2_fw.c` (not
deposited).

## `effcc_task1/effcc.c`: the closed form below the switch line

For every stratum with p_CC,CC = 1 (census code d_0 = 2), `effcc.c` evaluates the closed form of efficiency below the
switch line from integer error costs alone: with the closed classes of the self-play graph at ε = 0, K is the least
number of errors that lead from the class of mutual cooperation into another closed class, and E + m - 1 the cost of
the cheapest in-tree that leads the other m - 1 closed classes into it (a dynamic programme over sets of classes); the
form holds when K >= E + 2, and always when mutual cooperation is the only closed class. Its first argument is the
list of the strata in which the limit puts all weight on mutual cooperation (class C), its stripe and number of
stripes follow, and it counts the strata the form accepts and those on which it disagrees with the list. That list
was regenerated in a scratch directory for the original run and not kept; it is the class-C strata of the
support-based classes in `PRFull/release/census_class/census_support_class.u8.gz` (6 206 077 strata). From this
directory, on 2026-11-03,

    python3 -c "import gzip, numpy as np; b = np.frombuffer(gzip.open('../../../../PRFull/release/census_class/census_support_class.u8.gz').read(), np.uint8); np.savetxt('cc.txt', np.flatnonzero(b == 1), fmt='%d')"
    cc -O3 -o effcc effcc.c
    for s in 0 1 2 3 4 5 6 7 8 9 10 11; do ./effcc cc.txt $s 12 2> rerun_$s.log; cmp rerun_$s.log ver_$s.log; done   # 20 s

reproduced the twelve deposited logs of the original run byte for byte: the form accepts 6 206 077 strata (stripes
2, 5, 8 and 11; the others hold no stratum with p_CC,CC = 1), with no false positive and no false negative.

## `t2_altform2.c`: the characterisation above the switch line

`t2_altform2.c` reads the support of the ε → 0 limit of every stratum from the cheapest escapes of its closed classes
alone (a closed class is kept if it is recurrent in the graph of these escapes), without the arborescences of
`t2_fw.c`, and compares it, as `t2_fw.c` does, with the support that `freqs16` computes at one interior draw
(`suppMis`). It then tests the characterisation of efficiency above the switch line: a closed class inside the states
4 to 11 exists, and no closed class outside them can be reached without an error from the states 1 and 13, to which a
single error out of alternation leads (`altMis`, against the limit's support inside the states 4 to 11, `truthA`); and
a reachability rule for efficiency below the line (`ccMis`).

    cc -O3 -o t2_altform2 t2_altform2.c ../../cannon/m2f.c -lm
    for s in 0 1 2; do ./t2_altform2 $s 3 2> alt_$s.err & done; wait      # 25 s on 3 cores

The output of the original run was not kept. The run above, made on 2026-11-03, wrote

    stripe 0: suppMis 0 | altMis 0 (truthA 1122706) | ccMis 0 (truthC 0) truthL 1449522 | haveAlt 3188646 isolatedNonCore 0
    stripe 1: suppMis 0 | altMis 0 (truthA 1439858) | ccMis 0 (truthC 0) truthL 1880298 | haveAlt 3188646 isolatedNonCore 0
    stripe 2: suppMis 540 | altMis 0 (truthA 225822) | ccMis 40870 (truthC 6206077) truthL 8136209 | haveAlt 3188646 isolatedNonCore 0

(stripe 2 also lists the first ten strata of each kind of mismatch). The characterisation agrees with the limit on all
strata: 2 788 386 strata efficient above the switch line, no disagreement. The two other tests are cruder rules that
the paper does not use: the cheapest escapes alone give another support at 540 strata, where the arborescences of
`t2_fw.c` are needed, and the rule below the line errs at 40 870 strata (the paper's closed form there is that of
`effcc.c`).

## `verify_task1/`: the strata of W of dimension 8 to 15

Every stratum of dimension d = 8 to 15 is enumerated, its efficiency read at two interior draws, and each stratum whose
self-play is mutual cooperation (class C) is sampled: 512 draws (2048 for d = 12 to 15) where p15 is pinned to 0, the
only strata in which a friendly rival of W can lie, and 64 (256) elsewhere. Every draw is screened against 21
co-players (the 19 of `../witnesses.npy`, ALLD and ALLC); every draw that passes is scanned against all 65 536 binary
co-players. The programs import `../m2f.py`, the driver of `../../cannon/libm2f.dylib` (its build is in
`../README.md`), and ran on 12 cores of the author's laptop.

| run | stages | log |
|---|---|---|
| `python3 dim8W.py 8` | d = 8: enumeration and class; writes `dim8_C.npy`, the 491 127 class-C strata (3.9 MB, not deposited: the same run writes it again). Its screening stage did not finish: the log ends in `BrokenPipeError`s of the workers | `dim8.log` |
| `python3 dim8b.py 8 512 64` | d = 8: screen and scan of the strata of `dim8_C.npy` | `dim8b.log` |
| `python3 dimNW.py d 512 64`, d = 9, 10, 11 | all stages | `dimN.log` |
| `python3 dimNW.py d 2048 256`, d = 12 to 15 | all stages | not kept |

| d | strata | class C | points | confirmed friendly rivals of W |
|---|---|---|---|---|
| 8 | 3 294 720 | 491 127 | 42 617 792 | 0 |
| 9 | 1 464 320 | 212 767 | 16 116 928 | 0 |
| 10 | 512 512 | 70 542 | 4 867 264 | 0 |
| 11 | 139 776 | 17 526 | 1 149 440 | 0 |
| 12 to 15 | 34 112 | 3 569 | 917 248 | 0 |
| all | 5 445 440 | 795 531 | 65 668 672 | 0 |

The rows d = 8 to 11 are read from the logs. The run of d = 12 to 15 was repeated on 2026-11-03 with the pool of
`dimNW.py` reduced from 12 to 3 processes, which changes the draws (each worker's seed is fixed by its share) but not
the counts: class C 3152, 387, 29 and 1 strata; 810 496, 99 072, 7424 and 256 points; no draw passed the screen.

## `task2/`: the re-censuses of W and N

`task2/kit/` is the Cannon kit of the re-censuses, with its own README. On Cannon it ran as `~/PartnersRivalsM2f_task2`;
that directory is also archived whole, with its output, in `../../cannon/archive/cannon_2026-09-26.tgz`, and every file
deposited here is byte-identical to its copy in the archive. `t2.c` (mode `ncensus`, `ncensusP.sbatch`) re-censused N
and `t2b.c` (mode `ccensus`, `wcensus.sbatch`) W, each with up to 64 draws per candidate stratum and up to 4000 more
where the first ones found no friendly rival; the other arrays are probes, the stratified list and the test of the
census's seeding, and the re-census of E (`ncensusM.sbatch`), which was cancelled before it started. The seven
`*_JOBID` files hold the arrays 45785043, 45785044, 45785196, 45785345, 45785348, 45786256 and 45789243, the
re-censuses of CS Table 2. `compile.sh` builds `t2`; `t2b` is built in the same way from `t2b.c` and `m2f.c`.
On 2026-11-03 a build with `cc -O3` reproduced the archived lines of the strata 154 560 to 154 629 of W exactly
(`./t2b ccensus 1 154560 154630 64 4000 20260909`).

`task2/rec_W.txt` and `task2/rec_cannon.txt` collect the lines of `out/wc_*.txt` (W) and `out/ncP_*.txt` (N) of the
archived kit, one line per stratum with a confirmed friendly rival, `g nf h1 r1 h2 K2` (`t2b`) or `g nf h1 h2 K1 K2`
(`t2`): 364 and 5022 strata. Sorted, each is identical to the lines of the archived output.
