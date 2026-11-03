# cw.f vs ref_sf1.f: every change

`ref_sf1.f` is a byte copy of `DiskM1WF/Cannon/mw/sf1.f`
(sha256 `7d4649e312080253d1d8acd9019da54a5d5c1ef7ef22595796bd3a821063da7a`).
`diff ref_sf1.f cw.f` gives 48 hunks. They are listed below in order, with the spec item each one implements. In
`cw.f` every changed region carries a `[cw n]` comment:

| tag | spec item |
|---|---|
| hdr | header comment only |
| 1 | strategy storage: lineage slots (m1m2.F class table) |
| 2 | payoffs: ordered-pair cache keyed by lineage serials, 4-state GTH solve on a miss |
| 3 | mutation and initial population: the exact-alpha biased kernel |
| 4 | sampling: `iab(ic) += icnt(ic)` per sampled generation |
| 5 | classification: at slot death and at the end of the replicate; the invariant |
| 6 | outputs |
| 7 | `-DBINMODE` (the 16 corners, sf1's draws and payoffs, run through the machinery of 1, 2 and 4) |
| 8 | `-DITEND=` itend override, for tests only |

Everything the diff does not show is sf1.f unchanged: N, m, beta, u, eps (the `ie = 3` slot), nrep = nper = 10,
games.dat, the task decode, the interaction loop and its draw order, Fermi with `be = beta/(2m)` and
`thr = 25` and the sign-split formula, the RNG (`mtfill` into `rb(8192)`, `init_genrand` per replicate), the
per-individual `ef` loop, `pay = ef*rrn*rrm/z` and `Emax = max(R,(S+T)/2)`.

## Hunks

`ref` = line range in ref_sf1.f, `cw` = line range in cw.f.

| # | ref | cw | tag | change |
|---|---|---|---|---|
| 1 | 1 | 1-2 | hdr | title line: "CONTINUOUS space" |
| 2 | 3-7 | 4-21 | hdr | what cw is, and the list of `[cw 1..8]` items |
| 3 | 10-15 | 24-27 | hdr | convention paragraph restated for the effective vector s in [eps,1-eps]^4. The 16 corners are sf1's 16 strategies |
| 4 | 17-25 | 29-32 | hdr | sf1's payoff-table paragraph replaced by a pointer to ref_sf1.f's header; seed bases stated |
| 5 | 27-66 | 34-46 | hdr, 8 | sf1's history paragraphs dropped (still readable in ref_sf1.f); new usage/output description; `#define ITEND 10000000` default |
| 6 | 81 | 62-63 | 7 | `#ifdef BINMODE` around `iseedb=10000000` |
| 7 | 82 | 65-67 | 7 | continuum `iseedb=40000000` |
| 8 | 97 | 83-88 | 1, 2 | `maxc=2*n+2`; `ktbig=2^31` (integer*8) for the pair key |
| 9 | 104-105 | 95-98 | 8, 1 | `data itv /10000000/` becomes `/ITEND/`. `s(n),sn(n),si,sj` removed from the integer declaration (they become `icls/iclsn`) |
| 10 | 112 | 106-128 | 1, 2, 4, 5 | declarations: `icls, iclsn, icnt, ifree, istg, gid, gcount, key, strat(4,maxc), sstr(4,n), kid, kstg`; `pkey, pc1, pc2`; `iab`; `atomcnt, strcnt, corncnt, ncnt, atot, strtot, corntot, nexp, nsum, nsolve, nlook, icatom, icnear, sv, share` |
| 11 | 115 | 132-135 | 1, 2 | `common /bigcache/`, `/bigtab/` (the large tables live in COMMON, as in m1m2.F) |
| 12 | 193 | 214-223 | 5 | classifier game coordinates `gamu = cS`, `gamv = cT-1` from the exact columns. Asserts `cR = 1` and warns if the display columns differ (they do not) |
| 13 | 198 | 229-237 | 6 | the w-file header: sf1's lines under BINMODE, the continuum diagnostics header otherwise |
| 14 | 207 | 247-248 | 6 | `#ifdef BINMODE` before sf1's `# code count` h-header line |
| 15 | 208 | 250 | 6 | `#endif` |
| 16 | 222 | 265-266 | 7 | comment: the vlt/vht table is BINMODE-only |
| 17 | 229 | 274 | 7 | `#ifdef BINMODE` around sf1's 16x16 table build |
| 18 | 238 | 284-298 | 5 | `#endif`; the per-task corner cache begins: 16 corners classified once, returned stratum checked |
| 19 | 239 | 300-301 | 5 | (corner cache, continued: sf1's `enddo` lines interleave in the diff) |
| 20 | 240 | 303-315 | 5 | task totals zeroed; the continuum h-header gets `# corner atoms (code 0..15)` |
| 21 | 241 | 317-326 | 2 | pair cache emptied once per task; `gcount = 0`. gid increases for the whole task, so no per-replicate clearing is needed |
| 22 | 255 | 341-360 | 1, 4, 5 | per-replicate reset of `icnt, iab`, class counts, `sdiv, nsolve, nlook` |
| 23 | 258 | 364-365 | 1, 3 | comment: the continuum's initial population is n kernel draws, one slot each |
| 24 | 260 | 368-370 | 7 | `#ifdef BINMODE` and comment around sf1's 4 fair coins |
| 25 | 270 | 380-391 | 1, 3, 7 | `s(i)=is` becomes `kid(i)=is` plus the corner's pmin/pmax vector (BINMODE), or `call umut4` (continuum); then `gid(i)`, `icls(i)=i`, `icnt(i)=1` |
| 26 | 271 | 393-395 | 1 | `nhi = nlive = nhw = n` |
| 27 | 284 | 408-409 | 1 | `si=s(i)` becomes `ia=icls(i)` |
| 28 | 293 | 418-437 | 2, 7 | `sj=s(j)` becomes `ib=icls(j)`, plus the pair-cache lookup. On a miss: BINMODE reads `vlt/vht(kid(ia),kid(ib))`, the continuum calls `pay4` and forms sf1's two lines `rr*w0+ss*w1+tt*w2`, `rr*w0+tt*w1+ss*w2` |
| 29 | 295-296 | 439-441 | 2 | `f(i)+=vlt(si,sj)`, `f(j)+=vht(si,sj)` become `pc1(ia,ib)`, `pc2(ia,ib)`. Same order of accumulation |
| 30 | 307-308 | 451 | 4 | per-individual `hist(s(i))+=1` removed; the `efg` loop is otherwise sf1's |
| 31 | 310 | 454-460 | 4 | `iab(is) += icnt(is)` over live slots; `sdiv += nlive` |
| 32 | 313 | 464-465 | 1 | comment: mutants are staged |
| 33 | 314 | 467-468 | 1 | `nstg = 0` |
| 34 | 325-344 | 479-487 | 1, 7 | the mutant is staged (`istg`, `iclsn(i) = -nstg`). BINMODE keeps sf1's `imut = 1` single draw `int(rb*16)` into `kstg`; sf1's dead `imut /= 1` branch is not carried |
| 35 | 345 | 489-494 | 3 | continuum: `call umut4(rb,irb,sstr(1,nstg))` |
| 36 | 380 | 530 | 1 | comment |
| 37 | 382 | 532 | 1 | `sn(i)=s(j)` becomes `iclsn(i)=icls(j)` |
| 38 | 384 | 534 | 1 | `sn(i)=s(i)` becomes `iclsn(i)=icls(i)` |
| 39 | 391-392 | 541-547 | 1, 5 | `s(i)=sn(i)` is replaced by m1m2.F's new-generation step, which begins with the recount |
| 40 | 393 | 549-552 | 1 | (recount, continued) |
| 41 | 394 | 554-611 | 1, 5 | recycle dead slots. A slot that dies with `iab > 0` is classified before reuse (BINMODE also adds `iab` into `hist(kid)`). The staged mutants are then placed: free slots first, then above `nhi`, with an overflow stop and a check that the slot has `iab = 0`. Then `nhi`, `nlive`, the high-water mark `nhw`, and `icls = iclsn` |
| 42 | 396 | 614-643 | 5 | end of replicate: every slot with `iab > 0` is classified. Asserted: `sum(atomcnt) == sum(strcnt) == (itend-init)*n` and `gcount < 2^31` |
| 43 | 403 | 651-653 | 6, 7 | `#ifdef BINMODE` around sf1's output block, which is otherwise verbatim (unit-3 line, winners, w-file) |
| 44 | 467 | 718-750 | 5, 6, 7 | BINMODE diagnostics on stdout. Continuum: the unit-3 line `gid ie u v Emax pay ef meanlive solves/gen share(0..7)` and the w-file diagnostics line; per-task totals accumulated (both modes) |
| 45 | 477 | 761-762 | 6, 7 | `#ifdef BINMODE` around sf1's h-file tail |
| 46 | 485 | 771-788 | 6 | continuum h-file tail: `# atoms 8`, `# strata 81`, `# corners 16` sections of `index count` lines, then `# total N` |
| 47 | 488 | 792-830 | 3 | `#endif`, then new subroutine `umut4`: m1m2.F's exact-alpha kernel on sf1's buffer |
| 48 | 489 | 832-982 | 2, 5 | new subroutines `pay4` (4-state GTH), `clsab` (stratum, corner cache or `m1cls`, counts) and `cwcls1` (the single `bind(C)` call site of `m1cls`) |

Rows 19, 36, 40 and a few others are one-to-few-line hunks. They exist only because unchanged sf1 lines split a
tagged region in the diff.

## Behavioural differences from sf1 that are *not* bugs

* **BINMODE stdout.** BINMODE writes one `BIN ...` diagnostics line per replicate to stdout. sf1 writes
  nothing there. The three output files are byte-identical to sf1's (see README, V3).
* **The dead `imut /= 1` branch** (sf1's legacy 4-uniform mutation) is gone. `imut = 1` is a parameter in sf1,
  so that branch never ran.
* **The continuum w-file** is diagnostics, not sf1's winner list: a continuum strategy has no finite winner table.
* **`mid` and `payf3`** are still built in the continuum. They are used only by BINMODE's w-file. They cost 16 integer ops.

## Integration revision (2026-09-24, after the classifier was linked)

The hunk table above was made against the pre-integration `cw.f` (sha256 `254b4e9a...12e9d6`). The integrated
`cw.f` (sha256 `37ec579b248ffb57716aebfa4d96f871d3852eaf5810874874baba4c729a2ebb`) gives **47** hunks against
`ref_sf1.f`, because two new header lines merged two old hunks. Its line numbers are shifted by the edits below,
and nothing else changed. `diff` of the two `cw.f` versions shows exactly these edits:

| where | tag | change |
|---|---|---|
| header | 9 | a line for `[cw 9] -DCLSREP=<k>` |
| `ncnt` declaration and reset | 5 | `ncnt(4)` becomes `ncnt(6)`. The new entries count the **exact ties** (m1cls neartie = 2) and their abundance, separately from the **near ties** (neartie = 1) |
| `clsab` | 5 | `inr.ne.0` becomes `inr.eq.1` for the near-tie counter and `inr.eq.2` for the new exact-tie counter. Any other value stops the run |
| continuum w-file header and line | 6 | two more fields, `nexacttie exacttie_abund` |
| BINMODE stdout line | 7 | prints the 6 counters (its format widened to match). This is stdout only, so sf1's files are unaffected (V3) |
| `cwcls1` | 9 | under `-DCLSREP=k` only, the `m1cls` call is repeated k-1 extra times, for timing. Outputs are unchanged (checked: `cwr.x` output = `cwt.x` output, byte for byte) |

**Why split the counter.** The classifier returns 2 when the Nash margin is exactly 0.0. That is a structural
tie, such as pi(ALLC,p) = pi(p,p) = 1 on the slab p1 = 1, decided exactly by convention A's weak inequality. It
is common: 71 % of all sampled abundance in the 512-cell local run. It returns 1 only for a margin that is
rounding-sensitive (0 < |m| <= 1e-9 scale). Before this edit, one counter mixed the two, and the rare type 1
was invisible under the common type 2.

## Run F2 (2026-09-25): cw2.f = mkcw2.py(cw.f), and the kit changes that carry it

This is run F2, the process of `DiskM1WF/Cannon/mw/sf2.f`
over the continuous space. `ref_sf2.f` is a byte copy of sf2.f (sha256 `e227b4c8472ace1d7727481632a96a025c3d5d045bef402b6c849e0cf2061edf`).
sf2.f differs from sf1.f in exactly five parameter lines: n, iseedb, uval, bval, itv.

**cw2.f is never edited by hand.** `mkcw2.py` generates it from cw.f (sha256 asserted: `37ec579b...c729a2ebb`). It
applies 13 substitutions, each asserted to match exactly once, and asserts the diff size: 14 changed lines, 7 code
and 7 comment, with the line count unchanged. `diff cw.f cw2.f` is the audit.

| cw.f line | tag | change |
|---|---|---|
| 20 | hdr | comment: BINMODE bit-identical to `sf2.x` (was sf1.x) |
| 21 | hdr, 8 | comment: itend default `sf2's 1e8` (was sf1's 1e7) |
| 33-34 | hdr | comment: continuum seeds `42000000 + ...`; BINMODE keeps `sf2's 11000000`, compared with `sf2.x` |
| 44 | hdr | comment: BINMODE output = `sf2.x`'s three files |
| 47 | 8 | `#define ITEND 100000000` (sf2's `itv`) |
| 56 | — | `parameter (n=100,m=1)` (sf2's n) |
| 64 | 7 | comment: seed base `sf2's` in BINMODE |
| 66 | 7 | BINMODE `iseedb=11000000` (sf2's) |
| 68 | 7 | continuum `iseedb=42000000`: disjoint from F1 40M, kappa 41M, and 6M 8M 9M 10M 11M 12M 13M 20-25M 30-32M 39M |
| 94 | — | `data uval /1.d-04/` (sf2's) |
| 95 | — | `data bval /3.d+00/` (sf2's) |
| 97 | 8 | comment: `the literal 100000000 is now the macro ITEND (default 1e8)` |
| 240 | 6 | w-file header string `(maxc 202)` (was 2002): maxc = 2n+2 |

Nothing else differs. eps, nrep = nper = 10, games.dat, the kernel, the classifier, the task decode, the loop,
the RNG discipline and the output formats are F1's. The kernel measure `kappa` is shared.

**Kit changes** (none to cw.f, task.sh, prod.slurm, probe.slurm, kappa.slurm, games.dat, m1cls.c, m1cls.h, kappa.c,
ref_sf1.f, mt.f, mtb.f, payf2.f or payf3.f; their sha256 were verified unchanged at the end of this revision):

| file | change |
|---|---|
| `mkcw2.py` | new. `--check` writes nothing and exits 1 unless cw2.f is exactly its output |
| `cw2.f` `ref_sf2.f` | new (generated / byte copy) |
| `compile.sh` | new targets `cw2 cwb2 cwp2 cwt2 cwbt2 sf2t`; groups `f2`, `probe2` (ITP2 = 1e7) and `test2`; every F2 target first runs `mkcw2.py --check`; the stub refusal also lists cw2, cwb2 and cwp2. The only removed line is the old stub-refusal `case`. All 9 old targets rebuild byte-identically, and the default set is still `cw cwb kappa` |
| `probe2.slurm` `prod2.slurm` | new. test / sapphire, `PartnersRivalsM1f.cwprobe2` / `.cwprod2`, logs `{o,e}_f2probe_%a` / `{o,e}_f2prod_%a`, calling the unchanged `task.sh` with `probe2 cwp2.x` / `run2 cw2.x` |
| `submit.sh` | tags `probe2` and `prod2`, with their own JOBID_ and IDX_ files and ledger notes naming F2 and N=100 beta=3 u=1e-4 itend=1e8. For F2 tags it refuses a binary whose strings lack `maxc 202)`, and a cw2.f that fails `mkcw2.py --check`. The plain `JOBID` is still written only for `prod` (F1's 48366968). The RESUBMIT and QUEUE_OK guards are unchanged and apply to F2 |
| `verify.sh` | tags `probe2` (probe2/, logs f2probe) and `prod2` (run2/, logs f2prod). The probe2 rate uses 1e8 generations per task and names prod2.slurm. The F1 output is byte-identical (tested with fake sacct) |
| `pack.py` | `--run=F1` (the default) or `--run=F2`. F2 writes `cw2_e4.tsv`, `cw2_e4.rep.tsv` and `cw2_e4.diag.tsv`. It shares `kappa_e4.tsv`: it writes the file only if absent, and otherwise requires identical data lines. Every check is kept, parametrised by N, itend, umut and beta. The F1 outputs are byte-identical to the pre-F2 pack.py's (all four files `cmp`-equal) |
| `check_v4.py` | optional `--n N` (default 1000), so that the invariant uses N and the table bound is 2N+2 |
| `README.md` | status line and section 11 (F2) |

**The rule while 48366968 (F1 production) is live: no relink of cw.x on Cannon.** Build only `bash compile.sh f2`
and `bash compile.sh probe2` in the remote kit. Never run plain `bash compile.sh` or any F1 target there.

## Post-audit fixes, 2026-09-25 (main session, before any F2 submission)

| file | change | why |
|---|---|---|
| `submit.sh` | the F2 build guard matches the fragment `axc 202)` and requires `axc 2002)` to be absent | ifx stores the w-file header's string literals **interleaved** (`...high water (m` … `axc 202); nearti...`), so the original `grep 'maxc 202)'` refused a correct ifx build of cw2.f; it had only been tested on gfortran binaries. Tested on Cannon: cwp2.x, cw2.x accepted; cwp.x, cw.x, cwb.x refused |
| `pack.py` | F2 writes the SHARED `kappa_e4.tsv` only from a complete, consistent pack; the "identical, left untouched" message only when it is | audit F2 nonblocking 2 and 3: a partial F2 pack left a 511-cell kappa file that would block a later clean pack; the message printed even on a mismatch. F1 outputs re-checked byte-identical |
| `mkcw2.py` | one more comment edit: cw2.f's `[cw 7]` header names sf2's process (7 code + 8 comment lines) | audit F2 nonblocking 6 |
| `README.md`, `CHANGES.md`, `probe2.slurm` | cw.f's sha256 abbreviation (`...c729a2ebb`, not task.sh's `...5ebd74`); seed bases 6M 12M 13M added to the lists; probe2's header names the maxc string edit | audit F2 nonblocking 1, 7, 9 |
| `../../Figures/fig4c.py` | row labels carry N (`Memory-1 continuous, $N = 100$`), as the originals do | audit F2 nonblocking 4: SIFig4c was indistinguishable from Fig4c on the page |

Cannon-side checks with the Intel build (the agents' checks were gfortran-only): BINMODE-F2 vs sf2t at isl 243, 1203,
itend 2e5: 6/6 files identical. F1: BINMODE vs sf1t 6/6 identical, and the linked classifier's 16 corner atoms at all
512 cells 0/8192 mismatches against DiskM1WF/Data/m1_atoms_512.npz.
