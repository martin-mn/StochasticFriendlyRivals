# `census_class/`: the efficiency class of every stratum (internal review, round 2, additional comment 3)

The census (`PartnersRivalsM2f/cannon/census.c`, Cannon array 45751662, 512 tasks) wrote one line `g cls nP nM K`
for every one of the 3^16 support strata in which at least one of its K draws was confirmed competitive against all
65 536 binary co-players. `cls` is the efficiency class read from the self-play weights at the stratum's first draw,
with the tolerance 1e-12 (C: w_CC > 1 - 1e-12; A: w_CD + w_DC > 1 - 1e-12; Z: w_DD < 1e-12; `-` otherwise).
`PartnersRivalsM2f/check/st_eff.txt` kept only the 1 524 791 lines with an efficient class.

The class is also decided by the support of the self-play limit alone (lemma L1 of SI section "Certificates": the
support follows from exact integer exponents and is constant on the open stratum), with no tolerance. The files here
give that class for every stratum and put it beside the census's own reading.

| file | contents |
|---|---|
| `census_st_45751662.tgz` | the stratum output of the census as the cluster wrote it: `st_0.txt` ... `st_511.txt` (3 254 063 lines, the 1 729 272 lines of class `-` included) and the task logs `st_45751662_<t>.out/.err`; harvested read-only from `cannon:~/CL1/PartnersRivalsM2f/ca1/out` on 2026-09-27 |
| `census_st_45751662.md5` | md5 of each of the 1536 files, taken on Cannon; every file in the archive matches |
| `census_support_class.u8.gz` | 43 046 721 bytes after gunzip: byte g is the support-based class of the stratum with census code g = sum_j d_j 3^j (d_j = 0: p_j = 0, 1: p_j free in (0,1), 2: p_j = 1); 0 none, 1 C, 2 A, 3 Z |
| `census_all_support_class.txt.gz` | every census line, tasks 0 to 511 in order, with the support-based class as a sixth column |
| `st_eff_support_class.txt.gz` | the same for the lines of `st_eff.txt` |
| `census_class.json` | the counts, the checks and the SHA-256 of the uncompressed files |
| `mk_census_class.py` | writes the class of all strata (compiles `PRFull/check/ubound/ub.c`, mode `classes`, in a temporary directory; about 30 s on 3 threads) and the join with `st_eff.txt` |
| `join_census.py` | checks the archive against the md5 list and joins the class to every census line (about 15 s) |

Checks, all passed: the support-based counts are those of the Methods (C 6 206 077, A 2 788 386, Z 2 471 566, none
31 580 692); every class-C stratum has p_0 = 1; all 512 task files are present, no stratum is listed twice, and each
task log holds its summary line (task 192's also holds 14 lines "srun: error: wrap_on_data ... Insane message length"
from Slurm; the task COMPLETED, `accounting/sacct_table7.psv`); on every one of the 3 254 063 census lines the class
the census read from the weights equals the support-based class (`-/-` 1 729 272, `C/C` 1 318 785, `A/A` 9 802, `Z/Z`
196 204; 0 mismatches).

Not here: the census's exemplar points (`pt_*.txt`, 570 MB), which are still only on Cannon; see `../RELEASE.md`.
