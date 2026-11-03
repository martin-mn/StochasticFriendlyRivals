#!/usr/bin/env python3
"""assemble_release.py -- gather the code and data release of "Friendly rivals exist but rarely evolve among stochastic strategies of direct reciprocity" into
one tree, and write its SHA-256 manifest.

    python3 PRFull/release/assemble_release.py --dry-run            list what would be copied (per item and per
                                                                     directory, with sizes and the total) and every
                                                                     promised item that cannot be found; copies nothing
    python3 PRFull/release/assemble_release.py --dry-run --files     the same, and every file with its size
    python3 PRFull/release/assemble_release.py --dest DIR            copy into DIR (which must not exist yet) and write
                                                                     DIR/MANIFEST.sha256 and DIR/CONTENTS.tsv; refused,
                                                                     before anything is copied, while the manuscript
                                                                     guard flags a file (or, in the author's tree,
                                                                     finds no manuscript to compare with)
    python3 PRFull/release/assemble_release.py --harvest-cannon      copy, read-only, the census output that existed only
                                                                     on Cannon (cannon:~/CL1/PartnersRivalsM2f/ca1/out:
                                                                     exemplar points, dimension probe, sweep; about
                                                                     840 MB) into PRFull/release/harvest/, and check
                                                                     every file's md5 against md5sum run on Cannon.
                                                                     Nothing on Cannon is changed.  Run it only when
                                                                     Martin has agreed to the transfer.  (Done on
                                                                     2026-09-27: 2440 files, every md5 matches; the
                                                                     release puts them in PartnersRivalsM2f/cannon/
                                                                     ca1_out/, with the md5 list ca1_out.md5 beside.)

The release keeps the author's layout: the directory that holds PRFull/ becomes the release root, with PRFull/ and,
read-only and in part, the projects whose programs and data the paper uses (PartnersRivalsM1f, PartnersRivalsM2f,
PartnersRivals1, PartnersRivals, DiskM1WF, DiskM2WF).  Nothing is changed in any of them: files are only read.  The
promises are those of the paper's "Code and data availability" statement, SI section "Implementations,
conventions and cost" and the Computational Supplement (Text/cs.tex); ITEMS below lists, for each promise, where its files are.  Three sources decide what is taken
from the other projects: the rules of ITEMS (whole directories of programs and their output), the traced inputs of
the figure scripts (figure_inputs.json, written by reproduce_figures.sh --record), and every file named by an
absolute path in an included program (found by scanning them).  Left out everywhere: compiled code (Mach-O and ELF
files: the release ships the sources and build scripts), __pycache__, .DS_Store, LaTeX build files, working copies
(PRFull/release/work), the internal referee reports (PRFull/Referee), the scratch of the M2f derivation pass
(PartnersRivalsM2f/check/working, 196 MB), except the programs and outputs of five checks the paper describes that
ran there (CHECKS, taken by the item checks), PartnersRivalsM1f/cannon/cw/audit (synthetic
runs with made-up job IDs 60020849 and 60032178; they must not be deposited), virtual environments and user-level
package installs (venv/, pylib/), PRFull/cannon/evoaudit/data (643 MB: a local extraction of the strategy dumps in
PartnersRivalsM2f/cannon/archive/*.tgz, which are deposited; reproduce_certificates.sh extracts them again), the
rebuilt intermediates of PRFull/cannon/stratcert/prep (cls.bin, c2.idx, c2.qm, c2.qp, 173 MB; its README says how to
regenerate them in 2.5 min) and the scratch directory PRFull/cannon/sweep13/validate/build (the outputs of
validate_local.sh and rcheck.sh).  The harvested census output in PRFull/release/harvest is taken only by the items that
name it (census, sweeps), and goes to PartnersRivalsM2f/cannon/ca1_out in the release.

Release 1.0.0 holds the code and data; the paper's sources are added in release 1.1.0, when the article is published (run with --with-paper).
Without that option the paper's sources and PDFs (Text/ms.tex, refs.tex, reorder_bib.py, ms.pdf; the Computational
Supplement Text/cs.tex, cs_bib.py, cs.pdf), SI Table 4 and SI Box 3 with their legends, and their LaTeX by-products are
left out (PAPER_EXCLUDE; the two are deposited without their legends, as Text/si_table4_body.tex and
check/ref3/eq3_boundary/box/si_box3_strata_body.tex).  Manuscripts and
internal notes are not part of any release (MS_EXCLUDE), and a guard (tools/manuscript_guard.py) compares every
selected text file with the manuscripts listed in a local file and flags any copy of their text.  Where a project README carries text of the paper or its drafts, or internal notes,
the release ships a release-only copy of it from PRFull/release/overrides/<same path> (the project itself is not
changed).

Top level of the release: README.md (= PRFull/release/RELEASE.md), LICENSE, reproduce_all.sh (calls
PRFull/release/reproduce_all.sh), MANIFEST.sha256 (sha256sum format: "hash  path", every other file of the tree,
`shasum -a 256 -c MANIFEST.sha256` checks it), CONTENTS.tsv (path, bytes, item).
"""
import argparse
import collections
import fnmatch
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PRF = os.path.dirname(HERE)
ROOT = os.path.dirname(PRF)
AUTHOR_ROOT = "/Users/martin/Documents/CL1"
HARVEST = os.path.join(HERE, "harvest")
OVERRIDES = os.path.join(HERE, "overrides")   # release-only copies: overrides/<path> is shipped in place of <path>

# ------------------------------------------------------------------------------------------------ what is left out
EXCLUDE_ANY = ["*/.DS_Store", "*/__pycache__/*", "*.pyc", "*/.pycache/*", "*.o", "*.dSYM/*", "*/.git/*", "*~", "*.swp",
               "PRFull/release/work/*", "PRFull/Referee/*", "*/venv/*", "*/pylib/*",
               "PRFull/cannon/evoaudit/data/*", "PRFull/cannon/stratcert/prep/cls.bin", "PRFull/cannon/stratcert/prep/c2.idx",
               "PRFull/cannon/stratcert/prep/c2.qm", "PRFull/cannon/stratcert/prep/c2.qp",
               "PRFull/cannon/sweep13/validate/build/*",
               "PRFull/Text/ms.aux", "PRFull/Text/ms.log", "PRFull/Text/ms.out", "PRFull/Text/ms.fls",
               "PRFull/Text/ms.fdb_latexmk", "PRFull/Text/ms.synctex.gz", "*.bak",
               "PRFull/Text/cs.aux", "PRFull/Text/cs.log", "PRFull/Text/cs.out", "PRFull/Text/cs.fls",
               "PRFull/Text/cs.fdb_latexmk", "PRFull/Text/cs.synctex.gz",
               "PartnersRivalsM2f/check/working/*", "PartnersRivalsM1f/cannon/cw/audit/*"]
# manuscripts and internal notes are not part of any release: the companion paper's manuscript, the drafts of the merged projects,
# the agents' derivation notes and manuscript drafts of M2f (no program reads them), the project notes and dated
# backups of PRFull, PRFull/README.md (the author's working notes; the release's README is RELEASE.md), a third-party
# article kept as a reference text, and the sentences fig1data.py --verify-sources checks in the manuscripts
MS_EXCLUDE = ["PartnersRivals1/Revised2/*", "PartnersRivals1/Submitted/*", "PartnersRivals1/Text/*",
              "PartnersRivalsM1f/Text/*", "PartnersRivalsM2f/Text/*", "PartnersRivalsM2f/check/derivations/*",
              "PRFull/Text/old/*", "PRFull/notes/*", "PRFull/README.md", "PRFull/Referee/*",
              "PartnersRivalsM1f/check/working/ref/*", "PRFull/Figures/fig1_source_quotes.json",
              "PRFull/release/tools/manuscript_sources.json", "PRFull/release/overrides/*",
              "PRFull/release/downloads/*"]
# the paper and the Computational Supplement: left out of release 1.0.0, added by --with-paper (release 1.1.0)
PAPER_EXCLUDE = ["PRFull/Text/ms.*", "PRFull/Text/cs.*", "PRFull/Text/refs.tex", "PRFull/Text/reorder_bib.py",
                 "PRFull/Text/cs_bib.py", "PRFull/Text/*.pdf", "PRFull/Text/pdflatex*",
                 "PRFull/Text/si_table4.tex", "PRFull/Figures/sitable4_legend.tex",   # SI Table 4 with its legend
                 "PRFull/Text/si_box3_strata.tex", "PRFull/Text/si_box3_legend.tex",   # SI Box 3 with its legend
                 "PRFull/check/ref3/eq3_boundary/box/si_box3_strata.tex"]
PAPER_ITEMS = ("paper", "supplement")
WITH_PAPER = False                              # set by --with-paper
GUARD_WORDS = 60                                # a selected file sharing this many words (in runs of >= 10) with a manuscript is flagged
COMPANION_RUN = 15                              # ... and a file with a run of this many words from the companion paper's manuscript
# files above GUARD_WORDS that were read and found to share only labels, formulas or descriptions of this release (not
# manuscript prose); each is accepted up to the count found on review (2026-10-17; the last two 2026-10-26), and flagged
# again if it grows
GUARD_REVIEWED = {
    "PRFull/Figures/sitable4.py": (70, "the row labels of SI Table 4 (the table itself is deposited, without its legend)"),
    "PRFull/release/RELEASE.md": (140, "the release's own description and title, which the paper's availability statement repeats"),
    "PRFull/release/assemble_release.py": (110, "the title and the item descriptions, which the availability statement repeats"),
    "PRFull/check/cert/witnesses/index.json": (80, "one label, 'the eight binary friendly rivals 111*000101*1*010 of', 8 times"),
    "PRFull/Text/si_table4_body.tex": (120, "the row labels of SI Table 4 (the table is deposited without its legend)"),
    "PRFull/Figures/sitable4.out": (170, "the row labels of SI Table 4, twice (runs F1 and F2), with their values"),
}
EXCLUDE_ANY += MS_EXCLUDE
# the programs and outputs of five checks that the paper and the Computational Supplement describe and that ran in the
# scratch of the M2f derivation pass, which is otherwise left out (added 2026-11-03; the README beside them is new):
# t2_fw.c, the efficient sets of all 3^16 strata by minimum arborescences (CS Table 1); effcc.c with its logs and
# t2_altform2.c, the closed forms of efficiency below and above the switch line on all strata; the second further
# campaign in W, every stratum of dimension 8 to 15 (CS Section 2); and the kit of the re-censuses of the Methods (t2.c,
# t2b.c), as it ran on Cannon (~/PartnersRivalsM2f_task2, also inside
# PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz), with the collected lines of its W and N arrays. Taken one by
# one, by the item checks; the rest of working/ stays out.
WORKING = "PartnersRivalsM2f/check/working/"
CHECKS = [WORKING + f for f in (
    "README.md", "t2_fw.c", "t2_altform2.c", "effcc_task1/effcc.c")] + [
    WORKING + "effcc_task1/ver_%d.log" % s for s in range(12)] + [WORKING + f for f in (
    "verify_task1/dim8W.py", "verify_task1/dim8b.py", "verify_task1/dimNW.py",
    "verify_task1/dim8.log", "verify_task1/dim8b.log", "verify_task1/dimN.log",
    "task2/rec_W.txt", "task2/rec_cannon.txt")] + [WORKING + "task2/kit/" + f for f in (
    "README.md", "compile.sh", "task.sh", "t2.c", "t2b.c", "m2f.c", "witnesses.h", "list.txt",
    "nprobe.sbatch", "list.sbatch", "seed.sbatch", "ncensus.sbatch", "ncensusP.sbatch", "ncensusM.sbatch",
    "cprobe.sbatch", "wcensus.sbatch", "NPROBE_JOBID", "LIST_JOBID", "SEED_JOBID", "NCP_JOBID", "NCM_JOBID",
    "CPROBE_JOBID", "WC_JOBID")]
WALK_SKIP = ["PRFull/release/harvest/*"]      # not taken by walking a tree: only by the items that name these files
MAGIC = (b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xfe\xed\xfa\xce", b"\xca\xfe\xba\xbe", b"\x7fELF")
PROGRAM_EXT = (".py", ".sh", ".sbatch", ".slurm", ".c", ".h", ".f", ".f90")


# ------------------------------------------------------------------------------------------------ the promises
class Item:
    def __init__(self, key, promise, trees=(), files=(), globs=(), require=(), missing=()):
        self.key, self.promise = key, promise
        self.trees = list(trees)          # directories taken whole (relative to ROOT), minus the exclusions
        self.files = list(files)          # single files that must exist
        self.globs = list(globs)          # glob patterns (relative to ROOT)
        self.require = list(require)      # (glob, n_min, what): must match at least n_min files, else flagged
        self.missing = list(missing)      # promised files known not to be in the tree: (what, where, size)


CANNON_CA1 = "cannon:~/CL1/PartnersRivalsM2f/ca1/out"
ITEMS = [
    Item("census",
         "the census output for every stratum, with the confirmed friendly-rival points at full double precision; "
         "and the support-based efficiency class of every stratum (referee 2, additional comment 3)",
         trees=["PRFull/release/census_class"],
         files=["PartnersRivalsM2f/check/st_eff.txt", "PartnersRivalsM2f/check/census.py", "PartnersRivalsM2f/check/census.out",
                "PartnersRivalsM2f/check/dim.out", "PRFull/release/harvest/ca1_out.md5"],
         globs=["PRFull/release/harvest/ca1_out/pt_*.txt", "PRFull/release/harvest/ca1_out/probe_*"],
         require=[("PRFull/release/census_class/census_st_45751662.tgz", 1, "the stratum lines of all 512 census tasks"),
                  ("PRFull/release/census_class/census_support_class.u8.gz", 1, "the support-based class of all 3^16 strata"),
                  ("PRFull/release/census_class/census_all_support_class.txt.gz", 1, "every census line with its support-based class"),
                  ("PRFull/release/harvest/ca1_out/pt_*.txt", 512, "the confirmed points at full double precision"),
                  ("PRFull/release/harvest/ca1_out/probe_45749304*", 1, "the census probe output")],
         missing=[("the confirmed friendly-rival points of the census at full double precision (pt_0..511.txt, up to 4 "
                   "points per stratum and sign)", CANNON_CA1 + "/pt_*.txt", "827 MB"),
                  ("the census probe output (array 45749304)", CANNON_CA1 + "/probe_45749304*", "small")]),
    Item("sweeps",
         "the output of the sweeps (the predicate-free sweep, the dimension probe, the rarity sweeps at alpha = 0.333 "
         "and 1/3 and the alpha grid, the memory-one sweeps)",
         files=["PartnersRivalsM2f/check/sweep.out", "PartnersRivalsM2f/check/sweep_a333.tsv"],
         trees=["PRFull/cannon/sweep13/out_13", "PartnersRivalsM1f/cannon/out"],
         globs=["PRFull/release/harvest/ca1_out/sw_*", "PRFull/release/harvest/ca1_out/dim_*"],
         require=[("PRFull/cannon/sweep13/out_13/a_1o3_*.tsv", 512, "the 512 tasks of the sweep at alpha = 1/3"),
                  ("PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz", 1, "the archived alpha = 0.333 sweep and alpha grid"),
                  ("PRFull/release/harvest/ca1_out/sw_*", 256, "the raw output of the predicate-free sweep"),
                  ("PRFull/release/harvest/ca1_out/dim_*", 384, "the raw output of the dimension probe")],
         missing=[("the raw output of the predicate-free sweep (array 45751663, 256 tasks; summary in check/sweep.out)",
                   CANNON_CA1 + "/sw_*", "0.3 MB"),
                  ("the raw output of the dimension probe (array 45752239, 384 tasks; summary in check/dim.out)",
                   CANNON_CA1 + "/dim_*", "15 MB")]),
    Item("evolution",
         "the per-replicate evolutionary output, with the sampled strategies",
         trees=["PartnersRivalsM2f/check/wfc_a333", "PartnersRivalsM2f/check/wfc_f2", "PartnersRivalsM2f/cannon/archive",
                "PartnersRivalsM1f/cannon/cw/run1", "PartnersRivalsM1f/cannon/cw/run2", "PartnersRivalsM1f/cannon/cw/packed",
                "PartnersRivalsM1f/cannon/cw/kappa", "PRFull/cannon/cw2/out", "PRFull/cannon/wfc2/out"],
         require=[("PartnersRivalsM2f/check/wfc_a333/r_1o3_*.txt", 512, "memory two, N = 1000: per-replicate output"),
                  ("PartnersRivalsM2f/check/wfc_f2/r_1o3_*.txt", 512, "memory two, N = 100: per-replicate output"),
                  ("PartnersRivalsM2f/cannon/archive/cannon_out_f2_2026-09-26.tgz", 1, "memory two, N = 100: raw output with the sampled strategies"),
                  ("PartnersRivalsM1f/cannon/cw/packed/cw1_e4*.tsv", 1, "memory one, N = 1000"),
                  ("PartnersRivalsM1f/cannon/cw/packed/cw2_e4*.tsv", 1, "memory one, N = 100")]),
    Item("audit",
         "the exact audit of a sample of the atoms (SI Table 6), with its sample and its output; the strategy dumps it "
         "samples are those of the evolution item (the archives), not duplicated",
         trees=["PRFull/cannon/evoaudit"],
         require=[("PRFull/cannon/evoaudit/out_probe/chunk_*.jsonl", 4, "the probe's chunks (array 48975385)"),
                  ("PRFull/cannon/evoaudit/out/chunk_*.jsonl", 124, "the full run's chunks (array 48981982)"),
                  ("PRFull/cannon/evoaudit/out*/chunk_*.done", 128, "the completion markers of the 128 chunks"),
                  ("PRFull/cannon/evoaudit/sample/sample.npy", 1, "the sample of 27 261 strategies"),
                  ("PRFull/cannon/evoaudit/sample/design.json", 1, "the design of the sample"),
                  ("PRFull/cannon/evoaudit/results/summary.json", 1, "the audit's estimates"),
                  ("PRFull/cannon/evoaudit/results/confusion.txt", 1, "the confusion tables"),
                  ("PRFull/cannon/evoaudit/results/records.tsv.gz", 1, "one line per audited strategy"),
                  ("PRFull/cannon/evoaudit/results/analyze.out", 1, "the output of analyze.py"),
                  ("PartnersRivalsM2f/cannon/archive/cannon_2026-09-26.tgz", 1, "the strategy dumps it samples, N = 1000"),
                  ("PartnersRivalsM2f/cannon/archive/cannon_out_f2_2026-09-26.tgz", 1, "the strategy dumps it samples, N = 100")]),
    Item("sidata",
         "the values of the evolutionary maps at every game (SI Data 1 to 3)",
         files=["PRFull/Text/SIData1_games512.csv", "PRFull/Text/SIData2_sensitivity.csv", "PRFull/Text/SIData3_runlength.csv"],
         require=[("PRFull/Figures/mk_sidata1.py", 1, "the program that writes SI Data 1"),
                  ("PRFull/release/tools/mk_sidata23.py", 1, "the program that writes SI Data 2 and 3 (written for the "
                   "release; it reproduces the deposited files byte for byte)")]),
    Item("certificates",
         "the exact certificates of the friendly rivals of the text and of SI Table 2, and the certificates of the SI "
         "section Certificates with the programs that check them; the 25 witnesses",
         trees=["PRFull/check", "PRFull/cannon/stratcert"],
         require=[("PRFull/check/cert/witnesses/*.json", 26, "the 25 witnesses and their index"),
                  ("PRFull/check/cert/witnesses/margins/*.tsv.gz", 25, "the exact margins of the 25 witnesses"),
                  ("PRFull/check/cert/results/*_box.json", 5, "the certified boxes"),
                  ("PRFull/check/ubound/results/ub_certificate.tsv.gz", 1, "the exclusion certificate"),
                  ("PRFull/check/wsep/out/*.tsv", 10, "the separation certificate"),
                  ("PRFull/cannon/bingap/out/task_*.txt", 64, "the census of binary pairs"),
                  ("PRFull/cannon/stratcert/out/t_*.tsv", 64, "the certified points of the stratum counts: the 64 task "
                   "outputs (arrays 48975432, 48981870)"),
                  ("PRFull/cannon/stratcert/out/t_*.done", 64, "the completion markers of the 64 tasks"),
                  ("PRFull/cannon/stratcert/out/stratcert_all.tsv", 1, "the 12 748 certified points, one line per unit"),
                  ("PRFull/cannon/stratcert/collect.out", 1, "the output of the checker collect.py (verdict True)"),
                  ("PRFull/cannon/stratcert/units.tsv", 1, "the retained strata (units) of the stratum counts"),
                  ("PRFull/cannon/stratcert/points.tsv", 1, "the census points the certified points start from")]),
    Item("accounting",
         "the Slurm accounting of every array (Computational Supplement, Table 2; SI Table 8 before it, and SI Table 7 "
         "before the second referee round)",
         trees=["PRFull/release/accounting"],
         globs=["PartnersRivalsM1f/cannon/JOBID*", "PartnersRivalsM1f/cannon/cw/JOBID*", "PartnersRivalsM2f/cannon/JOBID*",
                "PartnersRivalsM2f/cannon/*/JOBID*", "PRFull/cannon/*/JOBID*"],
         require=[("PRFull/release/accounting/sacct_table7.psv", 1, "sacct of every task and step"),
                  ("PRFull/release/accounting/jobs_tsv_lines.tsv", 1, "the submission ledger lines"),
                  ("PRFull/cannon/stratcert/JOBID_*", 2, "the JOBID files of the certified points (probe, full run)"),
                  ("PRFull/cannon/evoaudit/JOBID_*", 2, "the JOBID files of the audit (probe, full run)")]),
    Item("figures",
         "the figure scripts (with the figures, their statistics and SI Table 4 as deposited)",
         trees=["PRFull/Figures"],
         files=["PRFull/Text/si_table4_body.tex"],
         require=[("PRFull/Figures/fig*.py", 5, "the figure scripts"), ("PRFull/Figures/*.pdf", 8, "the figures")]),
    Item("paper",
         "the paper's sources (not promised; included so that reproduce_all.sh can rebuild it)",
         files=["PRFull/Text/ms.tex", "PRFull/Text/refs.tex", "PRFull/Text/reorder_bib.py", "PRFull/Text/ms.pdf",
                "PRFull/Text/si_table4.tex", "PRFull/Figures/sitable4_legend.tex", "PRFull/Text/si_box3_strata.tex",
                "PRFull/Text/si_box3_legend.tex", "PRFull/check/ref3/eq3_boundary/box/si_box3_strata.tex"]),
    Item("supplement",
         "the Computational Supplement (cs.pdf), which the statement and the SI say is included in the release, with "
         "its source and the script that fills its reference list from refs.tex (so that reproduce_all.sh can rebuild it)",
         files=["PRFull/Text/cs.tex", "PRFull/Text/cs_bib.py", "PRFull/Text/cs.pdf"]),
    Item("release",
         "the reproduction scripts, the SHA-256 manifest and this assembly",
         trees=["PRFull/release"],
         require=[("PRFull/release/reproduce_all.sh", 1, "reproduce_all.sh"),
                  ("PRFull/release/reproduce_certificates.sh", 1, "reproduce_certificates.sh"),
                  ("PRFull/release/reproduce_figures.sh", 1, "reproduce_figures.sh"),
                  ("PRFull/release/RELEASE.md", 1, "RELEASE.md"), ("PRFull/release/LICENSE", 1, "LICENSE"),
                  ("PRFull/release/figure_inputs.json", 1, "the traced inputs of the figure scripts")]),
    Item("programs",
         "the source of every program with its build script, compiler and flags, and the Slurm scripts and task files "
         "that fix the seed of every run (Computational Supplement, Table 1)",
         trees=["PartnersRivalsM1f/cannon", "PartnersRivalsM2f/cannon", "PRFull/cannon"],
         files=["PartnersRivalsM1f/README.md", "PartnersRivalsM2f/README.md", "PRFull/README.md"],
         require=[("PRFull/cannon/*/README.md", 6, "a README per Cannon kit of the revisions"),
                  ("PRFull/cannon/*/compile.sh", 6, "a build script per kit"),
                  ("PRFull/cannon/*/*.sbatch", 25, "the Slurm scripts of the revisions"),
                  ("PRFull/cannon/stratcert/*.sbatch", 2, "the Slurm scripts of the certified points (probe, full run)"),
                  ("PRFull/cannon/evoaudit/*.sbatch", 2, "the Slurm scripts of the audit (probe, full run)"),
                  ("PartnersRivalsM2f/cannon/*.sbatch", 3, "the Slurm scripts of the memory-two census"),
                  ("PartnersRivalsM1f/cannon/cw/*.slurm", 4, "the Slurm scripts of the memory-one evolution"),
                  ("PartnersRivalsM2f/cannon/wfc/tasks_*.txt", 3, "the task files of the memory-two evolution"),
                  ("PRFull/cannon/wfc2/tasks_*.txt", 10, "the task files of the sensitivity runs"),
                  ("PartnersRivalsM2f/cannon/census.c", 1, "the census program"),
                  ("PartnersRivalsM1f/cannon/cw/cw.f", 1, "the memory-one evolutionary program"),
                  ("PartnersRivalsM2f/cannon/wfc/wfc.c", 1, "the memory-two evolutionary program")]),
    Item("checks",
         "the programs, with their logs or collected output, of five checks that ran in the scratch of the memory-two "
         "derivation pass: the efficient sets by minimum arborescences (t2_fw.c), the closed forms of efficiency "
         "(effcc.c, t2_altform2.c), the strata of W of dimension 8 to 15, and the kit of the re-censuses (t2.c, t2b.c)",
         files=CHECKS),
    Item("inputs",
         "the other files the figure scripts, the kit analyses and the certificate programs read (traced, or named "
         "by an absolute path in an included program)",
         files=["PartnersRivals/Data/f1_e4.rep", "PartnersRivals/Data/f2_e4.rep"]),
]
PROJECT_TREES = [("PartnersRivalsM1f/check", "programs"), ("PartnersRivalsM2f/check", "programs")]


# ------------------------------------------------------------------------------------------------ helpers
def excluded(rel):
    if rel in CHECKS:                            # named one by one, in a directory that is otherwise left out
        return False
    pats = EXCLUDE_ANY if WITH_PAPER else EXCLUDE_ANY + PAPER_EXCLUDE
    return any(fnmatch.fnmatch(rel, p) or fnmatch.fnmatch("/" + rel, p) for p in pats)


def source_of(rel):
    """the file shipped at rel: its release-only copy in OVERRIDES if there is one, else the file itself"""
    o = os.path.join(OVERRIDES, rel)
    return o if os.path.isfile(o) else os.path.join(ROOT, rel)


def is_program_binary(path):
    try:
        with open(path, "rb") as fh:
            return fh.read(4) in MAGIC
    except OSError:
        return False


def walk(rel_top, skipped):
    top = os.path.join(ROOT, rel_top)
    if os.path.isfile(top):
        yield rel_top
        return
    for d, dirs, files in os.walk(top):
        reld = os.path.relpath(d, ROOT)
        keep = []
        for x in sorted(dirs):
            p = os.path.join(d, x)
            r = os.path.join(reld, x)
            if os.path.islink(p):
                skipped["symbolic link (not followed)"].append(r)
            elif excluded(r + "/x"):
                skipped["excluded directory"].append(r)
            elif any(fnmatch.fnmatch(r + "/x", p) for p in WALK_SKIP):
                continue                                         # taken by the items that name its files
            else:
                keep.append(x)
        dirs[:] = keep
        for f in sorted(files):
            r = os.path.join(reld, f)
            p = os.path.join(ROOT, r)
            if os.path.islink(p):
                skipped["symbolic link (not followed)"].append(r)
            elif excluded(r):
                continue
            elif is_program_binary(p):
                skipped["compiled code (the sources and build scripts are deposited)"].append(r)
            else:
                yield r


def glob_rel(pattern):
    import glob
    return sorted(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(ROOT, pattern)) if os.path.isfile(p))


def referenced_paths(rel):
    """absolute paths under the author's root named in a program file, as paths relative to the root"""
    try:
        s = open(os.path.join(ROOT, rel), encoding="utf-8", errors="replace").read()
    except OSError:
        return []
    out = []
    for m in re.finditer(re.escape(AUTHOR_ROOT) + r"/([A-Za-z0-9_./+-]+)", s):
        out.append(m.group(1).rstrip("./"))
    return out


def human(n):
    for u in ("B", "kB", "MB", "GB"):
        if n < 1000 or u == "GB":
            return ("%d %s" % (n, u)) if u == "B" else ("%.1f %s" % (n, u))
        n /= 1000.0


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ the selection
def select():
    chosen = collections.OrderedDict()          # rel -> item key (first item that claims it)
    skipped = collections.defaultdict(list)
    flags = []
    def add(rel, key):
        if rel not in chosen:
            chosen[rel] = key
    for it in ITEMS:
        if it.key in PAPER_ITEMS and not WITH_PAPER:
            continue
        for t in it.trees:
            if not os.path.exists(os.path.join(ROOT, t)):
                flags.append((it.key, "directory not found: %s" % t))
                continue
            for r in walk(t, skipped):
                add(r, it.key)
        for f in it.files:
            if os.path.isfile(os.path.join(ROOT, f)) and excluded(f):
                skipped["listed but excluded by rule"].append(f)
            elif os.path.isfile(os.path.join(ROOT, f)):
                add(f, it.key)
            else:
                flags.append((it.key, "file not found: %s" % f))
        for g in it.globs:
            for r in glob_rel(g):
                if not excluded(r) and not is_program_binary(os.path.join(ROOT, r)):
                    add(r, it.key)
    for t, key in PROJECT_TREES:
        for r in walk(t, skipped):
            add(r, key)
    # the traced inputs of the figure scripts
    fi = os.path.join(HERE, "figure_inputs.json")
    traced = 0
    if os.path.exists(fi):
        J = json.load(open(fi))
        for step, rec in sorted(J.get("steps", {}).items()):
            for r in rec.get("reads", []):
                if "/__pycache__/" in r:                       # a module read from its bytecode: take the source
                    r = re.sub(r"/__pycache__/([^/.]+)\.[^/]*\.pyc$", r"/\1.py", r)
                if r.startswith("PRFull/release/work/") or excluded(r):
                    continue
                if os.path.isfile(os.path.join(ROOT, r)):
                    if r not in chosen:
                        traced += 1
                    add(r, "inputs")
                else:
                    flags.append(("inputs", "traced input not found: %s (step %s)" % (r, step)))
    else:
        flags.append(("inputs", "figure_inputs.json not found: run reproduce_figures.sh --record"))
    # files named by absolute paths in the included programs
    byref, refdirs = 0, set()
    for r in list(chosen):
        if not r.endswith(PROGRAM_EXT):
            continue
        for q in referenced_paths(r):
            p = os.path.join(ROOT, q)
            if os.path.isfile(p):
                if not excluded(q) and not is_program_binary(p) and q not in chosen:
                    add(q, "inputs")
                    byref += 1
            elif os.path.isdir(p):
                refdirs.add(q)
    # promised items
    for it in ITEMS:
        if it.key in PAPER_ITEMS and not WITH_PAPER:
            continue
        for g, nmin, what in it.require:
            n = len([r for r in glob_rel(g) if r in chosen])
            if n < nmin and g.startswith("PRFull/release/harvest/"):
                continue                                         # said below, with where the files are
            if n < nmin:
                flags.append((it.key, "PROMISED BUT NOT FOUND: %s (%s: %d of at least %d)" % (what, g, n, nmin)))
        for what, where, size in it.missing:
            g = {"pt_": "pt_*.txt", "sw_": "sw_*", "dim_": "dim_*", "probe_": "probe_*"}
            known = [k for k in g if ("/" + k) in where]
            if known and any(r.startswith("PRFull/release/harvest/ca1_out/" + known[0]) for r in chosen):
                continue
            flags.append((it.key, "NOT IN THE TREE: %s; where: %s; size %s" % (what, where, size)))
    # every array whose JOBID file is in the release must be in the Slurm accounting (CS Table 2)
    acc = os.path.join(HERE, "accounting", "table7_jobs.tsv")
    ids = {l.split("\t")[0] for l in open(acc) if l.strip() and not l.startswith("#")} if os.path.exists(acc) else set()
    for r in chosen:
        if re.search(r"/([A-Z0-9]+_)?JOBID[^/]*$", r):              # JOBID, JOBID_PROBE, ... and WC_JOBID, ...
            for jid in re.findall(r"\b(\d{8})\b", open(os.path.join(ROOT, r)).read()):
                if jid not in ids:
                    flags.append(("accounting", "array %s (%s) is not in accounting/table7_jobs.tsv, so not in the deposited "
                                  "sacct dump or CS Table 2" % (jid, r)))
    # the manuscript guard: no manuscript and no copy of manuscript text (other than the paper itself with --with-paper)
    sys.path.insert(0, os.path.join(HERE, "tools"))
    import manuscript_guard
    # the paper's files, also where another item claims them first (the figures item walks PRFull/Figures, the
    # certificates item PRFull/check)
    paper_files = {f for it in ITEMS if it.key in PAPER_ITEMS for f in it.files}
    paper = {r for r in chosen if WITH_PAPER and (chosen[r] in PAPER_ITEMS or r in paper_files)}
    shipped = {os.path.relpath(source_of(r), ROOT): r for r in chosen if r not in paper}
    hits, refs = manuscript_guard.scan(ROOT, list(shipped), min_run=10)
    hits = [(shipped[h[0]],) + tuple(h[1:]) for h in hits]
    guard = {"refs": refs, "hits": hits}
    for r, kind, d, info in hits:
        if kind == "text":
            n = sum(v["words"] for v in info.values())
            comp = info.get("companion", {}).get("longest", 0)
            if n < GUARD_WORDS and comp < COMPANION_RUN:
                continue
            if r in GUARD_REVIEWED and n <= GUARD_REVIEWED[r][0]:
                continue
        flags.append(("guard", "MANUSCRIPT TEXT: %s (%s: %s)" % (r, kind, d)))
    if not refs:   # an error in the author's tree; in an assembled release, which never holds the sources (MS_EXCLUDE), a note
        flags.append(("guard" if os.path.realpath(ROOT) == AUTHOR_ROOT else "guard-off",
                      "no manuscript found to compare with: the guard checked names and PDF titles only"))
    return chosen, skipped, flags, traced, byref, sorted(refdirs), guard


def dest_path(rel):
    """where a file goes in the release tree (the harvested Cannon output goes beside the census programs)"""
    if rel.startswith("PRFull/release/harvest/ca1_out/"):
        return "PartnersRivalsM2f/cannon/ca1_out/" + rel[len("PRFull/release/harvest/ca1_out/"):]
    if rel == "PRFull/release/harvest/ca1_out.md5":
        return "PartnersRivalsM2f/cannon/ca1_out.md5"
    return rel


# ------------------------------------------------------------------------------------------------ harvest (Cannon)
def harvest():
    dst = os.path.join(HARVEST, "ca1_out")
    os.makedirs(dst, exist_ok=True)
    remote = "~/CL1/PartnersRivalsM2f/ca1/out"
    pats = "pt_*.txt sw_* dim_* probe_45749304*"
    print("md5sum on Cannon (read-only) ...", flush=True)
    md5 = subprocess.check_output(["ssh", "-o", "BatchMode=yes", "cannon", "cd %s && md5sum %s" % (remote, pats)], text=True)
    open(os.path.join(HARVEST, "ca1_out.md5"), "w").write(md5)
    print("rsync %d files from cannon:%s ..." % (len(md5.splitlines()), remote), flush=True)
    lst = "\n".join(l.split()[1] for l in md5.splitlines() if l.strip()) + "\n"
    p = subprocess.run(["rsync", "-a", "--files-from=-", "cannon:%s/" % remote, dst + "/"], input=lst, text=True)
    if p.returncode:
        sys.exit("rsync failed (%d)" % p.returncode)
    bad = 0
    for l in md5.splitlines():
        h, f = l.split()
        loc = hashlib.md5(open(os.path.join(dst, f), "rb").read()).hexdigest()
        bad += loc != h
    print("harvested %d files into %s; md5 mismatches: %d" % (len(md5.splitlines()), dst, bad))
    sys.exit(1 if bad else 0)


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--files", action="store_true", help="with --dry-run: list every file")
    ap.add_argument("--dest", default=None)
    ap.add_argument("--harvest-cannon", action="store_true")
    ap.add_argument("--items", default=None, help="with --dest: copy only these items (comma-separated keys), to test")
    ap.add_argument("--with-paper", action="store_true", help="include the paper and the Computational Supplement "
                    "(release 1.1.0, when the article is published); without it they are left out (release 1.0.0)")
    a = ap.parse_args()
    global WITH_PAPER
    WITH_PAPER = a.with_paper
    if a.harvest_cannon:
        return harvest()
    if not a.dry_run and not a.dest:
        ap.error("give --dry-run or --dest DIR")
    chosen, skipped, flags, traced, byref, refdirs, guard = select()
    size = {r: os.path.getsize(source_of(r)) for r in chosen}
    overridden = [r for r in chosen if source_of(r) != os.path.join(ROOT, r)]
    total = sum(size.values())
    by_item = collections.OrderedDict((it.key, [0, 0]) for it in ITEMS)
    for r, k in chosen.items():
        by_item[k][0] += 1
        by_item[k][1] += size[r]
    groups = collections.defaultdict(lambda: [0, 0])
    for r in chosen:
        parts = r.split("/")
        g = "/".join(parts[:3]) if len(parts) > 3 else "/".join(parts[:-1]) or r
        groups[g][0] += 1
        groups[g][1] += size[r]
    print("assemble_release.py: the release of %s (root %s), %s" % ("PRFull", ROOT, "WITH the paper (--with-paper)"
          if WITH_PAPER else "without the paper and the Computational Supplement (release 1.0.0)"))
    print("\nper item (a file counts once, under the first item that claims it):")
    for it in ITEMS:
        n, b = by_item[it.key]
        print("  %-13s %7d files %10s   %s" % (it.key, n, human(b), it.promise[:110]))
    print("\nper directory:")
    for g in sorted(groups):
        print("  %-64s %7d files %10s" % (g, groups[g][0], human(groups[g][1])))
    print("\ntotal: %d files, %s (%d bytes); %d added from the traced figure inputs, %d named by absolute paths in programs"
          % (len(chosen), human(total), total, traced, byref))
    if refdirs:
        print("directories named by absolute paths in programs (their files are included only through the rules above):")
        for d in refdirs:
            print("    " + d)
    print("\nshipped as release-only copies (PRFull/release/overrides): %d" % len(overridden))
    for r in overridden:
        print("    %s" % r)
    print("\nleft out:")
    for why, lst in sorted(skipped.items()):
        b = sum(os.path.getsize(os.path.join(ROOT, r)) for r in lst if os.path.isfile(os.path.join(ROOT, r)) and not os.path.islink(os.path.join(ROOT, r)))
        print("  %-62s %6d %s  e.g. %s" % (why, len(lst), human(b) if b else "", ", ".join(lst[:3])))
    print("  by rule: " + ", ".join(p for p in EXCLUDE_ANY if not p.startswith("*")))
    print("\nmanuscript guard: compared with %d manuscript and draft files%s" % (len(guard["refs"]),
          "" if guard["refs"] else " (tools/manuscript_sources.json not on this machine: names only)"))
    def _n(h):
        return sum(v["words"] for v in h[3].values())
    def _flagged(h):
        return _n(h) >= GUARD_WORDS or h[3].get("companion", {}).get("longest", 0) >= COMPANION_RUN
    small = [h for h in guard["hits"] if h[1] == "text" and not _flagged(h)]
    print("  %d selected files share short runs (< %d words in runs of >= 10, and no run of %d from the companion) with a "
          "manuscript, mostly definitions and formulas quoted in comments; listed for review, not flagged:"
          % (len(small), GUARD_WORDS, COMPANION_RUN))
    for r, _, d, _ in sorted(small, key=lambda h: -_n(h)):
        print("    %s: %s" % (r, d))
    rev = [h for h in guard["hits"] if h[1] == "text" and _flagged(h) and h[0] in GUARD_REVIEWED and _n(h) <= GUARD_REVIEWED[h[0]][0]]
    print("  %d files above that, read and accepted (GUARD_REVIEWED):" % len(rev))
    for r, _, d, _ in rev:
        print("    %s: %s; %s" % (r, d, GUARD_REVIEWED[r][1]))
    print("\nflags (%d):" % len(flags))
    for k, f in flags:
        print("  [%s] %s" % (k, f))
    if a.dry_run:
        if a.files:
            print("\nfiles:")
            for r in chosen:
                print("%12d  %-13s %s" % (size[r], chosen[r], dest_path(r)))
        print("\ndry run: nothing copied")
        return 0
    # ---------------------------------------------------------------------------------------- copy and manifest
    nguard = sum(1 for k, _ in flags if k == "guard")
    if nguard:
        sys.exit("\nrefused: the manuscript guard has %d flag(s), listed above; fix them, or accept them after review in "
                 "GUARD_REVIEWED, before assembling. Nothing was copied." % nguard)
    dest = os.path.abspath(a.dest)
    if os.path.exists(dest):
        sys.exit("%s exists; give a new directory" % dest)
    if dest.startswith(ROOT + os.sep + "PRFull"):
        sys.exit("the release tree must not lie inside PRFull")
    os.makedirs(dest)
    rows = []
    if a.items:
        keep = set(a.items.split(","))
        chosen = collections.OrderedDict((r, k) for r, k in chosen.items() if k in keep)
        total = sum(size[r] for r in chosen)
        print("\n--items %s: a PARTIAL tree of %d files, for testing only" % (a.items, len(chosen)))
    for r in chosen:
        d = os.path.join(dest, dest_path(r))
        os.makedirs(os.path.dirname(d), exist_ok=True)
        shutil.copy2(source_of(r), d)
        rows.append((dest_path(r), size[r], chosen[r]))
    for src, name in (("RELEASE.md", "README.md"), ("LICENSE", "LICENSE")):
        shutil.copy2(os.path.join(HERE, src), os.path.join(dest, name))
        rows.append((name, os.path.getsize(os.path.join(dest, name)), "release"))
    with open(os.path.join(dest, "reproduce_all.sh"), "w") as fh:
        fh.write('#!/usr/bin/env bash\n# the release root: run PRFull/release/reproduce_all.sh (see README.md)\n'
                 'exec bash "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/PRFull/release/reproduce_all.sh" "$@"\n')
    os.chmod(os.path.join(dest, "reproduce_all.sh"), 0o755)
    rows.append(("reproduce_all.sh", os.path.getsize(os.path.join(dest, "reproduce_all.sh")), "release"))
    with open(os.path.join(dest, "CONTENTS.tsv"), "w") as fh:
        fh.write("path\tbytes\titem\n")
        for p, n, k in sorted(rows):
            fh.write("%s\t%d\t%s\n" % (p, n, k))
    rows.append(("CONTENTS.tsv", os.path.getsize(os.path.join(dest, "CONTENTS.tsv")), "release"))
    with open(os.path.join(dest, "MANIFEST.sha256"), "w") as fh:
        for p, _, _ in sorted(rows):
            fh.write("%s  %s\n" % (sha256(os.path.join(dest, p)), p))
    print("\ncopied %d files (%s) into %s; wrote MANIFEST.sha256 and CONTENTS.tsv" % (len(rows), human(total), dest))
    print("check it with:  cd %s && shasum -a 256 -c MANIFEST.sha256 | grep -v ': OK$'" % dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
