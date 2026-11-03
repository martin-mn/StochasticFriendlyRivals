#!/usr/bin/env python3
"""hosting_plan.py -- split an assembled release between the GitHub repository and the Zenodo record.

    python3 PRFull/release/hosting/hosting_plan.py TREE [OUTDIR]

TREE is a release tree written by assemble_release.py --dest (it reads TREE/CONTENTS.tsv).  Writes to OUTDIR (default:
this directory):
    github_files.txt            the files of the GitHub repository (every file of CONTENTS.tsv that no rule below sends
                                to Zenodo only; GitHub refuses files over 100 MB, and the repository is kept small)
    zenodo_archive_of_file.tsv  path, bytes, archive: each file of CONTENTS.tsv in exactly one archive of the Zenodo
                                record (the code archive is the GitHub tree)
    hosting_plan.txt            the sizes and counts of both, the largest files, and which archives the quick
                                verification (reproduce_all.sh --quick), and the full run, need beyond the GitHub
                                repository

Rules: the first that matches a path decides.  The INCLUDE rules keep on GitHub, inside directories that otherwise go
to Zenodo, the 80 run-length series that SI Figure 3 reads and the raw output of the dimension probe (dim_*, 15 MB)
that SI Table 1 (Figures/sitable1.py) reads.  pack_release.sh packs the archives; fetch_zenodo.sh fetches and unpacks
them into a clone of the repository.
"""
import collections
import os
import re
import sys

PREFIX = "PRFull-1.0.0"
INCLUDE = [r"PRFull/cannon/cw2/out/TS/p(029|171|415|500)_r\d\d/t[^/]*$",
           r"PRFull/cannon/wfc2/out/TS/r_1o3_(30|172|416|501)_[^/]*\.txt$",
           r"PartnersRivalsM2f/cannon/ca1_out/dim_"]
# (regex, archive) -- every path matched goes to Zenodo only, in that archive
ZENODO = [
    (r"DiskM2WF/Data/[^/]*\.bin$", "binary-paper-dumps"),
    (r"PartnersRivals1/Data/[^/]*\.bin$", "binary-paper-dumps"),
    (r"PartnersRivals1?/Cannon/[^/]*\.bin$", "binary-paper-dumps"),
    (r"PartnersRivalsM2f/cannon/archive/", "m2-archives-and-audit"),
    (r"PRFull/cannon/evoaudit/(out|out_probe|logs)/", "m2-archives-and-audit"),
    (r"PartnersRivalsM2f/cannon/ca1_out/", "census"),
    (r"PartnersRivalsM2f/cannon/ca1_out\.md5$", "census"),
    (r"PRFull/release/census_class/[^/]*\.tgz$", "census"),
    (r"PRFull/release/census_class/[^/]*support_class\.txt\.gz$", "census"),
    (r"PRFull/cannon/wfc2/out/", "evolution-runs"),
    (r"PRFull/cannon/cw2/out/", "evolution-runs"),
    (r"PartnersRivalsM1f/cannon/cw/(run1|run2|run_local|run_local2|vcheck|vcheck2|kappa|kappa_local|probe|probe2|"
     r"packed_local|packed_local2|quarantine|log)/", "evolution-runs"),
    (r"PartnersRivalsM1f/cannon/out/", "evolution-runs"),
    (r"PRFull/check/wsep/out/(cover|extra)/", "certificates"),
    (r"PRFull/check/wsep/out/(c6|c8|c68|c9b)\.exact\.tsv$", "certificates"),
    (r"PRFull/check/wsep/out/corner_", "certificates"),
    (r"PRFull/cannon/stratcert/prep/census_pts_PA\.txt$", "certificates"),
    (r"(^|/)logs/", "certificates"),
    (r"PartnersRivalsM1f/check/working/census\.pkl$", "certificates"),
]
# what each quick step needs beyond the GitHub repository (reproduce_figures.sh, reproduce_certificates.sh), and what
# the full scripts need beyond that
QUICK_NEEDS = {"fig45": "binary-paper-dumps", "sitable4": "binary-paper-dumps", "evoaudit": "m2-archives-and-audit"}
FULL_NEEDS = {"table5": "evolution-runs", "separation": "certificates"}


def place(path):
    if any(re.match(rx, path) for rx in INCLUDE):
        return "code"
    for rx, arc in ZENODO:
        if re.match(rx, path) or (rx.startswith("(^|/)") and re.search(rx, path)):
            return arc
    return "code"


def mb(n):
    return "%.1f MB" % (n / 1e6)


def main():
    tree = os.path.abspath(sys.argv[1])
    out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(__file__))
    os.makedirs(out, exist_ok=True)
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(tree, "CONTENTS.tsv"))][1:]
    rows = [(p, int(n)) for p, n, *_ in rows]
    arc = {p: place(p) for p, _ in rows}
    per = collections.OrderedDict((a, [0, 0, ("", 0)]) for a in
                                  ["code", "binary-paper-dumps", "m2-archives-and-audit", "census", "evolution-runs", "certificates"])
    for p, n in rows:
        r = per[arc[p]]
        r[0] += 1; r[1] += n
        if n > r[2][1]:
            r[2] = (p, n)
    with open(os.path.join(out, "github_files.txt"), "w") as fh:
        for p, _ in sorted(rows):
            if arc[p] == "code":
                fh.write(p + "\n")
    with open(os.path.join(out, "zenodo_archive_of_file.tsv"), "w") as fh:
        fh.write("path\tbytes\tarchive\n")
        for p, n in sorted(rows):
            fh.write("%s\t%d\t%s_%s.tar.gz\n" % (p, n, PREFIX, arc[p]))
    gh = [(p, n) for p, n in rows if arc[p] == "code"]
    L = ["# PRFull 1.0.0 hosting plan, from the CONTENTS.tsv of the assembled tree: %d files, %d bytes (%s); "
         "MANIFEST.sha256 and CONTENTS.tsv come on top" % (len(rows), sum(n for _, n in rows), mb(sum(n for _, n in rows))), "",
         "## GitHub repository (the code archive): %d files, %d bytes (%s); largest %s (%s); files over 50 MB: %d, over "
         "100 MB: %d" % (len(gh), sum(n for _, n in gh), mb(sum(n for _, n in gh)), per["code"][2][0], mb(per["code"][2][1]),
                         sum(n > 50e6 for _, n in gh), sum(n > 100e6 for _, n in gh)), "",
         "## Zenodo record: README.md, LICENSE, MANIFEST.sha256, CONTENTS.tsv, and the archives (each prefixed %s/, so that "
         "unpacking all of them in one directory rebuilds the tree):" % PREFIX]
    for a, (k, n, big) in per.items():
        L.append("  %-45s %6d files %14d bytes %10s   largest %s (%s)" % ("%s_%s.tar.gz" % (PREFIX, a), k, n, mb(n), big[0], mb(big[1])))
    L.append("  sum: %d files, %d bytes (CONTENTS.tsv: %d files, %d bytes)" % (sum(v[0] for v in per.values()),
             sum(v[1] for v in per.values()), len(rows), sum(n for _, n in rows)))
    L += ["", "## the quick verification from the GitHub repository alone",
          "  every step runs except: " + ", ".join("%s (needs %s_%s.tar.gz)" % (s, PREFIX, a) for s, a in QUICK_NEEDS.items()),
          "  those steps skip and say so; bash PRFull/release/fetch_zenodo.sh fetches the two archives, unpacks them into the "
          "tree and checks every file against MANIFEST.sha256.", "",
          "## the full run from the GitHub repository alone",
          "  the full scripts also skip: " + ", ".join("%s (needs %s_%s.tar.gz)" % (s, PREFIX, a) for s, a in FULL_NEEDS.items()),
          "  bash PRFull/release/fetch_zenodo.sh --all fetches all five data archives.", "",
          "## the 30 largest files"]
    for p, n in sorted(rows, key=lambda r: -r[1])[:30]:
        L.append("  %12d %9s  %-22s %s" % (n, mb(n), arc[p], p))
    top = collections.defaultdict(lambda: [0, 0, 0, 0])
    for p, n in rows:
        t = top[p.split("/")[0] if "/" in p else "(top level)"]
        t[0] += 1; t[1] += n
        if arc[p] == "code":
            t[2] += 1; t[3] += n
    L += ["", "## per top-level directory: files, size | on GitHub: files, size"]
    for d in sorted(top):
        t = top[d]
        L.append("  %-20s %6d %10s | %6d %10s" % (d, t[0], mb(t[1]), t[2], mb(t[3])))
    open(os.path.join(out, "hosting_plan.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:14]))
    print("wrote github_files.txt, zenodo_archive_of_file.tsv and hosting_plan.txt in " + out)


if __name__ == "__main__":
    main()
