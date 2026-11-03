#!/usr/bin/env python3
"""wsep_same.py -- is a regenerated W-separation certificate (check/wsep) the deposited one, as a certificate?

    python3 wsep_same.py NEWDIR DEPDIR CLASS [CLASS ...]      e.g.  out ../../../check/wsep/out c9a extra/x4

wsep excludes each stratum with the first co-player that works, trying a move-to-front pool before the full search;
with several threads the pool's order, and hence the co-player named in a row, the split between "pool" and "full"
and the list of co-players used, depend on the thread schedule.  So a rerun is not byte-identical.  What must agree
is the certificate's content: for each CLASS, the summary's counts (total, efficient, excluded by ALLD, SI ALLD
disagreements, survivors, NOT_EXCLUDED, lemma violations, the per-dimension counts, and pool + full), and the set of
strata in CLASS.tsv (every row a stratum excluded by some co-player, "how" P or F).  A class other than the negative
control must have NOT_EXCLUDED 0 and lemma_violations 0.  Exit status 1 if anything differs or fails.
"""
import os
import sys

KEYS = ("class", "efficient", "excluded_by_ALLD", "SI_alld_disagreements", "survivors", "NOT_EXCLUDED", "lemma_violations",
        "by_d")


def summary(path):
    d = {}
    for l in open(path):
        k, _, v = l.strip().partition(" ")
        if k in KEYS or k in ("excluded_pool", "excluded_full"):
            d[k] = v
    d["excluded"] = int(d.get("excluded_pool", 0)) + int(d.get("excluded_full", 0))
    d.pop("excluded_pool", None)
    d.pop("excluded_full", None)
    return d


def rows(path):
    pats, how = [], set()
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) > 5:
            pats.append(f[1])
            how.add(f[5])
    return pats, how


def main():
    new, dep, classes = sys.argv[1], sys.argv[2], sys.argv[3:]
    bad = 0
    for c in classes:
        sn, sd = summary(os.path.join(new, c + ".summary")), summary(os.path.join(dep, c + ".summary"))
        pn, hn = rows(os.path.join(new, c + ".tsv"))
        pd, _ = rows(os.path.join(dep, c + ".tsv"))
        msg = []
        if sn != sd:
            msg.append("summary differs: " + ", ".join("%s %s/%s" % (k, sn.get(k), sd.get(k)) for k in sorted(set(sn) | set(sd))
                                                       if sn.get(k) != sd.get(k)))
        if sorted(pn) != sorted(pd):
            msg.append("strata differ: %d rows here, %d deposited, %d in common" % (len(pn), len(pd), len(set(pn) & set(pd))))
        if not hn <= {"P", "F"}:
            msg.append("unknown 'how' %s" % sorted(hn - {"P", "F"}))
        if sn.get("NOT_EXCLUDED") != "0" or sn.get("lemma_violations") != "0":
            msg.append("NOT_EXCLUDED %s, lemma_violations %s" % (sn.get("NOT_EXCLUDED"), sn.get("lemma_violations")))
        bad += bool(msg)
        print("%-12s %7d strata, %7d excluded: %s" % (c, len(pn), sn["excluded"],
              "; ".join(msg) if msg else "same strata and verdicts as the deposit (co-players may differ)"))
    print("wsep classes compared: %d, differing: %d" % (len(classes), bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
