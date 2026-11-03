#!/usr/bin/env python3
"""check_table7.py -- recompute CS Table 2 (the cost of the computations, in the Computational Supplement cs.tex)
from the deposited Slurm accounting.  (CS Table 2 was SI Table 8, and SI Table 7 before the second referee round;
the file names keep that old number.)

    python3 check_table7.py            reads sacct_table7.psv and table7_jobs.tsv here; prints one line per row of
                                       CS Table 2 and a final verdict; exit status 1 if a row differs

Per array: the tasks (lines of sacct whose JobID has no step suffix), their states, the wall-clock time per task
(ElapsedRaw, shortest to longest), the span from the first start to the last end, the core-hours (CPUTimeRAW summed
over the tasks, as the Computational Supplement, Section 6, says), TotalCPU summed, and the peak resident set (MaxRSS, largest over
the job steps).  Per row of CS Table 2 the core-hours of its arrays are summed and compared with the table, to the
table's rounding (0.05 core-hours, 0.5 for the total).  Needs only the Python standard library.
"""
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PSV = os.path.join(HERE, "sacct_table7.psv")
JOBS = os.path.join(HERE, "table7_jobs.tsv")
TABLE = "CS Table 2"
TABLE_TOTAL = 8569                                   # the total of CS Table 2, core-hours


def secs(t):
    """sacct TotalCPU: [DD-][HH:]MM:SS[.mmm]"""
    if not t:
        return 0.0
    d = 0
    if "-" in t:
        d, t = t.split("-")
        d = int(d)
    parts = [float(x) for x in t.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    return d * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def mib(r):
    if not r:
        return 0.0
    m = re.match(r"([\d.]+)([KMGT]?)", r)
    v, u = float(m.group(1)), m.group(2)
    return v * {"": 1 / 1048576, "K": 1 / 1024, "M": 1, "G": 1024, "T": 1048576}[u]


def hm(s):
    s = int(round(s / 60.0))
    return "%d:%02d" % (s // 60, s % 60)


rows = [l.rstrip("\n").split("\t") for l in open(JOBS) if l.strip() and not l.startswith("#")]
lines = [l.rstrip("\n").split("|") for l in open(PSV)]
hdr, lines = lines[0], lines[1:]
ix = {k: i for i, k in enumerate(hdr)}
arr = collections.defaultdict(lambda: dict(tasks=0, states=collections.Counter(), el=[], cpu=0, tot=0.0,
                                            start=None, end=None, rss=0.0))
for L in lines:
    jid = L[ix["JobID"]]
    base = jid.split("_")[0].split(".")[0]
    a = arr[base]
    a["rss"] = max(a["rss"], mib(L[ix["MaxRSS"]]))
    if "." in jid:                                   # a job step: only its MaxRSS is used
        continue
    a["tasks"] += 1
    a["states"][L[ix["State"]].split()[0]] += 1
    a["el"].append(int(L[ix["ElapsedRaw"]] or 0))
    a["cpu"] += int(L[ix["CPUTimeRAW"]] or 0)
    a["tot"] += secs(L[ix["TotalCPU"]])
    for k, f in (("start", min), ("end", max)):
        v = L[ix[k.capitalize()]]
        if v and v not in ("Unknown", "None"):
            a[k] = v if a[k] is None else f(a[k], v)

bad = 0
tot_ch = tot_cpu = 0.0
by_row = collections.OrderedDict()
for jid, row, ch, kit in rows:
    by_row.setdefault(row, dict(jobs=[], table=None))
    by_row[row]["jobs"].append(jid)
    if ch != "-":
        by_row[row]["table"] = float(ch)
print("%-10s %6s %-26s %-13s %9s %9s %8s  %s" % ("job", "tasks", "states", "per task", "core-h", "TotalCPU", "MaxRSS", "row"))
for jid, row, ch, kit in rows:
    a = arr.get(jid)
    if a is None or a["tasks"] == 0:
        print("%-10s  NOT IN THE DUMP  (%s)" % (jid, row))
        bad += 1
        continue
    st = ",".join("%s %d" % kv for kv in sorted(a["states"].items()))
    print("%-10s %6d %-26s %5s-%-7s %9.1f %9.1f %6.0f MiB  %s" % (jid, a["tasks"], st[:26], hm(min(a["el"])), hm(max(a["el"])),
                                                              a["cpu"] / 3600, a["tot"] / 3600, a["rss"], row))
print()
for row, r in by_row.items():
    ch = sum(arr[j]["cpu"] for j in r["jobs"] if j in arr) / 3600.0
    tot_ch += ch
    tot_cpu += sum(arr[j]["tot"] for j in r["jobs"] if j in arr) / 3600.0
    ok = r["table"] is not None and abs(ch - r["table"]) <= 0.05 + 1e-9
    bad += (not ok)
    print("%-32s %2d arrays  core-h %8.2f  %s %8s  %s" % (row, len(r["jobs"]), ch, TABLE, r["table"], "ok" if ok else "DIFFERS"))
ok = abs(tot_ch - TABLE_TOTAL) <= 0.5
bad += (not ok)
print("%-32s %2d arrays  core-h %8.2f  %s %8s  %s" % ("total", len(rows), tot_ch, TABLE, TABLE_TOTAL, "ok" if ok else "DIFFERS"))
print("TotalCPU is %.2f%% below CPUTimeRAW in total (the Computational Supplement: 0.6%%); largest MaxRSS %.0f MiB"
      % (100 * (1 - tot_cpu / tot_ch), max(a["rss"] for a in arr.values())))
print("rows or arrays that differ:", bad)
sys.exit(1 if bad else 0)
