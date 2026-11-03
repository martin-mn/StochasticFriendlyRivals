#!/usr/bin/env python3
"""join_census.py -- the support-based efficiency class joined to EVERY line of the archived census output
(referee 2, additional comment 3), not only to the efficient lines of PartnersRivalsM2f/check/st_eff.txt.

    nice -n 10 python3 join_census.py [--tgz census_st_45751662.tgz] [--md5 census_st_45751662.md5]
                                                                         (about 1 min; writes only in this directory)

Input.
  census_st_45751662.tgz   the stratum output of the census, Cannon array 45751662 (PartnersRivalsM2f/cannon/census.c,
                           mode "strata", 512 tasks): st_<t>.txt, one line "g cls nP nM K" for every stratum in which at
                           least one of its K draws was confirmed competitive against all 65536 binary co-players
                           (cls read from the weights at the stratum's first draw, tolerance 1e-12), and the task logs
                           st_45751662_<t>.out/.err.  Harvested on 2026-09-27 from cannon:~/CL1/PartnersRivalsM2f/ca1/out
                           (tar over ssh); every member's md5 equals the md5sum taken on Cannon (census_st_45751662.md5).
  census_support_class.u8.gz  the class decided by the support alone for all 3^16 strata (mk_census_class.py).
Output.
  census_all_support_class.txt.gz   every census line, tasks in order 0..511, with a sixth column: the support-based
                                    class (-, C, A or Z) of that stratum
  census_class.json                 updated with the counts of this join (key "full_census_join")
Checks (exit status 1 if one fails): every md5 matches; 512 task files and 512 empty .err-free logs are present
(each .err holds the task's one-line summary; task 192's also 14 lines "srun: error: wrap_on_data ... Insane message length"
from Slurm, the task COMPLETED); 3 254 063 lines; no stratum twice; the census class and the support class
agree on every line (none of the 1 729 272 lines the census read as not efficient is efficient by its support, and none
of the 1 524 791 efficient lines differs; the latter are the lines of st_eff.txt).
"""
import argparse
import collections
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tarfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = {0: "-", 1: "C", 2: "A", 3: "Z"}
INV = {v: k for k, v in CODE.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tgz", default=os.path.join(HERE, "census_st_45751662.tgz"))
    ap.add_argument("--md5", default=os.path.join(HERE, "census_st_45751662.md5"))
    a = ap.parse_args()
    t0 = time.time()
    want = {}
    for l in open(a.md5):
        if l.strip():
            h, f = l.split()
            want[f] = h
    got, data = {}, {}
    with tarfile.open(a.tgz, "r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            b = tf.extractfile(m).read()
            got[m.name] = hashlib.md5(b).hexdigest()
            data[m.name] = b
    bad_md5 = sorted(f for f in want if got.get(f) != want[f]) + sorted(f for f in got if f not in want)
    tasks = sorted((int(re.match(r"st_(\d+)\.txt$", f).group(1)) for f in data if re.match(r"st_\d+\.txt$", f)))
    ok_tasks = tasks == list(range(512))
    errs = [f for f in data if f.endswith(".err")]
    srun = collections.Counter()                  # Slurm's own messages in a task's stderr (task 192: 14 lines of
    err_ok = True                                 # "srun: error: wrap_on_data ... Insane message length"; COMPLETED)
    for f in errs:
        L = [l for l in data[f].decode().split("\n") if l.strip()]
        st = [l for l in L if re.match(r"^stripe \d+/512 KEFF=64 KOTH=8: emitted \d+ strata, \d+ full 65536-scans$", l)]
        other = [l for l in L if l not in st]
        err_ok &= len(st) == 1 and all(l.startswith("srun: error: wrap_on_data") for l in other)
        if other:
            srun[f] = len(other)
    outs_empty = all(len(data[f]) == 0 for f in data if f.endswith(".out"))
    sup = np.frombuffer(gzip.open(os.path.join(HERE, "census_support_class.u8.gz")).read(), dtype=np.uint8)
    assert sup.size == 3 ** 16
    rows = []
    for t in range(512):
        for l in data["st_%d.txt" % t].decode().split("\n"):
            if l.strip():
                rows.append(l.split())
    g = np.array([int(r[0]) for r in rows], dtype=np.int64)
    cen = np.array([INV[r[1]] for r in rows], dtype=np.uint8)
    s = sup[g]
    dup = int(g.size - np.unique(g).size)
    tab = collections.Counter(zip((CODE[int(x)] for x in cen), (CODE[int(x)] for x in s)))
    mism = int((cen != s).sum())
    joined = "".join("%s %s\n" % (" ".join(r), CODE[int(c)]) for r, c in zip(rows, s)).encode()
    with gzip.GzipFile(os.path.join(HERE, "census_all_support_class.txt.gz"), "wb", mtime=0) as fh:
        fh.write(joined)
    by_cen = {k: int((cen == v).sum()) for k, v in INV.items()}
    rep_path = os.path.join(HERE, "census_class.json")
    rep = json.load(open(rep_path)) if os.path.exists(rep_path) else {}
    rep["full_census_join"] = dict(
        source="cannon:~/CL1/PartnersRivalsM2f/ca1/out/st_*.txt (array 45751662), harvested 2026-09-27, md5 checked",
        files=len(data), md5_mismatches=len(bad_md5), task_files_0_to_511=ok_tasks, err_logs_ok=err_ok,
        out_logs_empty=outs_empty, slurm_messages_in_err_logs=dict(srun), lines=len(rows), duplicate_strata=dup, lines_by_census_class=by_cen,
        census_class_vs_support_class={"%s/%s" % k: v for k, v in sorted(tab.items())},
        mismatches=mism,
        sha256_census_all_support_class_txt=hashlib.sha256(joined).hexdigest(),
        sha256_tgz=hashlib.sha256(open(a.tgz, "rb").read()).hexdigest())
    with open(rep_path, "w") as fh:
        json.dump(rep, fh, indent=1)
    print("census output: %d files, md5 mismatches %d, tasks 0..511 %s, logs %s%s" %
          (len(data), len(bad_md5), "ok" if ok_tasks else "MISSING", "ok" if (err_ok and outs_empty) else "UNEXPECTED",
           "".join("; %s holds %d srun messages besides its summary" % kv for kv in sorted(srun.items()))))
    print("census lines %d (by the census's class %s), strata listed twice %d" % (len(rows), by_cen, dup))
    print("census class / support class: %s" % ", ".join("%s/%s %d" % (k[0], k[1], v) for k, v in sorted(tab.items())))
    print("lines whose two classes differ: %d" % mism)
    print("wrote census_all_support_class.txt.gz and census_class.json in %.1f s" % (time.time() - t0))
    ok = not bad_md5 and ok_tasks and err_ok and outs_empty and dup == 0 and mism == 0 and len(rows) == 3254063
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
