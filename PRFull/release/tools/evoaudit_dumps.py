#!/usr/bin/env python3
"""evoaudit_dumps.py -- the memory-two strategy dumps that cannon/evoaudit samples, taken from the deposited archives,
and the checks of the evoaudit step of reproduce_certificates.sh.

    python3 evoaudit_dumps.py extract ARCHIVE_DIR KIT     write KIT/data/F1/S_1o3_<task>_1.bin and KIT/data/F2/... from
                                                          ARCHIVE_DIR/cannon_2026-09-26.tgz (PartnersRivalsM2f_wfc/out/,
                                                          N = 1000) and ARCHIVE_DIR/cannon_out_f2_2026-09-26.tgz
                                                          (out_f2/, N = 100); check that each regime has its 512 games
                                                          and that every sampled strategy of KIT/sample/sample.npy is
                                                          its record in the dump (irep, count, atom and the 16 doubles,
                                                          exactly)
    python3 evoaudit_dumps.py summary RESULTS             the audit's result from RESULTS/summary.json (analyze.py):
                                                          audited, errors, missing, unresolved, changed atoms and bits

The release does not duplicate the dumps (0.67 GB): the kit's data/ is a local extraction of the same files, and the
archives are deposited in PartnersRivalsM2f/cannon/archive.  Only the Python standard library and numpy are needed.
"""
import collections
import json
import os
import re
import sys
import tarfile

import numpy as np

REC = np.dtype([("irep", "<i4"), ("count", "<i4"), ("atom", "<i4"), ("p", "<f8", (16,))])      # = src/load.py
SOURCES = [("F1", "cannon_2026-09-26.tgz", r"^(\./)?PartnersRivalsM2f_wfc/out/S_1o3_(\d+)_1\.bin$"),
           ("F2", "cannon_out_f2_2026-09-26.tgz", r"^(\./)?out_f2/S_1o3_(\d+)_1\.bin$")]
NGAMES, SPR = 512, 40000


def extract(arch, kit):
    nrec = {}
    for reg, name, rx in SOURCES:
        dst = os.path.join(kit, "data", reg)
        os.makedirs(dst, exist_ok=True)
        got = set()
        with tarfile.open(os.path.join(arch, name), "r|gz") as tf:      # one pass over the stream
            for m in tf:
                mm = re.match(rx, m.name)
                if not mm or not m.isfile():
                    continue
                task = int(mm.group(2))
                with open(os.path.join(dst, "S_1o3_%d_1.bin" % task), "wb") as fh:
                    fh.write(tf.extractfile(m).read())
                got.add(task)
        if got != set(range(1, NGAMES + 1)):
            sys.exit("FAIL: %s: %d of the %d dumps in %s" % (reg, len(got), NGAMES, name))
        n = 0
        for t in range(1, NGAMES + 1):
            d = np.fromfile(os.path.join(dst, "S_1o3_%d_1.bin" % t), dtype=REC)
            if d["count"].sum() != SPR:
                sys.exit("FAIL: %s task %d: %d samples, not %d" % (reg, t, d["count"].sum(), SPR))
            n += len(d)
        nrec[reg] = n
    print("dumps: %d files, %d records (F1 %d, F2 %d) from %s" % (2 * NGAMES, sum(nrec.values()), nrec["F1"], nrec["F2"],
                                                                 " and ".join(s[1] for s in SOURCES)))
    S = np.load(os.path.join(kit, "sample", "sample.npy"))
    ok = 0
    for reg in (1, 2):
        for task in np.unique(S["task"][S["reg"] == reg]):
            m = np.where((S["reg"] == reg) & (S["task"] == task))[0]
            d = np.fromfile(os.path.join(kit, "data", "F%d" % reg, "S_1o3_%d_1.bin" % task), dtype=REC)[S["row"][m]]
            same = ((d["irep"] == S["irep"][m]) & (d["count"] == S["count"][m]) & (d["atom"] == S["atom"][m])
                    & (d["p"] == S["p"][m]).all(axis=1))
            ok += int(same.sum())
    print("sampled strategies equal their records in the dumps: %d of %d" % (ok, len(S)))
    return 0 if ok == len(S) else 1


def summary(res):
    J = json.load(open(os.path.join(res, "summary.json")))
    unres = {k: g["props"]["unresolved (not audited or error)"]["n_y"] for k, g in J["groups"].items()}
    bits = {b: sum(J["groups"][r]["props"]["changed " + b]["n_y"] for r in ("F1", "F2")) for b in "ECS"}
    flips = collections.Counter("%s->%s" % (c["prod"], c["val"]) for c in J["changed_records"])
    print("audited %d, errors %d, missing %d, unresolved %d; changed atoms %d (%s); changed bits: E %d, C %d, S %d; "
          "validated 101: %d; production atom not reproduced: %d"
          % (J["audited"], J["errors"], J["missing"], unres["F1"] + unres["F2"], len(J["changed_records"]),
             ", ".join("%s %d" % kv for kv in sorted(flips.items())) or "none", bits["E"], bits["C"], bits["S"],
             J["validated_101"], len(J["not_reproduced"])))
    for c in J["changed_records"]:
        print("  changed: sid %d F%d %s task %d: %s -> %s (bit %s)" % (c["sid"], c["reg"], c["wedge"], c["task"], c["prod"],
                                                                       c["val"], c["changed"]))
    return 0


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "extract":
        sys.exit(extract(sys.argv[2], sys.argv[3]))
    if len(sys.argv) == 3 and sys.argv[1] == "summary":
        sys.exit(summary(sys.argv[2]))
    sys.exit(__doc__)
