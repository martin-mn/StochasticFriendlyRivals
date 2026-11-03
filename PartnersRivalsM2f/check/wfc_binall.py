#!/usr/bin/env python3
"""wfc_binall.py <dir-with-B<task>.bin>  --  the continuum classifier of cannon/wfc/wfc.c, fed every one of the 65536
binary genomes, against Revised2's exact atom table PartnersRivals/Cannon/atoms512.bin (512 x 65536 bytes, ipt-major,
code = 4 efficient + 2 Nash + competitive).  Task t is line t of games.dat, whose ipt is t-1 (gid 2000+ipt).
Any disagreement is printed genome by genome."""
import sys, glob, os, re
import numpy as np
ATOMS = "/Users/martin/Documents/CL1/PartnersRivals/Cannon/atoms512.bin"
A = np.fromfile(ATOMS, dtype=np.uint8).reshape(512, 65536)
d = sys.argv[1]
tot = bad = 0
for fn in sorted(glob.glob(os.path.join(d, "B*.bin")), key=lambda s: int(re.findall(r"B(\d+)", s)[0])):
    t = int(re.findall(r"B(\d+)", fn)[0]); ipt = t - 1
    b = np.fromfile(fn, dtype=np.uint8); assert b.size == 65536
    diff = np.nonzero(b != A[ipt])[0]
    tot += 65536; bad += diff.size
    cnt = np.bincount(b, minlength=8)
    print("ipt %3d  mismatches %5d   counts 000..111: %s" % (ipt, diff.size, " ".join("%5d" % c for c in cnt)))
    for g in diff[:10]:
        print("     genome %5d  wfc %d%d%d  atoms512 %d%d%d" % (g, b[g] >> 2, (b[g] >> 1) & 1, b[g] & 1,
                                                          A[ipt, g] >> 2, (A[ipt, g] >> 1) & 1, A[ipt, g] & 1))
print("TOTAL: %d genome-games compared, %d disagree" % (tot, bad))
