#!/usr/bin/env python3
"""Aggregate out/task_*.txt of the full bingap array.  Run ONLY after
   sacct -j $(cat JOBID) --format=JobID,State,ExitCode -X --noheader | awk '$2 != "COMPLETED"'
prints nothing, and check that every out/task_*.err is empty."""
import glob, re, sys
from fractions import Fraction as F
rows = []
for fn in sorted(glob.glob('out/task_*.txt')):
    s = open(fn).read()
    m = re.search(r'task (\d+) of (\d+): pairs (\d+) zero (\d+) flagged (\d+) symmetry_violations (\d+) below_1/37 (\d+) '
                  r'below_1e-9 (\d+) min_nonzero (\d+)/(\d+) at p=(-?\d+) q=(-?\d+) max_Dstar (\d+) at p=(-?\d+) q=(-?\d+)', s)
    if not m: sys.exit('unparsable ' + fn)
    rows.append([int(x) for x in m.groups()])
T = rows[0][1]
assert len(rows) == T and sorted(r[0] for r in rows) == list(range(T)), 'missing tasks'
tot = sum(r[2] for r in rows); expect = 65536*65537//2
print('pairs %d (expected %d: %s)' % (tot, expect, tot == expect))
print('exact zeros %d, flagged %d, symmetry violations %d' % (sum(r[3] for r in rows), sum(r[4] for r in rows), sum(r[5] for r in rows)))
print('pairs (q >= p) with 0 < |m| < 1/37: %d ; with 0 < |m| < 1e-9: %d' % (sum(r[6] for r in rows), sum(r[7] for r in rows)))
mn = min((F(r[8], r[9]), r[10], r[11]) for r in rows)
print('smallest nonzero |w_DC - w_CD| over all binary pairs: %s = %.6g at p=%d q=%d' % (mn[0], float(mn[0]), mn[1], mn[2]))
mx = max((r[12], r[13], r[14]) for r in rows)
print('largest D* (number of minimum-order spanning trees): %d at p=%d q=%d' % mx)
