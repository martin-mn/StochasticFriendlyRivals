#!/usr/bin/env python3
"""rcheck_compare.py -- compare a replay of the efficiency test (rcheck.c) with the deposited per-task counts.

    python3 validate/rcheck_compare.py [--replay DIR] [--ref DIR] [--sweep DIR] TASK [TASK ...]

For each TASK, DIR/replay_<TASK>.out is the output of `rcheck <TASK> 3190000 classify` (validate/rcheck.sh writes it,
default DIR validate/build).  It is compared with
  - the task's RC line of rc.tsv (--ref, default: this script's directory), field for field;
  - the task's FLAG lines of flags.tsv (same directory; the header line is skipped), in order;
  - the read counts n, effCC, effALT, effZ of the RC line against the ALPHA line of the Cannon run
    (array 48928780), --sweep/out_13/a_1o3_<TASK>.tsv (default: the kit, one directory up), which must carry alpha 1/3.
The columns of rc.tsv (no header): RC task n readC readA readZ exactC exactA exactZ flagged W S N E, as described
in rcheck.c and README.md.  Prints one line per difference and the summary
    replayed N tasks: RC lines differ 0, FLAG lines differ 0, ALPHA counts differ 0
Exit status 0 when nothing differs, 1 otherwise.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def lines_of(path):
    with open(path) as fh:
        return [l.rstrip('\n') for l in fh if l.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--replay', default=os.path.join(HERE, 'build'), help='directory of the replay_<task>.out files')
    ap.add_argument('--ref', default=HERE, help='directory of the deposited rc.tsv and flags.tsv')
    ap.add_argument('--sweep', default=os.path.dirname(HERE), help='the sweep13 kit, whose out_13/ holds the ALPHA lines')
    ap.add_argument('tasks', nargs='+', type=int)
    a = ap.parse_args()
    rc = {}
    for l in lines_of(os.path.join(a.ref, 'rc.tsv')):
        f = l.split('\t')
        assert f[0] == 'RC' and len(f) == 14, l
        rc[int(f[1])] = l
    flags = {}
    fl = lines_of(os.path.join(a.ref, 'flags.tsv'))
    assert fl[0].startswith('tag\ttask\t'), fl[0]
    for l in fl[1:]:
        f = l.split('\t')
        assert f[0] == 'FLAG', l
        flags.setdefault(int(f[1]), []).append(l)
    nrc = nflag = nalpha = 0
    for t in a.tasks:
        out = lines_of(os.path.join(a.replay, 'replay_%d.out' % t))
        mine = [l for l in out if l.startswith('RC\t')]
        assert len(mine) == 1 and len(out) == 1 + sum(l.startswith('FLAG\t') for l in out), 'task %d: malformed replay' % t
        if mine[0] != rc.get(t):
            nrc += 1
            print('task %d: RC line DIFFERS: replay %r, rc.tsv %r' % (t, mine[0], rc.get(t)))
        fr, fd = [l for l in out if l.startswith('FLAG\t')], flags.get(t, [])
        nd = sum(x != y for x, y in zip(fr, fd)) + abs(len(fr) - len(fd))
        if nd:
            nflag += nd
            print('task %d: %d FLAG line(s) DIFFER (replay %d, flags.tsv %d)' % (t, nd, len(fr), len(fd)))
        al = lines_of(os.path.join(a.sweep, 'out_13', 'a_1o3_%d.tsv' % t))
        assert len(al) == 1 and al[0].split('\t')[:2] == ['ALPHA', '1/3'], 'task %d: not one ALPHA line at 1/3' % t
        got, want = mine[0].split('\t')[2:6], al[0].split('\t')[2:6]       # n effCC effALT effZ
        if got != want:
            nalpha += 1
            print('task %d: read counts DIFFER from the ALPHA line: replay %s, Cannon %s' % (t, got, want))
    print('replayed %d tasks: RC lines differ %d, FLAG lines differ %d, ALPHA counts differ %d' % (len(a.tasks), nrc, nflag, nalpha))
    return 1 if nrc or nflag or nalpha else 0


if __name__ == '__main__':
    sys.exit(main())
