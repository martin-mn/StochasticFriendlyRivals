"""Inputs for cover.c: the co-players each run used (from out/<cls>.summary) and the survivor patterns."""
import os
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'out')
os.makedirs(os.path.join(OUT, 'cover'), exist_ok=True)
for tag, names in [('sep', ['c6', 'c8', 'c68']), ('p9a', ['c9a']), ('p9b', ['c9b'])]:
    Q = set(); P = []
    for n in names:
        for l in open(os.path.join(OUT, n + '.summary')):
            a = l.split()
            if a and a[0] == 'q': Q.add(int(a[1]))
        P += [l.split('\t')[1] for l in open(os.path.join(OUT, n + '.tsv'))]
    open(os.path.join(OUT, 'cover', tag + '.q'), 'w').write('\n'.join(map(str, sorted(Q))) + '\n')
    open(os.path.join(OUT, 'cover', tag + '.pat'), 'w').write('\n'.join(P) + '\n')
    print(tag, len(Q), 'co-players', len(P), 'strata')
