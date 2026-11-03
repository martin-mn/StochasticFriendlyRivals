"""Build the unit list of the stratcert kit: the retained friendly-rival support strata of W, N, E and the
ray T>S, with the census points stored for each.

Inputs (all regenerated or copied here, read-only sources):
  cls.bin            exact efficiency class of all 3^16 open strata  (./ub classes cls.bin 3;  ub.c = check/ubound/ub.c)
  c2.idx/.qp/.qm     the 1681-co-player exclusion cover  (./ub cover cls.bin check/ubound/results/pool_c2_q.txt c2 3)
                     -> a stratum is RETAINED for (sign) when no co-player of the pool excludes it (q = -1)
  census_pts_PA.txt  the census points of PartnersRivalsM2f (Cannon ~/CL1/PartnersRivalsM2f/ca1/out/pt_*.txt,
                     array 45751662), lines with sign P (T>S) of any class and sign M (T<S) of class A,
                     prefixed with file name and line number
  M2f/check/st_eff.txt  the census stratum list (for the cross-check only)

Output: ../units.tsv, one line per UNIT = (stratum, sign):
  uid  sign  pattern  ub_index  census_code  dim  class  loci  npts  src
  loci: comma list among W,N,ray+ (sign +) or E (sign -); a sign-+ unit of class C is in W and ray+, of class A in
  N and ray+, of class Z in ray+ only (competitiveness depends only on the sign of T-S; efficiency on the switch
  line is w_DD = 0, i.e. class C, A or Z).
and ../points.tsv: uid  k  source  p0 .. p15 (the %.17g doubles as printed by the census; exact binary rationals)
"""
import numpy as np, collections, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
N3 = 3**16
cls = np.fromfile(os.path.join(HERE, 'cls.bin'), dtype=np.uint8)
idx = np.fromfile(os.path.join(HERE, 'c2.idx'), dtype=np.uint32).astype(np.int64)
qp = np.fromfile(os.path.join(HERE, 'c2.qp'), dtype=np.int32)
qm = np.fromfile(os.path.join(HERE, 'c2.qm'), dtype=np.int32)
CL = {0: '-', 1: 'C', 2: 'A', 3: 'Z'}

def digits(i):
    d = []
    for _ in range(16): d.append(i % 3); i //= 3
    return d
def pattern(i): return ''.join('01*'[x] for x in digits(i))            # ub encoding: 0 pinned 0, 1 pinned 1, 2 free
def conv_census(g):                                                      # census digits 0 -> 0, 1 -> free, 2 -> one
    r = 0
    for j, dg in enumerate(digits(g)): r += {0: 0, 1: 2, 2: 1}[dg]*3**j
    return r
def census_code(i):
    r = 0
    for j, dg in enumerate(digits(i)): r += {0: 0, 2: 1, 1: 2}[dg]*3**j
    return r
def dim(i): return digits(i).count(2)

c = cls[idx]
plus = idx[(qp < 0) & np.isin(c, [1, 2, 3])]
minusA = idx[(qm < 0) & (c == 2)]
sets = {'W': [int(i) for i in plus if cls[i] == 1], 'N': [int(i) for i in plus if cls[i] == 2],
        'E': [int(i) for i in minusA], 'ray+': [int(i) for i in plus]}
print('retained strata (not excluded by the 1681 co-players of check/ubound/results/pool_c2_q.txt):')
for k, v in sets.items():
    h = np.bincount([dim(i) for i in v], minlength=9)
    print('  %-5s %5d  by dim %s' % (k, len(v), h.tolist()))
w = sets['W']
f00 = [i for i in w if digits(i)[6] == 0 and digits(i)[8] == 0]; f11 = [i for i in w if digits(i)[6] == 1 and digits(i)[8] == 1]
print('  W by (p6,p8): (0,0) %d  (1,1) %d  other %d' % (len(f00), len(f11), len(w) - len(f00) - len(f11)))
print('  ray+ = W u N u Z-only: %d = %d + %d + %d' % (len(sets['ray+']), len(w), len(sets['N']),
      sum(1 for i in sets['ray+'] if cls[i] == 3)))
nmir = 0
def mirror_idx(i):          # C/D mirror p_s -> 1 - p_{15-s}: pattern digit 0 <-> 1, free stays, positions reversed
    d = digits(i); r = 0
    for s in range(16):
        x = d[15 - s]; r += {0: 1, 1: 0, 2: 2}[x]*3**s
    return r
print('  mirror(N) == E as sets:', set(mirror_idx(i) for i in sets['N']) == set(sets['E']))

# census points
pts = collections.defaultdict(list)
for l in open(os.path.join(HERE, 'census_pts_PA.txt')):
    f = l.split(); fn, ln, sg, g, cl = f[0], f[1], f[2], int(f[3]), f[4]; p = f[5:]
    assert len(p) == 16
    i = conv_census(g); s = +1 if sg == 'P' else -1
    pts[(i, s)].append(('%s:%s' % (fn, ln), cl, p))
# st_eff cross-check
cen = {'W': set(), 'N': set(), 'E': set(), 'ray+': set()}
for l in open('/Users/martin/Documents/CL1/PartnersRivalsM2f/check/st_eff.txt'):
    if not l.strip(): continue
    g, cl, nP, nM, K = l.split(); i = conv_census(int(g))
    if cl == 'C' and int(nP) > 0: cen['W'].add(i)
    if cl == 'A' and int(nP) > 0: cen['N'].add(i)
    if cl == 'A' and int(nM) > 0: cen['E'].add(i)
    if int(nP) > 0: cen['ray+'].add(i)
for k in sets:
    S = set(sets[k]); print('  %-5s census FR strata %5d, all retained: %s, retained not in census: %s'
                           % (k, len(cen[k]), cen[k] <= S, [pattern(i) for i in sorted(S - cen[k])]))

units = [(i, +1) for i in sorted(sets['ray+'])] + [(i, -1) for i in sorted(sets['E'])]
nopt = []; clsbad = 0; npts = 0
with open(os.path.join(HERE, '..', 'units.tsv'), 'w') as fu, open(os.path.join(HERE, '..', 'points.tsv'), 'w') as fp:
    fu.write('#uid\tsign\tpattern\tub_index\tcensus_code\tdim\tclass\tloci\tnpts\n')
    fp.write('#uid\tk\tsource\tp0..p15 (census doubles, %.17g)\n')
    for u, (i, s) in enumerate(units):
        cl = CL[int(cls[i])]
        loci = ('W,ray+' if cl == 'C' else 'N,ray+' if cl == 'A' else 'ray+') if s > 0 else 'E'
        P = pts.get((i, s), [])
        for k, (src, pcl, p) in enumerate(P):
            if pcl != cl: clsbad += 1
            fp.write('%d\t%d\t%s\t%s\n' % (u, k, src, '\t'.join(p)))
        npts += len(P)
        if not P: nopt.append((u, pattern(i), s))
        fu.write('%d\t%+d\t%s\t%d\t%d\t%d\t%s\t%s\t%d\n' % (u, s, pattern(i), i, census_code(i), dim(i), cl, loci, len(P)))
print('units: %d (sign + %d, sign - %d); census points attached %d; class disagreements census vs exact %d'
      % (len(units), len(sets['ray+']), len(sets['E']), npts, clsbad))
print('units with no stored census point:', nopt)
unused = sum(len(v) for (i, s), v in pts.items() if not ((s > 0 and i in set(sets['ray+'])) or (s < 0 and i in set(sets['E']))))
print('census points of these kinds not attached to a unit:', unused)
