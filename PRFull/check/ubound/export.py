"""export the upper-bound certificate (12 co-players, per-stratum assignment) and the count tables"""
import gzip, numpy as np, status as st
idx, qp, qm = st.load('runs/x'); M, d = st.masks(idx); cls, dim = st.base()
T = np.fromfile('runs/tri_star.bin', dtype=np.int64).reshape(-1, 3)
res = np.fromfile('runs/ver_star.res', dtype=np.int8); assert len(res) == len(T) and (res == 1).all()
Q = [int(x) for x in open('qstar.txt').read().split()]
def gstr(g): return ''.join(str((g >> s) & 1) for s in range(16))
CL = '-CAZ'
lines = ['# upper-bound certificate: every (stratum, sign) below is EXCLUDED by the co-player q: the margin',
         '# sign*(w_DC - w_CD) against q is < 0 at every point of the open stratum (eps -> 0 limit).',
         '# columns: stratum (char s = p_s: 0, 1, * free), dim, exact efficiency class, sign, q (genome), q as 16 bits (char s = q_s)']
order = np.lexsort((T[:, 0], -T[:, 1]))
for i, sg, q in T[order]:
    lines.append('%s\t%d\t%s\t%+d\t%d\t%s' % (st.pattern(int(i)), dim[int(i)], CL[cls[int(i)]], sg, q, gstr(int(q))))
with gzip.open('results/ub_certificate.tsv.gz', 'wt') as f: f.write('\n'.join(lines) + '\n')
out = []
out.append('THE 12 CERTIFICATE CO-PLAYERS (genome; bits q_0..q_15, q_s = 1 means C at state s = 4x+y)')
for q in Q:
    n = {+1: int(((T[:, 2] == q) & (T[:, 1] == 1)).sum()), -1: int(((T[:, 2] == q) & (T[:, 1] == -1)).sum())}
    out.append('  %5d  %s   certifies %7d (sign +1)  %6d (sign -1)' % (q, gstr(q), n[1], n[-1]))
out.append('')
out.append('TARGETS PER LOCUS AND DIMENSION (all excluded; the certificate co-player of each is in ub_certificate.tsv.gz)')
for k in ('W', 'N', 'E', 'ray+', 'W68'):
    m, sg, b = M[k]; t = m & (d >= b); h = np.bincount(d[t], minlength=17)
    qq = qp if sg > 0 else qm; assert (qq[t] >= 0).all()
    out.append('  %-5s sign %+d, dim >= %2d: %7d strata  by dim %s' % (k, sg, b, t.sum(), {j: int(h[j]) for j in range(b, 17) if h[j]}))
open('results/certificate_summary.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
