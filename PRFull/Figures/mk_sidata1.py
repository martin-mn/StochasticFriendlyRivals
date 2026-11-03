"""SI Data 1: one row per sampled game, the per-game summaries of the four cube runs, from fig45.py's loaders (read-only)."""
import sys, os, io, time, csv
sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
t0 = time.time()
sys.path.insert(0, '/Users/martin/Documents/CL1/PRFull/Figures')
import fig45 as F
sys.stdout = sys.__stdout__
import numpy as np
ORDER = F.ORDER
buf = io.StringIO(); so = sys.stdout; sys.stdout = buf
out = {run: (F.load_m1(run), F.load_m2(run)) for run in ('F1', 'F2')}
sys.stdout = so
print('loaded in %.1f s' % (time.time() - t0))
UV = F.UV; u, v = UV[:, 0], UV[:, 1]
disk = UV / np.sqrt(16 + u * u + v * v)[:, None]
d_sw = (u + v - 1) / np.sqrt(2)      # signed distance to the switch line, > 0 above it
d_ts = (u - v - 1) / np.sqrt(2)      # signed distance to T = S, > 0 where T < S
# kernel masses: the same for both population sizes at a memory
for sp, k in (('m1', 0), ('m2', 1)):
    a, b = out['F1'][k]['meas'], out['F2'][k]['meas']
    assert np.array_equal(a, b), sp
RUNS = [('m1', 'N1000', out['F1'][0]), ('m2', 'N1000', out['F1'][1]), ('m1', 'N100', out['F2'][0]), ('m2', 'N100', out['F2'][1])]
hdr = ['game', 'u', 'v', 'quadrant', 'wedge', 'disk_x', 'disk_y', 'dist_switch_line', 'dist_T_eq_S_line']
for sp, n, d in RUNS:
    p = '%s_%s_' % (sp, n)
    hdr += [p + 'efficiency', p + 'most_abundant', p + 'reps_naming_most_abundant', p + 'most_enriched']
    if sp == 'm2':
        hdr += [p + 'most_enriched_determined']
    hdr += [p + 'share_' + a for a in ORDER] + [p + 'share_efficient', p + 'share_stable', p + 'share_competitive']
for sp, k in (('m1', 0), ('m2', 1)):
    hdr += ['%s_kernel_mass_%s' % (sp, a) for a in ORDER]
fn = sys.argv[1]
g = lambda x: '%.6g' % x
with open(fn, 'w', newline='') as fh:
    w = csv.writer(fh, lineterminator='\n'); w.writerow(hdr)
    for i in range(512):
        row = [i, repr(float(u[i])), repr(float(v[i])), F.QUAD[i], F.WEDGE[i], '%.6f' % disk[i, 0], '%.6f' % disk[i, 1],
               '%.6f' % d_sw[i], '%.6f' % d_ts[i]]
        for sp, n, d in RUNS:
            nrep = int((d['R'][i].argmax(1) == d['wa'][i]).sum())
            row += [g(d['E'][i]), ORDER[d['wa'][i]], nrep, ORDER[d['we'][i]]]
            if sp == 'm2':
                und = {k for k, _, _ in d['und']}; row += [int(i not in und)]
            row += [g(x) for x in d['S'][i]] + [g(x) for x in d['PS'][i]]
        for sp, k in (('m1', 0), ('m2', 1)):
            row += [g(x) for x in out['F1'][k]['meas'][i]]
        w.writerow(row)
print('wrote', fn, os.path.getsize(fn), 'bytes,', len(hdr), 'columns')
