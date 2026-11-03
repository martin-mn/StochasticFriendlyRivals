import csv, numpy as np, sys
rows = list(csv.DictReader(open('../Text/SIData1_games512.csv')))
assert len(rows) == 512 and [int(r['game']) for r in rows] == list(range(512))
ORDER = ['000','100','010','001','110','011','111']
reg = {'all': lambda r: True, 'PD': lambda r: r['quadrant']=='PD','SD': lambda r: r['quadrant']=='SD',
       'W': lambda r: r['wedge']=='W','S': lambda r: r['wedge']=='S','N': lambda r: r['wedge']=='N','E': lambda r: r['wedge']=='E',
       'above': lambda r: r['wedge'] in 'NE'}
for p in ['m1_N1000','m2_N1000','m1_N100','m2_N100']:
    E = np.array([float(r[p+'_efficiency']) for r in rows])
    print(p, 'eff min %.4f median %.4f max %.4f <0.9 %d' % (E.min(), np.median(E), E.max(), (E<0.9).sum()))
    for lab in ('most_abundant','most_enriched'):
        c = {a: sum(r[p+'_'+lab]==a for r in rows) for a in ORDER}
        print('   ', lab, {a:n for a,n in c.items() if n})
    for rn in ('PD','S','W','N','E','above'):
        sub = [r for r in rows if reg[rn](r)]
        s111 = np.array([float(r[p+'_share_111']) for r in sub]); comp = np.array([float(r[p+'_share_competitive']) for r in sub])
        lead = {a: sum(r[p+'_most_abundant']==a for r in sub) for a in ORDER}
        print('    %-5s n=%3d median s111 %.4g  median comp %.4g  lead %s' % (rn, len(sub), np.median(s111), np.median(comp), {a:n for a,n in lead.items() if n}))
    comp = np.array([float(r[p+'_share_competitive']) for r in rows]); print('    all median comp %.4g' % np.median(comp))
    nrep = np.array([int(r[p+'_reps_naming_most_abundant']) for r in rows]); print('    all ten replicates name leader at', (nrep==10).sum())
    if p.startswith('m2'):
        print('    undetermined', sum(r[p+'_most_enriched_determined']=='0' for r in rows))
    # shares sum
    tot = np.array([sum(float(r[p+'_share_'+a]) for a in ORDER) for r in rows]); print('    share sums %.6f..%.6f' % (tot.min(), tot.max()))
for i in (0, 117, 184, 191, 246, 337, 396, 237, 326, 415, 500):
    r = rows[i]
    print(i, r['u'], r['v'], r['quadrant'], r['wedge'], r['dist_switch_line'], r['dist_T_eq_S_line'], '| m1', r['m1_N1000_most_abundant'], r['m1_N1000_share_111'], '| m2', r['m2_N1000_most_abundant'], r['m2_N1000_most_enriched'], r['m2_N1000_share_111'], r['m2_N1000_most_enriched_determined'])
