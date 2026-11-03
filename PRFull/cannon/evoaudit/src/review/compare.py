import json, glob
KIT='/Users/martin/Documents/CL1/PRFull/cannon/evoaudit'
K={}
for fn in glob.glob(KIT+'/out_probe/chunk_*.jsonl'):
    for l in open(fn): r=json.loads(l); K[r['sid']]=r
for fn in ('indep_lib.jsonl','indep_full.jsonl'):
    for l in open(fn):
        o=json.loads(l); r=K[o['sid']]; C=r['C']; S=r['S']
        kC_neg = C['n_neg_i']+C['n_neg_x']; kC_zero=C['n_struct0']+C['n_zero_x']; Cfull = C['n_exact']==C['n_unres']
        kS_gain = S['n_gain_i']+S['n_struct_gain']+S['n_gain_x']; kS_tie=S['n_struct_tie']+S['n_tie_x']; Sfull = S['n_exact']==S['n_unres']
        line = '%5d %s %s prod %s kit %s indep %s' % (o['sid'], o['reg'], o['wedge'], o['prod'], r['val'], o['indep'])
        if 'full' in o: line += ' fullexact %s mirror_mismatch %d' % (o['full']['atom'], o['full']['mirror_mismatch'])
        line += ' | C neg kit %d%s indep %d, zero kit %d%s indep %d' % (kC_neg, '' if Cfull else '(partial)', o['C']['n_neg'], kC_zero, '' if Cfull else '(partial)', o['C']['n_zero'])
        line += ' | S gain kit %d%s indep %d, tie kit %d%s indep %d' % (kS_gain, '' if Sfull else '(partial)', o['S']['n_gain'], kS_tie, '' if Sfull else '(partial)', o['S']['n_tie'])
        line += ' | min kit %.6g indep %.6g, gmax kit %.6g indep %.6g' % (C['min_f'], o['C']['min_f'], S['max_f'], o['S']['gmax_f'])
        if 'full' in o: line += ' | full: nneg %d nzero %d ngain %d ntie %d mmin %.6g gmax %.6g' % tuple(o['full'][k] for k in ('nneg','nzero','ngain','ntie','mmin','gmax'))
        print(line)
