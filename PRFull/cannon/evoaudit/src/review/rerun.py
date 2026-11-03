import sys, json, glob, numpy as np
KIT='/Users/martin/Documents/CL1/PRFull/cannon/evoaudit'; sys.path.insert(0, KIT+'/src')
import audit as au, load
S=np.load(KIT+'/sample/sample.npy'); G=load.games(KIT+'/games.dat')
K={}
for fn in glob.glob(KIT+'/out_probe/chunk_*.jsonl'):
    for l in open(fn): r=json.loads(l); K[r['sid']]=r
for sid in map(int, sys.argv[1:]):
    r=S[sid]; g=G[int(r['task'])-1]
    o=au.audit(r['p'], g['cs'], g['ct'], prod_atom=int(r['atom']))
    k=K[sid]; diffs=[]
    for key in o:
        if key in ('sec',): continue
        if key=='pr':
            # production doubles depend on the compiler (clang here vs GCC 8.5 on Cannon)
            d=[f for f in o['pr'] if o['pr'][f]!=k['pr'][f]]
            if d: diffs.append('pr:'+','.join(d))
            continue
        if o[key]!=k[key]: diffs.append(key)
    print(sid, 'val', o['val'], 'cannon', k['val'], 'differing fields:', diffs or 'none')
