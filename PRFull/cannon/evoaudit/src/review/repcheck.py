import os, sys, numpy as np, filecmp
KIT='/Users/martin/Documents/CL1/PRFull/cannon/evoaudit'
sys.path.insert(0, KIT+'/src'); import load
M2='/Users/martin/Documents/CL1/PartnersRivalsM2f/check'
src={'F1':M2+'/wfc_a333','F2':M2+'/wfc_f2'}
for reg in ('F1','F2'):
    nsame=0; ndiff=0; badrep=0; nrep=0; fields=None
    for task in range(1,513):
        a=os.path.join(KIT,'data',reg,'r_1o3_%d_1.txt'%task); b=os.path.join(src[reg],'r_1o3_%d_1.txt'%task)
        if filecmp.cmp(a,b,shallow=False): nsame+=1
        else: ndiff+=1
        d=load.dump(reg,task)
        for line in open(b):
            if not line.startswith('REP'): continue
            t=line.split(); irep=int(t[9]); nsamp=int(t[10]); ndist=int(t[11]); acnt=[int(x) for x in t[12:20]]
            m=d['irep']==irep; nrep+=1
            cnt=np.bincount(d['atom'][m], weights=d['count'][m], minlength=8).astype(int)
            if list(cnt)!=acnt or m.sum()!=ndist or d['count'][m].sum()!=nsamp: badrep+=1; print(reg,task,irep,list(cnt),acnt,m.sum(),ndist)
    print(reg,'r files identical',nsame,'differ',ndiff,'replicates',nrep,'mismatching atom counts',badrep)
