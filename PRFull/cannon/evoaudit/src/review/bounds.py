import sys, math, numpy as np, collections, time
KIT='/Users/martin/Documents/CL1/PRFull/cannon/evoaudit'; sys.path.insert(0, KIT+'/src')
import sample as sm, load
t=time.time()
G, A, h, c, Cr, hs, Nh, Ch, piG, piA, piP, pi = sm.population()
print('population %.1fs'%(time.time()-t), len(h))
# independent recomputation of pi
wid=A['wedge'].astype(int); reg=A['reg'].astype(int); atm=A['atom'].astype(int)
ngw=collections.Counter(g['wedge'] for g in G)
bw=np.array([ {0:3000,1:1000,2:3000,3:3000}[w] for w in range(4)])[wid]
Cw=np.array([40000*ngw['WSNE'[w]] for w in range(4)])[wid]
pg=np.minimum(1, bw*c/Cw)
ah=np.where(np.isin(atm,[6,7,3,1]),150.,50.)
key=(reg*32+wid*8+atm)
u,inv,cnt=np.unique(key,return_inverse=True,return_counts=True)
Csum=np.bincount(inv,weights=c)
pa=np.minimum(1, ah/cnt[inv]); pp_=np.minimum(1, ah*c/Csum[inv])
pi2=1-(1-pg)*(1-pa)*(1-pp_)
print('max |pi - pi2|', np.abs(pi-pi2).max())
rng=np.random.default_rng(sm.SEED); sel=np.where(rng.random(len(h))<pi)[0]; print('drawn',len(sel))
S=np.load(KIT+'/sample/sample.npy'); print('sample matches', len(S)==len(sel) and (S['task']==A['task'][sel]).all() and (S['row']==A['row'][sel]).all() and (S['reg']==A['reg'][sel]).all())
# p in sample equals dump p
bad=0
for i in range(0,len(S),997):
    d=load.dump('F%d'%S['reg'][i], int(S['task'][i])); bad += not np.array_equal(d['p'][S['row'][i]], S['p'][i])
print('spot check p (every 997th):', 'bad', bad)
L95=-math.log(0.05)
lg=np.where(pi<1, -np.log1p(-np.minimum(pi, np.nextafter(1,0))), np.inf)
def kap(mask,Cg):
    m=mask&(pi<1); return (lg[m]*Cg/c[m]).min()
for rg in (1,2):
    k=kap(reg==rg,512*40000); print('F%d'%rg,'bound %.7g%%'%(100*L95/k))
    for w in 'WSNE':
        m=(reg==rg)&(wid=='WSNE'.index(w)); k=kap(m,40000*ngw[w]); print('  ',w,'bound %.7g%%'%(100*L95/k), 'games',ngw[w])
# W-111 records: p0, p15
for rg in (1,2):
    m=np.where((S['reg']==rg)&(S['wedge']==0)&(S['atom']==7))[0]
    print('F%d W-111 sampled'%rg, len(m), 'p0==1', int((S['p'][m,0]==1).sum()), 'p15==0', int((S['p'][m,15]==0).sum()), 'pi==1', int((S['pi'][m]==1).sum()))
    mp=(reg==rg)&(wid==0)&(atm==7); print('   population W-111', int(mp.sum()))
