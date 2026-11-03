import pickle, numpy as np
VROWS = pickle.load(open('vrows.pkl','rb'))
FULL8 = ['0000','1000','0110','1110','0001','1001','0111','1111']
rng = np.random.default_rng(0); M = 6_000_000
r = np.sqrt(rng.random(M)); a = 2*np.pi*rng.random(M)
k = 4.0
rho = k*r/np.sqrt(1-r**2)
U, V = rho*np.cos(a), rho*np.sin(a)
def nashmask(g):
    ok = np.ones(M, bool)
    for aa,bb,cc in VROWS[g]:
        if (aa,bb,cc) != (0,0,0):
            ok &= (float(aa)*(U-1) + float(bb)*V <= float(cc) + 1e-12)
    return ok
cnt = sum(nashmask(g).astype(np.int8) for g in FULL8)
print('disk (k=4) shares of the sixteens count, from the 8 full-dimensional regions:')
for j in range(9):
    m = (cnt == j)
    if m.any(): print(f'   N={j}: {100*m.mean():6.2f} %   (sibling README: see table)')
print('   N=8 set equals {u<=-1, v<=-1}? mismatch share %.4f %%'
      % (100*np.mean((cnt==8) ^ ((U<=-1)&(V<=-1)))))
print('   N=0 set equals {u>0,v>0} minus [0,1]^2 ? mismatch share %.4f %%'
      % (100*np.mean((cnt==0) ^ (((U>0)&(V>0)) & ~((U<=1)&(V<=1))))))
dim2 = (U >= 1) & (V >= 0)
print('\ncontinuum: dim 3 on %.2f %% of the disk, dim 2 on %.2f %%' % (100*(~dim2).mean(), 100*dim2.mean()))
print('   dim-2 region inside the sixteens no-equilibrium region?  dim2 \\ {N=0} = %.4f %% of the disk'
      % (100*np.mean(dim2 & (cnt != 0))))
print('   {N=0} \\ dim2 = %.2f %% of the disk' % (100*np.mean((cnt == 0) & ~dim2)))
w = dim2 & (cnt != 0)
if w.any():
    print('   the exception: u range %.3f..%.3f  v range %.3f..%.3f'
          % (U[w].min(), U[w].max(), V[w].min(), V[w].max()))
sat = (U <= 1) & (V <= 0)
print('\n   saturation of the graded field V1=1 on {u<=1,v<=0}: %.2f %% of the disk' % (100*sat.mean()))
print('   the sixteens maximum count 8 on {u<=-1,v<=-1}:        %.2f %% of the disk' % (100*np.mean((U<=-1)&(V<=-1))))
print('   both facets carry a 3-dim piece, {u<0}:               %.2f %% of the disk' % (100*np.mean(U<0)))
