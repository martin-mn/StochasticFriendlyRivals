import pickle, numpy as np
from fractions import Fraction as F
VROWS = pickle.load(open('vrows.pkl','rb'))
KEYS = list(VROWS)
FULL8 = ['0000','1000','0110','1110','0001','1001','0111','1111']
def masks(L, N):
    uu = np.linspace(-L, L, N); U, V = np.meshgrid(uu, uu, indexing='ij')
    M = {}
    for g in KEYS:
        ok = np.ones_like(U, bool)
        for a,b,c in VROWS[g]:
            if (a,b,c) != (0,0,0):
                ok &= (float(a)*(U-1) + float(b)*V <= float(c) + 1e-12)
        M[g] = ok
    return U, V, M
U, V, M = masks(8.0, 1601)
cnt8 = sum(M[g] for g in FULL8)
cntall = sum(M[g] for g in KEYS)
print('per-strategy share of [-8,8]^2 where limit-Nash (weak):')
for g in KEYS: print(f'   {g} {100*M[g].mean():7.3f} %')
print('\nN=8 (of the 8 full-dimensional) equals {u<=0, v<=0}?',
      np.array_equal(cnt8==8, (U<=1e-12)&(V<=1e-12)))
for cond, nm in [(((U<=1e-12)&(V<=1e-12)), 'u<=0,v<=0'),
                 (((U<=1+1e-12)&(V<=1e-12)), 'u<=1,v<=0')]:
    print(f'   N=8 vs {nm}: agree {100*np.mean((cnt8==8)==cond):7.3f} %')
d = (cnt8==8) ^ ((U<=1e-12)&(V<=1e-12))
if d.any():
    print('   disagreement sample (u,v):', list(zip(U[d][:5].round(3), V[d][:5].round(3))))
dim2 = (U >= 1-1e-12) & (V >= -1e-12)
zero = (cntall == 0)
print('\ncontinuum dim-2 region vs the sixteen no-equilibrium region:')
print('   dim2 area %  %.3f   N=0 (all 16) area %%  %.3f   N=0 (full-dim 8) area %%  %.3f'
      % (100*dim2.mean(), 100*zero.mean(), 100*(cnt8==0).mean()))
print('   dim2 \\ N=0 :', 100*np.mean(dim2 & ~zero), '%   N=0 \\ dim2 :', 100*np.mean(zero & ~dim2), '%')
w = dim2 & ~zero
if w.any():
    print('   dim2 but some vertex is Nash, sample:', sorted(set(zip(U[w].round(2), V[w].round(2))))[:8])
# projective-disk shares (k=4 as in the sibling figure)
th = np.random.default_rng(0).random((4000000,2))
r = np.sqrt(th[:,0]); a = 2*np.pi*th[:,1]
k = 4.0
rho = k*r/(1-r+1e-300)           # disk radius r -> plane radius rho (k=4 projective disk)
Ud, Vd = rho*np.cos(a), rho*np.sin(a)
dd2 = (Ud >= 1) & (Vd >= 0)
print('\ndisk (k=4) shares:  continuum dim2 %.2f %%   dim3 %.2f %%' % (100*dd2.mean(), 100*(~dd2).mean()))
