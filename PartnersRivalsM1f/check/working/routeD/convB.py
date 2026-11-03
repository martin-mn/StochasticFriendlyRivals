"""Convention B: families p(eps).  Report both deficits at a ladder of eps."""
import mp as M, mpmath as mp
mp.mp.dps=60
uW,vW=mp.mpf(-1)/2,mp.mpf(1)/2          # W:  s = 1-u-v = 1,  T>S
uS,vS=mp.mpf(1)/2,mp.mpf(-1)/2          # S:  s = 1,          T<S
uN,vN=mp.mpf('0.2'),mp.mpf('1.5')       # N
uE,vE=mp.mpf('1.5'),mp.mpf('-0.2')      # E
KS=[4,6,8,10,12,16,20]

def report(tag,fam,u,v,sign,ks=KS):
    print(tag)
    for k in ks:
        eps=mp.mpf(10)**(-k)
        p=fam(eps)
        e=M.edef(p,u,v,eps); r=M.rm(p,sign,eps)
        print(f'   eps=1e-{k:<3d} edef={mp.nstr(e,6):>12}  rm={mp.nstr(r,6):>12}  '
              f'p=[{", ".join(mp.nstr(mp.mpf(t),5) for t in p)}]')
    print()

print('############ W  (u,v)=(-1/2,1/2), s=1, sign=+1')
for c in ('0','0.1','0.3','0.5'):
    for th in ('0.25','0.5','0.75','1','1.5'):
        C=mp.mpf(c); T=mp.mpf(th)
        report(f'--- B1 c={c} theta={th}:  p=(1-eps, c, c, eps^theta)',
               lambda eps,C=C,T=T: [1-eps, C, C, eps**T], uW,vW,1)
print('############ W  ratio family p=(1-kappa*eps, 1/2,1/2, lambda*eps)  [theta=phi=1]')
for kap,lam in (('1','1'),('0','1'),('1','1000'),('1000','1')):
    K=mp.mpf(kap); L=mp.mpf(lam)
    report(f'--- B2 kappa={kap} lambda={lam}',
           lambda eps,K=K,L=L: [1-K*eps, mp.mpf(1)/2, mp.mpf(1)/2, L*eps], uW,vW,1, [6,8,10,12])
print('############ W  c>1/2 must fail rivalry')
report('--- B1 c=0.7 theta=0.5', lambda eps: [1-eps, mp.mpf('0.7'), mp.mpf('0.7'), mp.sqrt(eps)], uW,vW,1,[6,10,14])
print('############ W  the p2=1-c*eps^theta edge family of rivals.md 6(c) (rivalry only)')
for th in ('0.5','1','2'):
    T=mp.mpf(th)
    report(f'--- B3 theta={th}: p=(1/2, 1-eps^theta, 0, 0)',
           lambda eps,T=T: [mp.mpf(1)/2, 1-eps**T, 0, 0], uW,vW,1,[6,10,14])
print('############ S: the limit points of FR_S that convention A misses')
for ps in ('0.5','1','2'):
    P=mp.mpf(ps)
    report(f'--- B4 psi={ps}: p=(1, eps^psi, 1, 1/2)   -> limit point (1,0,1,1/2)',
           lambda eps,P=P: [1, eps**P, 1, mp.mpf(1)/2], uS,vS,-1,[6,10,14])
    report(f'--- B5 psi={ps}: p=(1, 1, eps^psi, 1/2)   -> limit point (1,1,0,1/2)',
           lambda eps,P=P: [1, 1, eps**P, mp.mpf(1)/2], uS,vS,-1,[6,10,14])
print('############ N and E: the convention-B friendly rivals')
for th in ('0.5','1','2'):
    T=mp.mpf(th)
    report(f'--- B6 N theta={th}: p=(1/2, 0, 1, eps^theta) -> limit (1/2,0,1,0) in Lambda_N',
           lambda eps,T=T: [mp.mpf(1)/2, 0, 1, eps**T], uN,vN,1,[6,10,14])
    report(f'--- B7 E theta={th}: p=(1-eps^theta, 0, 1, 1/2) -> limit (1,0,1,1/2) in Lambda_E',
           lambda eps,T=T: [1-eps**T, 0, 1, mp.mpf(1)/2], uE,vE,-1,[6,10,14])
