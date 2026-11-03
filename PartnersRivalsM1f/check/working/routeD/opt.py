"""D(eps) = sup_p d_eps(p): how close to a friendly rival can one get at error rate eps?"""
import mp as M, mpmath as mp, random
mp.mp.dps=50
u,v=mp.mpf(-1)/2, mp.mpf(1)/2      # W, s = 1-u-v = 1

def f(x,eps):
    p=[min(max(t,mp.mpf(0)),mp.mpf(1)) for t in x]
    return M.d_eps(p,u,v,eps)

def nelder(x0,eps,it=250,scale=None):
    n=4
    if scale is None: scale=[mp.mpf('0.05')]*n
    S=[list(x0)]
    for i in range(n):
        y=list(x0); y[i]=min(max(y[i]+scale[i],mp.mpf(0)),mp.mpf(1)); S.append(y)
    F=[f(s,eps) for s in S]
    for _ in range(it):
        o=sorted(range(n+1),key=lambda i:-F[i]); S=[S[i] for i in o]; F=[F[i] for i in o]
        cen=[sum(S[i][j] for i in range(n))/n for j in range(n)]
        xr=[cen[j]+(cen[j]-S[n][j]) for j in range(n)]; fr=f(xr,eps)
        if fr>F[0]:
            xe=[cen[j]+2*(cen[j]-S[n][j]) for j in range(n)]; fe=f(xe,eps)
            S[n],F[n]=(xe,fe) if fe>fr else (xr,fr)
        elif fr>F[n-1]: S[n],F[n]=xr,fr
        else:
            xc=[cen[j]+mp.mpf('0.5')*(S[n][j]-cen[j]) for j in range(n)]; fc=f(xc,eps)
            if fc>F[n]: S[n],F[n]=xc,fc
            else:
                for i in range(1,n+1):
                    S[i]=[(S[i][j]+S[0][j])/2 for j in range(n)]; F[i]=f(S[i],eps)
    o=max(range(n+1),key=lambda i:F[i]); return S[o],F[o]

random.seed(3)
print(' eps        D(eps)          D/sqrt(eps)     argmax p')
for k in (4,6,8,10,12,14):
    eps=mp.mpf(10)**(-k)
    best=None
    seeds=[[mp.mpf(1),mp.mpf('0.5'),mp.mpf('0.5'),mp.sqrt(eps)],
           [mp.mpf(1),mp.sqrt(eps),mp.mpf(1),mp.mpf(0)],
           [mp.mpf(1),mp.mpf('0.1'),mp.mpf('0.9'),mp.sqrt(eps)],
           [mp.mpf(1),mp.mpf(0),mp.mpf(0),mp.sqrt(eps)]]
    for s in seeds:
        for sc in (mp.mpf('0.2'), mp.sqrt(eps), eps**mp.mpf('0.75')):
            x,val=nelder(s,eps,200,[sc]*4)
            if best is None or val>best[1]: best=(x,val)
        s=best[0]
    x,val=best
    print(f' 1e-{k:<3d}  {mp.nstr(val,8):>14}  {mp.nstr(val/mp.sqrt(eps),8):>14}  '
          f'[{", ".join(mp.nstr(t,6) for t in x)}]')
