import mp as M, mpmath as mp
mp.mp.dps=50
# D(eps) along the predicted optimal family p=(1, p2, 1, p4), p2 = sqrt(2 s eps)
print('  s      eps      best d_eps      d/sqrt(eps)   -sqrt(s/2)')
for (u,v) in [(mp.mpf(-1)/2,mp.mpf(1)/2), (mp.mpf(-1),mp.mpf(0)), (mp.mpf(0),mp.mpf(1)/2),
              (mp.mpf(-2),mp.mpf(0))]:
    s=1-u-v
    for k in (8,12,16):
        eps=mp.mpf(10)**(-k)
        best=None
        for fp2 in [mp.mpf(t)/20 for t in range(4,40)]:
            p2=fp2*mp.sqrt(2*s*eps)
            for fp4 in ('0.01','0.05','0.1','0.2','0.3','0.32','0.4','0.6'):
                p4=mp.mpf(fp4)*p2
                val=M.d_eps([1,p2,1,p4],u,v,eps)
                if best is None or val>best[0]: best=(val,p2,p4)
        val,p2,p4=best
        print(f' {mp.nstr(s,3):>5} 1e-{k:<3d} {mp.nstr(val,8):>14} {mp.nstr(val/mp.sqrt(eps),8):>13}'
              f'  {mp.nstr(-mp.sqrt(s/2),8):>11}   p2/sqrt(2 s eps)={mp.nstr(p2/mp.sqrt(2*s*eps),4)}'
              f' p4/p2={mp.nstr(p4/p2,3)}')
