import sys, random
from fractions import Fraction as F
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib as L
from routeC_verify import GEN, pred_gap, exact_gap
random.seed(11)

print("[A] on T=S (u-v=1) the payoff difference vanishes exactly at every eps: alternators are rivals")
for corner in [(0,1),(1,0)]:
    for (a,d) in [(F(3,7),F(2,5)),(F(0),F(1,3)),(F(9,10),F(1))]:
        p=(a,F(corner[0]),F(corner[1]),d)
        m=max(abs(L.stationary_exact(p,q,F(1,1000))[2]-L.stationary_exact(p,q,F(1,1000))[1])
              for q in L.DET_EXACT)
        # T=S means the coefficient (T-S) multiplying (w_DC-w_CD) is zero, so difference is 0
        print(f"   corner={corner} (a,d)=({a},{d}): max|w_DC-w_CD| over 16 co-players = {float(m):.4f}"
              f"  -> pi diff = (T-S)*that = 0 on T=S")

print()
print("[B] boundary of the partner region is attained (weak inequality), exact Richardson")
cases=[((0,1),F(1,2),F(5,17),F(-3),F(9,2),'u<0: d = d* = (1+u+v)/(1+v-u) = 5/17'),
       ((0,1),F(1,2),F(6,17),F(-3),F(9,2),'u<0: d just above d*'),
       ((0,1),F(6,7),F(1,2),F(5),F(-3),'v<0: a = a* = -2v/(u-v-1) = 6/7'),
       ((0,1),F(5,7),F(1,2),F(5),F(-3),'v<0: a just below a*')]
for corner,a,d,u,v,why in cases:
    pg=pred_gap(corner,a,d,u,v)
    e8=exact_gap(corner,a,d,u,v,F(1,10**8)); e9=exact_gap(corner,a,d,u,v,F(1,10**9))
    print(f"   (a,d)=({a},{d}) at (u,v)=({u},{v}): predicted gap={pg} ({float(pg):+.6f}), "
          f"Richardson={float((10*e9-e8)/9):+.9f}   [{why}]")

print()
print("[C] no efficient non-alternator above the line (spot check, 4000 random p, (u,v)=(0.5,5))")
worst=-9
for _ in range(4000):
    p=[random.random() for _ in range(4)]
    if L.in_eff_alt(p,1e-12): continue
    worst=max(worst,L.efficiency_deficit(p,0.5,5.0,1e-8))
print("   max efficiency deficit over non-alternators =",f"{worst:.6f}","(0 would mean efficient)")

print()
print("[D] partner witnesses at extreme games above the line, exact Richardson gap")
for (u,v,a,d) in [(F(-10),F(12),F(1,2),F(1,10)),(F(20),F(-15),F(49,50),F(1,2)),
                  (F(0),F(1001,1000),F(0),F(1,1000)),(F(1001,1000),F(0),F(99,100),F(1))]:
    for corner in [(0,1)]:
        e8=exact_gap(corner,a,d,u,v,F(1,10**8)); e9=exact_gap(corner,a,d,u,v,F(1,10**9))
        print(f"   p=({a},0,1,{d}) at (u,v)=({u},{v}): pred={float(pred_gap(corner,a,d,u,v)):+.6f}"
              f"  Richardson={float((10*e9-e8)/9):+.3e}")
