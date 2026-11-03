import sys, numpy as np
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib as L

eps = 1e-7
worst = {}
rows=[]
for corner in [(0,1),(1,0)]:
    for p1 in [0.0,0.05,0.2,0.4,0.5,0.6,0.8,0.95,0.999]:
        for p4 in [0.001,0.05,0.2,0.4,0.5,0.6,0.8,0.95,1.0]:
            if p1==0.0 and p4==1.0: continue
            p=(p1,corner[0],corner[1],p4)
            assert L.in_eff_alt(p), p
            for u in [-1.0,-0.5,0.0,0.3,0.6,0.9,1.0,1.2,2.0,5.0]:
                for v in [-0.5,0.0,0.2,0.5,1.0,2.0,5.0,10.0]:
                    if u+v<=1.0+1e-12: continue
                    g=L.nash_gap(p,u,v,eps)
                    ed=L.efficiency_deficit(p,u,v,eps)
                    rows.append((corner,p1,p4,u,v,g,ed))
rows.sort(key=lambda r:r[5])
print("worst 25 nash gaps (most negative = worst violation):")
for r in rows[:25]:
    print(f"  corner={r[0]} p1={r[1]:6.3f} p4={r[2]:6.3f} u={r[3]:5.2f} v={r[4]:5.2f} gap={r[5]:+.6e} effdef={r[6]:+.2e}")
print("n rows", len(rows))
print("max |effdef| =", max(abs(r[6]) for r in rows))
print("count gap < -1e-5:", sum(1 for r in rows if r[5] < -1e-5))
print("count gap < -1e-9:", sum(1 for r in rows if r[5] < -1e-9))
