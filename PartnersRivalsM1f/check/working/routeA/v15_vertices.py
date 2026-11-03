import sys, numpy as np
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib
claims={'1001':lambda u,v:(u<=1e-12+1) and (v<=1+1e-12),          # WSLS
        '1111':lambda u,v:v<=1e-12,                                # ALLC
        '0000':lambda u,v:u<=1e-12}                                # ALLD
for g,cl in claims.items():
    p=np.array([float(c) for c in g]); bad=0
    for u in np.linspace(-4,4,161):
        for v in np.linspace(-4,4,161):
            if (m1lib.nash_gap(p,u,v,1e-9)>-1e-6)!=cl(u,v): bad+=1
    print(f"  {g} ({m1lib.NAMES.get(g,'')}): closed-form Nash region mismatches on 161^2 grid: {bad}")
# exact area of the 16-vertex no-Nash set on [-4,4]^2
print("  exact area share of {u>0,v>0}\\[0,1]^2 in [-4,4]^2 =", (16-1)/64)
