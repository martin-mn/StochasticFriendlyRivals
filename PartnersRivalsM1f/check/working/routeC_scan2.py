import sys, numpy as np
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib as L
eps=1e-7
P1=[0.0,0.1,0.3,0.5,0.7,0.9,0.99]
P4=[0.01,0.1,0.3,0.5,0.7,0.9,1.0]
US=[-2,-1,-0.5,0,0.25,0.5,0.75,1,1.25,1.5,2,3,5]
VS=[-1,-0.5,0,0.25,0.5,0.75,1,1.5,2,3,5,10]
for corner in [(0,1),(1,0)]:
    print("="*70); print("corner (p2,p3) =",corner)
    print("      " + "".join(f"{v:>6}" for v in VS) + "   <- v")
    for uu in US:
        line=f"u={uu:5.2f} "
        for vv in VS:
            if uu+vv<=1+1e-12: line+="     ."; continue
            # fraction of alternators that are Nash
            n=0; tot=0
            for p1 in P1:
                for p4 in P4:
                    if p1==0.0 and p4==1.0: continue
                    p=(p1,corner[0],corner[1],p4); tot+=1
                    if L.nash_gap(p,uu,vv,eps) > -1e-5: n+=1
            line+=f"{n:>3}/{tot:<2}"
        print(line)
