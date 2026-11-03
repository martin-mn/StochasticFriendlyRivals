"""Exact-Fraction verification of the Route-C closed forms.

Predicted eps->0 limit of pi(sigma,p) against an alternator p=(a,c2,c3,d), and hence the
predicted Nash gap  min_sigma [ (1+u+v)/2 - pi(sigma,p) ].  Compared with the exact-Fraction
pair chain at eps = 10^-8 and 10^-9 (ratio test separates 'limit 0' from 'limit < 0').
"""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib as L
from itertools import product

def piA(g,a,d,u,v):
    s1,s2,s3,s4=g
    E=(1+u+v)/F(2)
    key=''.join(map(str,g))
    if key in ('0010','0011','1010','1011'): return E
    if key in ('0000','0100','1000'): return d*(v+1)/(d+1)
    if key in ('1101','1110','1111'): return (u*(a-1)-1)/(a-2)
    if key=='0001': return (d-u*(d-1)+(v+1)*(a*d-d+1))/(a*d-d+3)
    if key=='0101': return (a*(v+1)-u*(d-1)+1)/(a-d+3)
    if key=='0110': return d*(u+v+2)/(-a+3*d+1)
    if key=='0111': return (a*(v+1)+u*(a*d-d+1)+1)/(a*d-d+3)
    if key=='1001': return (d-u*(a-1)-(a-1)*(v+1))/(-3*a+d+3)
    if key=='1100': return (d*(a-2)*(v+1)+(d+1)*(u*(a-1)-1))/(2*(a-2)*(d+1))
    raise KeyError(key)

def piB(g,a,d,u,v):
    E=(1+u+v)/F(2); key=''.join(map(str,g))
    if key in ('0000','0001','0100','1000'): return v+1
    if key in ('0111','1101','1110','1111'): return u
    if key=='1100': return E
    if key=='0010': return d*(v+2)/(-a+2*d+1)
    if key=='0011': return (a*d*(v+1)+d+u*(a*d-a-d+1))/(2*(a*d-a+1))
    if key=='0101': return (a*(d+1)*(v+1)+u*(a*d-a-2*d+2))/(2*(a*d-d+1))
    if key=='0110': return (2*d*(v+1)+2*d+u*(-a+d+1))/(-3*a+5*d+3)
    if key=='1001': return (2*d-2*u*(a-1)+(v+1)*(-a+d+1))/(-5*a+3*d+5)
    if key=='1010': return d*(u*(a-1)+(a-1)*(v+1)-1)/(2*a*d+a-3*d-1)
    if key=='1011': return (d-u*(a-1))/(-2*a+d+2)
    raise KeyError(key)

GEN=list(product([0,1],repeat=4))
def pred_gap(corner,a,d,u,v):
    f = piA if corner==(0,1) else piB
    E=(1+u+v)/F(2)
    return min(E-f(g,a,d,u,v) for g in GEN)

def exact_gap(corner,a,d,u,v,eps):
    p=(a,F(corner[0]),F(corner[1]),d)
    w=L.stationary_exact(p,p,eps); own=w[0]+w[1]*u+w[2]*(1+v)
    best=max((w2[0]+w2[1]*u+w2[2]*(1+v)) for w2 in
             (L.stationary_exact([F(x) for x in g],p,eps) for g in GEN))
    return own-best

AS=[F(0),F(1,7),F(1,3),F(1,2),F(3,4),F(9,10),F(99,100)]
DS=[F(1,100),F(1,7),F(1,3),F(1,2),F(3,4),F(9,10),F(1)]
GAMES=[]                              # (u,v) strictly above u+v=1, all sign patterns
for uu,vv in [(F(0),F(3,2)),(F(1,2),F(1)),(F(-1),F(5,2)),(F(-3),F(9,2)),(F(3),F(-1,2)),
              (F(5),F(-3)),(F(2),F(1)),(F(1,4),F(2)),(F(9,4),F(5,4)),(F(-1,2),F(2)),
              (F(4),F(-1)),(F(3,2),F(1,2)),(F(2),F(1)),(F(11,4),F(7,4)),(F(-2),F(4)),
              (F(6),F(-2)),(F(1),F(1)),(F(0),F(4)),(F(7,2),F(5,2))]:
    assert uu+vv>1
    GAMES.append((uu,vv))

bad=0; tested=0; maxerr=F(0)
for corner in [(0,1),(1,0)]:
    for a in AS:
        for d in DS:
            if a==0 and d==1: continue
            for (uu,vv) in GAMES:
                pg=pred_gap(corner,a,d,uu,vv)
                e8=exact_gap(corner,a,d,uu,vv,F(1,10**8))
                e9=exact_gap(corner,a,d,uu,vv,F(1,10**9))
                tested+=1
                # Richardson: limit ~ (10*e9 - e8)/9
                rich=(10*e9-e8)/9
                err=abs(rich-pg); maxerr=max(maxerr,err)
                if err>F(1,10**5):
                    bad+=1
                    if bad<8: print("MISMATCH",corner,a,d,uu,vv,"pred",float(pg),"rich",float(rich))
print(f"tested {tested} (corner,a,d,u,v) quadruples; mismatches {bad}; max |Richardson - predicted| = {float(maxerr):.3e}")

# --- the reduced criterion ---
def crit_A(a,d,u,v):
    X=u-v-1; Y=u+v-1
    return (d*X+Y+2 >= 0) and ((a-1)*X+Y >= 0)
def crit_B(a,d,u,v):
    return u-v-1 == 0
bad2=0
for corner in [(0,1),(1,0)]:
    for a in AS:
        for d in DS:
            if a==0 and d==1: continue
            for (uu,vv) in GAMES:
                pg=pred_gap(corner,a,d,uu,vv)
                c=(crit_A if corner==(0,1) else crit_B)(a,d,uu,vv)
                if (pg>=0)!=c:
                    bad2+=1
                    if bad2<8: print("CRIT MISMATCH",corner,a,d,uu,vv,float(pg),c)
print("reduced-criterion mismatches:",bad2)
