"""The other two 2-dimensional carriers.
{p1=0,p4=1} (the CC<->DD alternation trap, W=(1/2,0,0,1/2)):  claim the 16 rows reduce to
    (III) p3 <= 1-2u      (deviator ALLC)     (IV) p2 <= -2v   (deviator ALLD)
{p2=1,p3=0} (the other alternating sheet): limit-Nash iff T=S, and there for the whole face."""
from fractions import Fraction as F
import random
from limlib import nash_rows, DET
random.seed(31); rnd = lambda: F(random.randint(1,59),60)
G = {''.join(str(int(x)) for x in s): k for k,s in enumerate(DET)}

bad = ties = n = 0
for _ in range(30):
    p2, p3 = rnd(), rnd()
    rows, W = nash_rows((F(0), p2, p3, F(1)))
    assert W == [F(1,2), 0, 0, F(1,2)], W
    # (III) p3 + u <= (1+p3)/2  ->  master a(u-1)+bv<=c with a=1, b=0, c=(1-p3)/2 - 1 + 1 ... derive:
    # payoff gap = (p3+u)/(1+p3) - 1/2 <= 0  ->  multiply (1+p3): u + p3 - (1+p3)/2 <= 0
    #            = 1*(u-1) + 0*v - ((1-p3)/2)  <= 0
    III = (F(1), F(0), -(1+p3)/2)
    # (IV) (1+v)/(2-p2) - 1/2 <= 0 -> 1+v - (2-p2)/2 <= 0 -> 0*(u-1) + 1*v - (-p2/2) <= 0
    IV = (F(0), F(1), -p2/2)
    def ismult(r,g):
        return r != (0,0,0) and r[0]*g[1]==r[1]*g[0] and r[0]*g[2]==r[2]*g[0] and \
               ((r[0]>0)==(g[0]>0) if g[0] else (r[1]>0)==(g[1]>0))
    assert any(ismult(r,III) for r in rows), (p2,p3,III,rows)
    assert any(ismult(r,IV) for r in rows), (p2,p3,IV,rows)
    det = III[0]*IV[1]-III[1]*IV[0]
    for a,b,c in rows:
        n += 1
        if (a,b,c)==(0,0,0): ties += 1; continue
        lam=(a*IV[1]-b*IV[0])/det; mu=(III[0]*b-III[1]*a)/det; nu=lam*III[2]+mu*IV[2]-c
        if not (lam>=0 and mu>=0 and nu>=0):
            bad += 1
            if bad<4: print('  NO CERT', (p2,p3), (a,b,c), (lam,mu,nu))
print(f'face p1=0,p4=1 : {n} rows over 30 points, {ties} ties, {bad} without a certificate')
print('   -> strictly feasible with p2,p3 in (0,1) iff  u < 1/2 and v < 0')

# the (1,0) sheet
print('\nface p2=1,p3=0 (the other alternating sheet):')
for (u,v) in [(F(1),F(0)),(F(2),F(1)),(F(-1),F(-2)),(F(1,2),F(1,2)),(F(3),F(1))]:
    okw = okn = 0
    for _ in range(20):
        p1,p4 = rnd(), rnd()
        rows,W = nash_rows((p1,F(1),F(0),p4))
        assert W == [0,F(1,2),F(1,2),0]
        if all(a*(u-1)+b*v <= c for a,b,c in rows): okw += 1
        if all(a*(u-1)+b*v <  c for a,b,c in rows if (a,b,c)!=(0,0,0)): okn += 1
    print(f'   game ({u},{v})  T-S={1+v-u}: weakly Nash {okw}/20, strictly {okn}/20')
