"""V2 (item 3): the ZD-equaliser family, exact.

p_tilde = (p1-1,p2-1,p3,p4) = -b (S_Y - kappa 1),  S_Y = (R,T,S,P) = (1,1+v,u,0) read in the
EQUALISER'S OWN state labelling (its action first), so
    p = ( 1-b(1-kappa), 1-b(1+v-kappa), b(kappa-u), b*kappa ).

Everything below is in exact Fractions with the clamped chain of m1lib.
"""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/Users/martin/Documents/CL1/PartnersRivalsM1f/check')
import m1lib

def zd(u,v,kap,b):
    u,v,kap,b = map(F,(u,v,kap,b))
    return (1-b*(1-kap), 1-b*(1+v-kap), b*(kap-u), b*kap)

def coop_rate(p,q,eps):
    """c_X: long-run frequency with which the RESIDENT p (focal in this chain) plays C."""
    w = m1lib.stationary_exact(p,q,eps)
    return w[0]+w[1]

def pi_opponent(p,q,eps,u,v):
    """payoff of the co-player q against the equaliser p (p focal in the chain)."""
    w = m1lib.stationary_exact(p,q,eps)
    return w[0] + w[1]*(1+F(v)) + w[2]*F(u)

def check(u,v,kap,b,epslist):
    p = zd(u,v,kap,b)
    assert all(0<=x<=1 for x in p), (p,"infeasible")
    print(f"  (u,v)=({u},{v}) kappa={kap} b={b}  p={tuple(str(x) for x in p)}")
    for eps in epslist:
        e = F(eps)
        # (i) the exact identity  pi_Y = kappa + eps(1-2c_X)/(b(1-2eps))  for all 16 co-players
        worst = 0
        for q in m1lib.DET_EXACT:
            lhs = pi_opponent(p,q,e,u,v)
            c   = coop_rate(p,q,e)
            rhs = F(kap) + e*(1-2*c)/(F(b)*(1-2*e))
            worst = max(worst, abs(lhs-rhs))
        # (ii) the exact Nash gap identity
        cs   = [coop_rate(p,q,e) for q in m1lib.DET_EXACT]
        cself= coop_rate(p,p,e)
        gap_direct = pi_opponent(p,p,e,u,v) - max(pi_opponent(p,q,e,u,v) for q in m1lib.DET_EXACT)
        gap_form   = 2*e/(F(b)*(1-2*e))*(min(cs)-cself)
        print(f"    eps={eps}: max|identity residual| = {worst};  gap {float(gap_direct):.12e}"
              f"  formula-residual {gap_direct-gap_form};  gap/eps = {float(gap_direct/e):.6f}"
              f"  min_c={float(min(cs)):.6f} c_self={float(cself):.6f}")
    return p

print("=== the worked example of D1 ===")
check(F(1,2), 5, 1, F(9,50), ['1/1000','1/100000','1/10000000'])
print("   predicted leading coefficient 2/b*(min c - c_self) = ", float(2/F(9,50)*(F(1,6)-1)), " = -250/27 =", -250/27)

print()
print("=== a spread of feasible (u,v,kappa,b) ===")
for (u,v,kap,b) in [(F(1,2),5,1,F(9,50)), (F(-1,2),F(3,2),F(1,4),F(1,5)),
                    (F(1,4),F(-1,2),F(1,2),F(1,3)), (0,3,1,F(1,7)),
                    (1,F(1,2),1,F(1,4)), (F(-2),F(-1,2),F(1,3),F(1,6))]:
    check(u,v,kap,b,['1/100000'])
