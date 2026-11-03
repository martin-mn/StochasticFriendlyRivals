import sympy as sp
from itertools import product
a,d,u,v = sp.symbols('a d u v')
E = (1+u+v)/sp.Integer(2)

# corner (0,1)  --- limits copied from routeC_alt_nash.py output, re-derived there
PIA = {
 '0000': d*(v+1)/(d+1),
 '0001': (d - u*(d-1) + (v+1)*(a*d-d+1))/(a*d-d+3),
 '0010': E, '0011': E,
 '0100': d*(v+1)/(d+1),
 '0101': (a*(v+1) - u*(d-1) + 1)/(a-d+3),
 '0110': d*(u+v+2)/(-a+3*d+1),
 '0111': (a*(v+1) + u*(a*d-d+1) + 1)/(a*d-d+3),
 '1000': d*(v+1)/(d+1),
 '1001': (d - u*(a-1) - (a-1)*(v+1))/(-3*a+d+3),
 '1010': E, '1011': E,
 '1100': (d*(a-2)*(v+1) + (d+1)*(u*(a-1)-1))/(2*(a-2)*(d+1)),
 '1101': (u*(a-1)-1)/(a-2), '1110': (u*(a-1)-1)/(a-2), '1111': (u*(a-1)-1)/(a-2),
}
PIB = {
 '0000': v+1, '0001': v+1,
 '0010': d*(v+2)/(-a+2*d+1),
 '0011': (a*d*(v+1) + d + u*(a*d-a-d+1))/(2*(a*d-a+1)),
 '0100': v+1,
 '0101': (a*(d+1)*(v+1) + u*(a*d-a-2*d+2))/(2*(a*d-d+1)),
 '0110': (2*d*(v+1) + 2*d + u*(-a+d+1))/(-3*a+5*d+3),
 '0111': u,
 '1000': v+1,
 '1001': (2*d - 2*u*(a-1) + (v+1)*(-a+d+1))/(-5*a+3*d+5),
 '1010': d*(u*(a-1)+(a-1)*(v+1)-1)/(2*a*d+a-3*d-1),
 '1011': (d - u*(a-1))/(-2*a+d+2),
 '1100': E, '1101': u, '1110': u, '1111': u,
}

print("### denominators positive on 0<=a<1, 0<d<=1 ?")
for nm,D in [('A0001',a*d-d+3),('A0101',a-d+3),('A0110',-a+3*d+1),('A0111',a*d-d+3),
             ('A1001',-3*a+d+3),('A1100',2*(a-2)*(d+1)),('A1101',a-2),
             ('B0010',-a+2*d+1),('B0011',2*(a*d-a+1)),('B0101',2*(a*d-d+1)),
             ('B0110',-3*a+5*d+3),('B1001',-5*a+3*d+5),('B1010',2*a*d+a-3*d-1),('B1011',-2*a+d+2)]:
    mn = sp.minimum(sp.minimum(D,a,sp.Interval(0,1)),d,sp.Interval(0,1))
    mx = sp.maximum(sp.maximum(D,a,sp.Interval(0,1)),d,sp.Interval(0,1))
    print(f"  {nm}: range over [0,1]^2 = [{mn}, {mx}]")

print()
print("### deficits (E - pi), numerator after clearing a POSITIVE denominator")
for tag,PI in [('A',PIA),('B',PIB)]:
    print('--- corner', 'A=(0,1)' if tag=='A' else 'B=(1,0)')
    for g in sorted(PI):
        defc = sp.simplify(sp.together(E - PI[g]))
        n,dd = sp.fraction(sp.cancel(defc))
        print(f"  {g}: num = {sp.factor(sp.expand(n))}   den = {sp.factor(dd)}")
