"""Reduce the 16 Nash inequalities against an alternator to their irredundant core."""
import sympy as sp
a,d,u,v,X,Y = sp.symbols('a d u v X Y')
exec(open('/Users/martin/Documents/CL1/PartnersRivalsM1f/check/routeC_conditions.py'
          ).read().split('print("### denominators')[0])
# X = u-v-1 = S-T ;  Y = u+v-1 = S+T-2R.   u=(X+Y+2)/2, v=(Y-X)/2
sub = {u:(X+Y+2)/2, v:(Y-X)/2}
DENSIGN = {  # sign of the denominator on 0<=a<1, 0<d<=1 (verified from corner values)
 'A':{'0000':+1,'0001':+1,'0010':+1,'0011':+1,'0100':+1,'0101':+1,'0110':-1,'0111':+1,
      '1000':+1,'1001':-1,'1010':+1,'1011':+1,'1100':-1,'1101':-1,'1110':-1,'1111':-1},
 'B':{'0000':+1,'0001':+1,'0010':-1,'0011':+1,'0100':+1,'0101':+1,'0110':-1,'0111':+1,
      '1000':+1,'1001':-1,'1010':-1,'1011':-1,'1100':+1,'1101':+1,'1110':+1,'1111':+1}}
out={}
for tag,PI in [('A',PIA),('B',PIB)]:
    print('='*66); print('corner', 'A = (p2,p3)=(0,1)' if tag=='A' else 'B = (p2,p3)=(1,0)')
    for g in sorted(PI):
        defc=sp.cancel(sp.together(E-PI[g])); n,dd=sp.fraction(defc)
        n = sp.expand(sp.simplify(n.subs(sub)))
        # requirement: deficit >= 0  <=>  DENSIGN * n >= 0
        req = sp.expand(DENSIGN[tag][g]*n)
        req = sp.collect(req,[X,Y])
        out[(tag,g)]=req
        print(f"  sigma {g}:  need  {req} >= 0")
print()
print('### corner A: the two candidate generators')
C1 = out[('A','0000')]; C2 = out[('A','1101')]
print(' C1 (ALLD-type 0000/0100/1000):', sp.expand(C1), '>= 0')
print(' C2 (ALLC-type 1101/1110/1111):', sp.expand(C2), '>= 0')
for g in ['0001','0101','0110','0111','1001','1100']:
    r = out[('A',g)]
    # try r = alpha*C1 + beta*C2 with alpha,beta rational functions of a,d
    al,be = sp.symbols('al be')
    sol = sp.solve(sp.Poly(sp.expand(r-al*C1-be*C2),X,Y).coeffs(),[al,be],dict=True)
    print(f"  {g}: r = {sp.expand(r)}")
    if sol:
        s=sol[0]; print(f"        = ({sp.simplify(s[al])})*C1 + ({sp.simplify(s[be])})*C2")
    else:
        print("        NOT a combination of C1,C2")
