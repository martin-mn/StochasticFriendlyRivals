import sympy as sp
a,b,c,d = sp.symbols('a b c d')
M = sp.Matrix([[a**2, 2*a*(1-a), (1-a)**2],
               [b*c,  b*(1-c)+c*(1-b), (1-b)*(1-c)],
               [d**2, 2*d*(1-d), (1-d)**2]])
x,y,z = sp.symbols('x y z')
v = sp.Matrix([[x,y,z]])
eqs = list((v*M - v))[:2] + [x+y+z-1]
sol = sp.solve(eqs,[x,y,z],dict=True)[0]
X = sp.cancel(sp.together(sol[x])); Y = sp.cancel(sp.together(sol[y])); Z = sp.cancel(sp.together(sol[z]))
for nm,e in (('X',X),('Y',Y),('Z',Z)):
    n,dd = sp.fraction(sp.cancel(e))
    print(nm,'num =', sp.factor(sp.expand(n)))
    print(nm,'den =', sp.factor(sp.expand(dd)))
    print()
