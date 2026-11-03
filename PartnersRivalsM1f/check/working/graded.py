"""Graded Nash field: relative volume of the limit-Nash set inside each carrying subvariety.

facet {p1=1}:  V1 = vol{(p2,p3,p4) in [0,1]^3 : p4 v <= 1-p2,  p3 v <= (1-p2)(1-u)}
facet {p4=0}:  V4 = vol{(p1,p2,p3) in [0,1]^3 : p3 <= (1-p1)|u|, p3(1+v) <= (1-p2)|u|}, u<0
sheet {p2=0,p3=1}: A = area{(p1,p4) in [0,1]^2 : p4 X + Y + 2 >= 0, (p1-1) X + Y >= 0}
"""
import numpy as np

def V1(u, v):
    t = 1.0 - u
    if v <= 0:
        if u <= 1: return 1.0
        k = (u - 1) / (-v)
        return 1 - k/2 if k <= 1 else 1/(2*k)
    if u >= 1: return 0.0
    s = v / t
    b, c = max(v, s), min(v, s)
    if b <= 1: return 1 - b/2 - c*c/(6*b)
    if c <= 1: return (0.5 - c*c/6)/b
    return 1.0/(3*v*s)

def V1_mc(u, v, n=400000, seed=0):
    r = np.random.default_rng(seed).random((n,3))
    p2,p3,p4 = r[:,0],r[:,1],r[:,2]
    return float(((p4*v <= 1-p2) & (p3*v <= (1-p2)*(1-u))).mean())

def V4(u, v):
    if u >= 0: return 0.0
    U = -u
    A = 1.0/U
    B = max(0.0, 1+v)/U
    L = min(1.0, 1/A, (1/B if B > 0 else np.inf))
    return L - (A+B)*L*L/2 + A*B*L**3/3

def V4_mc(u, v, n=400000, seed=1):
    r = np.random.default_rng(seed).random((n,3))
    p1,p2,p3 = r[:,0],r[:,1],r[:,2]
    return float((((1-p1)*u + p3 <= 0) & ((1-p2)*u + p3*(1+v) <= 0)).mean())

def A01(u, v):
    """relative area of the Nash set inside the alternating sheet (p1,0,1,p4)."""
    X, Y = u-v-1, u+v-1
    cl = lambda z: max(0.0, min(1.0, z))
    if   X > 0: m4 = cl(1 + (Y+2)/X)
    elif X < 0: m4 = cl((Y+2)/(-X))
    else:       m4 = 1.0 if Y+2 >= 0 else 0.0
    if   X > 0: m1 = cl(Y/X)
    elif X < 0: m1 = cl(1 + Y/(-X))
    else:       m1 = 1.0 if Y >= 0 else 0.0
    return m4*m1

def A01_mc(u, v, n=400000, seed=2):
    r = np.random.default_rng(seed).random((n,2))
    p1,p4 = r[:,0],r[:,1]
    X, Y = u-v-1, u+v-1
    return float(((p4*X+Y+2 >= 0) & ((p1-1)*X+Y >= 0)).mean())

tests = [(-3,-3),(-1,-1),(-0.5,-2),(0.5,-1),(2,-1),(2,-3),(3,-0.5),(0.5,0.5),(-1,0.5),
         (0.9,0.5),(0.5,2),(-2,2),(0,3),(-3,0.5),(1.5,0.5),(2,2),(3,1),(1,-1),(0.25,0.25),
         (-0.5,0.75),(-4,-0.5),(-1,-0.5),(-0.2,-3)]
print(f'{"u":>6s}{"v":>6s} | {"V1":>9s}{"V1 mc":>9s} | {"V4":>9s}{"V4 mc":>9s} | {"A01":>9s}{"A01 mc":>9s}')
w1=w4=wa=0
for u,v in tests:
    a,b = V1(u,v), V1_mc(u,v); c,d = V4(u,v), V4_mc(u,v); e,f = A01(u,v), A01_mc(u,v)
    w1=max(w1,abs(a-b)); w4=max(w4,abs(c-d)); wa=max(wa,abs(e-f))
    print(f'{u:6.2f}{v:6.2f} | {a:9.5f}{b:9.5f} | {c:9.5f}{d:9.5f} | {e:9.5f}{f:9.5f}')
print(f'max |closed form - Monte Carlo| :  V1 {w1:.4f}   V4 {w4:.4f}   A01 {wa:.4f}')

# item 3's four branches, below/on the switch line only
def V_item3(u,v):
    t=1-u; s=v/t if t!=0 else np.inf
    if u<=1 and v<=0: return 1.0
    if u>1: return 1-(u-1)/(2*abs(v))
    if u<1 and 0<v<=1:
        a,b=min(s,v),max(s,v); return 1-b/2-a*a/(6*b)
    if u<=0 and v>1: return (0.5-s*s/6)/v
    return None
print('\nagreement with item 3 (partners.md) on u+v <= 1, where its four branches apply:')
bad=0; n=0
for u in np.linspace(-6,6,121):
    for v in np.linspace(-6,6,121):
        if u+v<=1:
            x=V_item3(u,v)
            if x is not None:
                n+=1
                if abs(x-V1(u,v))>1e-12: bad+=1
print(f'  {n} games with u+v<=1, {bad} disagreements between item 3s V and V1')
print('\nwhere item 3s branch labels break down ABOVE the line (they were only claimed below it):')
for u,v in [(0.9,0.5),(0.5,2),(0.8,0.9),(2,-1.5)]:
    x=V_item3(u,v); xs='   none' if x is None else f'{x:9.5f}'
    print(f'  u={u:5.2f} v={v:5.2f}   item3 -> {xs}   V1 -> {V1(u,v):9.5f}  mc {V1_mc(u,v):9.5f}')
