from fractions import Fraction as F
import pickle, random
import numpy as np
from limlib import nash_rows, DET

# ---- 1. targeted witness at u=3, v=-1/1000  (the coarse sampler missed it)
for (u, v, p) in [(F(3), F(-1,1000), (F(1), 1-F(1,10**5), F(1,2), F(1,2))),
                  (F(8), F(-1,10**4), (F(1), 1-F(1,10**7), F(1,2), F(1,2)))]:
    rows, _ = nash_rows(p)
    nz = [r for r in rows if r != (0,0,0)]
    print(f'witness u={u} v={v}  p={tuple(str(x) for x in p)}  strictly Nash:',
          all(a*(u-1)+b*v < c for a,b,c in nz))

# ---- 2. the alternating sheet (p1,0,1,p4): do C1,C2 exhaust the 16 at EVERY game?
#   C1: p4 X + Y + 2 >= 0,  C2: (p1-1) X + Y >= 0,  X=u-v-1, Y=u+v-1
#   in master form a(u-1)+b v <= c :  -C1 -> a=-p4,b=p4-1,c=2 ... derive from rows instead
print('\nalternating sheet (p1,0,1,p4):')
random.seed(4); rnd = lambda: F(random.randint(1,59),60)
bad = 0; tot = 0
for _ in range(25):
    p1, p4 = rnd(), rnd()
    rows, Wp = nash_rows((p1, F(0), F(1), p4))
    assert Wp == [0, F(1,2), F(1,2), 0], Wp
    # C1 as a(u-1)+bv<=c :  -(p4X+Y+2)<=0 -> -(p4(u-v-1)+(u+v-1)+2) <=0
    #   = -(p4+1)u + (p4-1)v + p4 -1 <= 0 -> (-(p4+1))(u-1) + (p4-1)v <= 1-p4+p4+1... expand:
    C1 = (-(p4+1), p4-1, F(2))          # -(p4 X + Y + 2)
    C2 = (-p1, p1-2, F(0))              # -((p1-1) X + Y)
    for _ in range(5):
        uu, vv = F(random.randint(-40,40),7), F(random.randint(-40,40),7)
        X, Y = uu-vv-1, uu+vv-1
        assert C1[0]*(uu-1)+C1[1]*vv-C1[2] == -(p4*X+Y+2)
        assert C2[0]*(uu-1)+C2[1]*vv-C2[2] == -((p1-1)*X+Y)
    def ismult(r, g):
        return r != (0,0,0) and r[0]*g[1]==r[1]*g[0] and r[0]*g[2]==r[2]*g[0] and \
               ((r[0]>0)==(g[0]>0) if g[0] else (r[1]>0)==(g[1]>0))
    assert any(ismult(r, C1) for r in rows), ('C1 not a row', p1, p4, C1, rows)
    assert any(ismult(r, C2) for r in rows), ('C2 not a row', p1, p4, C2, rows)
    for a,b,c in rows:
        tot += 1
        if (a,b,c) == (0,0,0): continue
        det = C1[0]*C2[1]-C1[1]*C2[0]
        lam = (a*C2[1]-b*C2[0])/det
        mu  = (C1[0]*b-C1[1]*a)/det
        nu  = lam*C1[2]+mu*C2[2]-c
        if not (lam>=0 and mu>=0 and nu>=0):
            bad += 1
            if bad < 4: print('  NO CERT', (p1,p4), (a,b,c), (lam,mu,nu))
print(f'  {tot} rows over 25 random (p1,p4), {bad} without a Farkas certificate against C1,C2')

# ---- 3. the sixteen over the plane
VROWS = pickle.load(open('vrows.pkl','rb'))
A = np.array([[float(a) for a,b,c in VROWS[g]] for g in VROWS])
B = np.array([[float(b) for a,b,c in VROWS[g]] for g in VROWS])
C = np.array([[float(c) for a,b,c in VROWS[g]] for g in VROWS])
NZ = ~((A==0)&(B==0)&(C==0))
N = 1601; L = 8.0
uu = np.linspace(-L, L, N); vv = np.linspace(-L, L, N)
U, V = np.meshgrid(uu, vv, indexing='ij')
cnt = np.zeros_like(U, dtype=int)
for i in range(16):
    ok = np.ones_like(U, dtype=bool)
    for j in range(16):
        if NZ[i,j]:
            ok &= (A[i,j]*(U-1) + B[i,j]*V <= C[i,j] + 1e-12)
    cnt += ok
print('\nthe sixteen, weak count on a 1601^2 grid over [-8,8]^2:')
for k in range(9):
    m = (cnt==k)
    if m.any(): print(f'  N={k}: {100*m.mean():6.3f} % of the square')
print('  N=8 set == {u<=0 and v<=0} ?', np.array_equal(cnt==8, (U<=1e-12)&(V<=1e-12)))
print('  N=0 set == {u>0 and v>0 and not(u<=1 and v<=1)} ?',
      np.array_equal(cnt==0, (U>1e-12)&(V>1e-12)&~((U<=1+1e-12)&(V<=1+1e-12))))
# continuum dimension field
dim3 = (U < 1-1e-12) | (V < -1e-12)
print(f'  continuum dim 3 on {100*dim3.mean():6.3f} %, dim 2 on {100*(~dim3).mean():6.3f} % of the square')
for k in range(9):
    m = cnt==k
    if m.any():
        print(f'    of the games with N={k}: {100*(dim3&m).sum()/m.sum():6.2f} % have dim 3')
