"""Exact renewal decomposition at eps=0 around the focal offer state <DC,CD> (s*=9)."""
import sympy as sp
from exact import EX, code2q, cst, x, margin_limit
p = list(EX); p[9] = x
SSTAR = 9
names = ['CC','CD','DC','DD']
def step(s, q):
    """eps=0 transitions from state s: list of (prob, next_state, outcome)."""
    a = s // 4
    out = []
    pc = sp.sympify(p[s]); qc = q[cst(s)]
    for fD in (0,1):
        pf = pc if fD == 0 else 1 - pc
        for qD in (0,1):
            pq = qc if qD == 0 else 1 - qc
            pr = sp.simplify(pf*pq)
            if pr != 0:
                o = 2*fD + qD
                out.append((pr, 4*o + a, o))
    return out
def gain(o):  # focal net exploitations: DC +1, CD -1
    return {2: 1, 1: -1}.get(o, 0)

for c in [52620, 40324]:
    q = code2q(c)
    # hitting probabilities of s* and expected net / length until hitting s*
    hs = sp.symbols('h0:16'); Ls = sp.symbols('L0:16'); As = sp.symbols('A0:16')
    eqs = []
    for s in range(16):
        tr = step(s, q)
        def val(v, ns):  # value after entering ns: 0 if ns==s*
            return 0 if ns == SSTAR else v[ns]
        eqs.append(sp.Eq(As[s], sum(pr*(1 if ns==SSTAR else As[ns]) for pr,ns,o in tr)))
    solA = sp.solve(eqs, As, dict=True)
    print(c, ''.join(map(str,q)), "hit prob of s* from each state:", [sp.simplify(solA[0][a]) if a in solA[0] else 'free' for a in As] if solA else 'nosol')
    eqh=[]; eqL=[]
    for s in range(16):
        tr = step(s, q)
        eqh.append(sp.Eq(hs[s], sum(pr*(gain(o) + (0 if ns==SSTAR else hs[ns])) for pr,ns,o in tr)))
        eqL.append(sp.Eq(Ls[s], sum(pr*(1 + (0 if ns==SSTAR else Ls[ns])) for pr,ns,o in tr)))
    sh = sp.solve(eqh, hs, dict=True)[0]; sL = sp.solve(eqL, Ls, dict=True)[0]
    # conditional on focal's action at s*
    for fD,label in [(0,'offer (focal C)'),(1,'no offer (focal D)')]:
        qD = 1 - q[cst(SSTAR)]
        o = 2*fD + qD; ns = 4*o + SSTAR//4
        net = gain(o) + (0 if ns==SSTAR else sh[hs[ns]])
        ln = 1 + (0 if ns==SSTAR else sL[Ls[ns]])
        print(f"   {label}: co-player acts {'D' if qD else 'C'}; outcome {names[o]}; focal net (DC-CD) until next offer state = {sp.simplify(net)}; expected rounds = {sp.simplify(ln)}")
    m_ren = sp.simplify((x*(-1) + (1-x)*sp.Rational(5,8)) / (x*3 + (1-x)*sp.Rational(25,4)))
    print("   renewal-reward margin:", sp.factor(m_ren), " equals exact limit:", sp.simplify(m_ren - (5-13*x)/(2*(25-13*x)))==0)
