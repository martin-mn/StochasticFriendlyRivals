"""The step of SI 'The dimension of the stable set' from alpha = beta = 0 to the plane of zero-determinant equalisers,
in exact algebra: the average-reward optimality equations of the deviator's decision problem against p at eps = 0,
with both actions optimal at every state (alpha = b_CC - b_DC = 0, beta = b_CD - b_DD = 0), b_t = r_t + h_t.
States from the deviator's view; the resident cooperates at s with q_s = p_{M(s)}; r = (1, u, 1+v, 0)."""
import sympy as sp
pCC, pCD, pDC, pDD, u, v, g = sp.symbols('pCC pCD pDC pDD u v g')
h = sp.symbols('hCC hCD hDC hDD')
r = (1, u, 1 + v, 0)
b = [r[t] + h[t] for t in range(4)]
q = (pCC, pDC, pCD, pDD)                       # q_s = p_{M(s)}, M swaps CD and DC
eqs = []
for s in range(4):
    VC = q[s] * b[0] + (1 - q[s]) * b[1]       # deviator plays C: next outcome CC or CD
    VD = q[s] * b[2] + (1 - q[s]) * b[3]       # deviator plays D: next outcome DC or DD
    eqs += [sp.Eq(g + h[s], VC), sp.Eq(g + h[s], VD)]
eqs.append(sp.Eq(h[3], 0))                     # normalisation of the relative values
sol = sp.solve(eqs, list(h) + [g], dict=True)
print('solutions (generic p):', sol)
# eliminate: the system is linear in (h, g); its solvability conditions on p are what alpha = beta = 0 imposes
M, rhs = sp.linear_eq_to_matrix([e.lhs - e.rhs for e in eqs], list(h) + [g])
aug = M.row_join(rhs)
print('rank of the coefficient matrix:', M.rank(), ' rank of the augmented matrix (generic p):', aug.rank())
# on the plane E the system is solvable:
k = 1 - pCC + pDD
onE = {pCD: pCC - v * k, pDC: pDD - u * k}
print('rank on E:', M.subs(onE).rank(), aug.subs(onE).applyfunc(sp.simplify).rank())
solE = sp.solve([e.subs(onE) for e in eqs], list(h) + [g], dict=True)
print('solution on E:', [{kk: sp.simplify(vv) for kk, vv in d.items()} for d in solE])
# all solutions with p_CD and p_DC also unknown: the solvability conditions are exactly E (when 1 - p_CC + p_DD != 0)
allsol = sp.solve(eqs, list(h) + [g, pCD, pDC], dict=True)
for d in allsol:
    print('p_CD =', sp.factor(d.get(pCD, pCD)), '  p_DC =', sp.factor(d.get(pDC, pDC)), '  g =', sp.factor(d.get(g, g)))
    assert sp.simplify(d[pCD] - onE[pCD]) == 0 and sp.simplify(d[pDC] - onE[pDC]) == 0
print('alpha = beta = 0 at all four states  <=>  p on the plane E (for 1 - p_CC + p_DD != 0): confirmed')
