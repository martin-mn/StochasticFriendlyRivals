"""Route 2 -- the DISCOUNTED Akin / Press-Dyson identity, checked in exact rationals.

Conventions are the project's (see ../m1lib.py): p = (p1,p2,p3,p4) = P(C | CC,CD,DC,DD)
with the FOCAL action first, states CC=0 CD=1 DC=2 DD=3, base game (R,S,T,P)=(1,u,1+v,0).

Discounted objects
------------------
M          row-stochastic 4x4 pair-chain matrix (same M as undiscounted)
x = v0     distribution over the four states in ROUND 0
v          = (1-delta) sum_{t>=0} delta^t x M^t = (1-delta) x (I - delta M)^{-1}
           the discounted occupancy vector; sum v = 1; pi_delta = v . (R,S,T,P)
c  = (1,1,0,0)   indicator that the FOCAL player cooperated in the current state
c' = (1,0,1,0)   indicator that the CO-player cooperated
ptilde = p - c

THE IDENTITIES CHECKED HERE
  (A1)  v . ptilde = ((1-delta)/delta) ( v.c - x.c )         [Martin's conjecture]
  (A2)  v . (c - delta p) = (1-delta) (x.c)                  [equivalent, delta=0 safe]
  (A3)  v . [ delta p - c + (1-delta)(x.c) 1 ] = 0           [LaPorte's 2nd column]
  (A4)  v . [ delta q o SWAP - c' + (1-delta)(x.c') 1 ] = 0  [LaPorte's 3rd column]
  (A5)  with p1 = 1 and x.c = 1:
          delta (1-p2) v_CD = (1-delta+delta p3) v_DC + (1-delta+delta p4) v_DD
  (A6)  the discounted Press-Dyson determinant of LaPorte_notes/Press_Dyson_Discounting_2.pdf
        equals v . (R,S,T,P).

Everything is Fraction arithmetic; residuals are reported as exact zero / non-zero.
"""

import random
from fractions import Fraction as F
from itertools import product

SWAP = (0, 2, 1, 3)
C_FOCAL = (F(1), F(1), F(0), F(0))   # c
C_CO    = (F(1), F(0), F(1), F(0))   # c'


# --------------------------------------------------------------------- exact linear algebra

def solve(A, b):
    """Exact Gaussian elimination with partial (first non-zero) pivoting; A is n x n."""
    n = len(b)
    A = [row[:] for row in A]
    b = b[:]
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None:
            raise ZeroDivisionError("singular")
        A[col], A[piv] = A[piv], A[col]
        b[col], b[piv] = b[piv], b[col]
        for r in range(col + 1, n):
            f = A[r][col] / A[col][col]
            if f:
                for k in range(col, n):
                    A[r][k] -= f * A[col][k]
                b[r] -= f * b[col]
    x = [F(0)] * n
    for r in reversed(range(n)):
        x[r] = (b[r] - sum(A[r][k] * x[k] for k in range(r + 1, n))) / A[r][r]
    return x


def det4(A):
    """Exact 4x4 determinant by cofactor expansion (small, so this is fine)."""
    from itertools import permutations
    tot = F(0)
    for perm in permutations(range(4)):
        sgn = 1
        pl = list(perm)
        for i in range(4):
            for j in range(i + 1, 4):
                if pl[i] > pl[j]:
                    sgn = -sgn
        term = F(1)
        for i in range(4):
            term *= A[i][perm[i]]
        tot += sgn * term
    return tot


# --------------------------------------------------------------------- the chain

def clamp(p, eps):
    eps = F(eps)
    return [eps + (1 - 2 * eps) * F(x) for x in p]


def transition(p, q, eps=0):
    """Row-stochastic pair chain, focal p vs co-player q, with convention-A error eps."""
    p, q = clamp(p, eps), clamp(q, eps)
    M = []
    for s in range(4):
        a, b = p[s], q[SWAP[s]]
        M.append([a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)])
    return M


def occupancy(p, q, delta, x, eps=0):
    """v = (1-delta) x (I - delta M)^{-1}, exactly.  Needs NO ergodicity: I - delta M is
    invertible for every 0 <= delta < 1 and every row-stochastic M."""
    M = transition(p, q, eps)
    delta = F(delta)
    # v (I - delta M) = (1-delta) x   <=>   (I - delta M)^T v^T = (1-delta) x^T
    A = [[(F(1) if i == j else F(0)) - delta * M[j][i] for j in range(4)] for i in range(4)]
    b = [(1 - delta) * F(xi) for xi in x]
    return solve(A, b)


def x_from_openings(p0, q0):
    """Round-0 state distribution from independent opening cooperation probabilities."""
    p0, q0 = F(p0), F(q0)
    return [p0 * q0, p0 * (1 - q0), (1 - p0) * q0, (1 - p0) * (1 - q0)]


def dot(a, b):
    return sum(ai * bi for ai, bi in zip(a, b))


# --------------------------------------------------------------------- the residuals

def residuals(p, q, delta, x, eps=0):
    delta = F(delta)
    v = occupancy(p, q, delta, x, eps)
    pc = clamp(p, eps)
    qc = clamp(q, eps)
    ptilde = [pc[s] - C_FOCAL[s] for s in range(4)]

    out = {}
    out["sum_v_minus_1"] = sum(v) - 1
    if delta != 0:
        out["A1"] = dot(v, ptilde) - (1 - delta) / delta * (dot(v, C_FOCAL) - dot(x, C_FOCAL))
    out["A2"] = dot(v, [C_FOCAL[s] - delta * pc[s] for s in range(4)]) \
        - (1 - delta) * dot(x, C_FOCAL)
    colP = [delta * pc[s] - C_FOCAL[s] + (1 - delta) * dot(x, C_FOCAL) for s in range(4)]
    out["A3"] = dot(v, colP)
    colQ = [delta * qc[SWAP[s]] - C_CO[s] + (1 - delta) * dot(x, C_CO) for s in range(4)]
    out["A4"] = dot(v, colQ)
    return v, out


def reduced_A5(p, q, delta, x, eps=0):
    """Residual of (A5); only meaningful when p1 = 1 (after clamping) and x.c = 1."""
    delta = F(delta)
    v = occupancy(p, q, delta, x, eps)
    pc = clamp(p, eps)
    lhs = delta * (1 - pc[1]) * v[1]
    rhs = (1 - delta + delta * pc[2]) * v[2] + (1 - delta + delta * pc[3]) * v[3]
    return lhs - rhs


# --------------------------------------------------------------------- Press-Dyson det (A6)

def pd_determinant_payoff(p, q, delta, x, u, vv, eps=0):
    """LaPorte_notes/Press_Dyson_Discounting_2.pdf, transcribed.  Rows are CC,CD,DC,DD;
    x = (x1,x2,x3,x4) is the round-0 state distribution; note x1+x2 = x.c and
    x1+x3 = x.c'.  Returns (numerator det, denominator det, ratio)."""
    d = F(delta)
    p1, p2, p3, p4 = clamp(p, eps)
    q1, q2, q3, q4 = clamp(q, eps)
    x1, x2, x3, x4 = [F(t) for t in x]
    R, S, T, P = F(1), F(u), 1 + F(vv), F(0)
    A = x1 + x2          # x.c   (focal cooperated in round 0)
    B = x1 + x3          # x.c'  (co-player cooperated in round 0)
    col1 = [-1 + (1 - d) * x1 + d * p1 * q1,
            (1 - d) * x1 + d * p2 * q3,
            (1 - d) * x1 + d * p3 * q2,
            (1 - d) * x1 + d * p4 * q4]
    col2 = [-1 + (1 - d) * A + d * p1,
            -1 + (1 - d) * A + d * p2,
            (1 - d) * A + d * p3,
            (1 - d) * A + d * p4]
    col3 = [-1 + (1 - d) * B + d * q1,
            (1 - d) * B + d * q3,
            -1 + (1 - d) * B + d * q2,
            (1 - d) * B + d * q4]
    num = det4([[col1[i], col2[i], col3[i], [R, S, T, P][i]] for i in range(4)])
    den = det4([[col1[i], col2[i], col3[i], F(1)] for i in range(4)])
    return num, den, (num / den if den != 0 else None)


# --------------------------------------------------------------------- random rationals

def rq(rng, maxden=9):
    den = rng.randint(1, maxden)
    return F(rng.randint(0, den), den)


def rand_case(rng, eps_zero=True):
    p = [rq(rng) for _ in range(4)]
    q = [rq(rng) for _ in range(4)]
    delta = F(rng.randint(1, 19), 20)          # 0.05 .. 0.95, never 0 or 1
    if rng.random() < 0.5:
        x = x_from_openings(rq(rng), rq(rng))  # independent openings
    else:                                      # arbitrary (even correlated) round-0 law
        raw = [F(rng.randint(0, 7)) for _ in range(4)]
        if sum(raw) == 0:
            raw = [F(1)] * 4
        x = [t / sum(raw) for t in raw]
    eps = F(0) if eps_zero else F(rng.randint(1, 9), 100)
    return p, q, delta, x, eps


# --------------------------------------------------------------------- main checks

def main():
    rng = random.Random(20260912)

    print("=" * 78)
    print("(c)  EXACT-RATIONAL CHECK OF THE DISCOUNTED AKIN IDENTITY")
    print("=" * 78)
    print("  A1 : v.ptilde - ((1-d)/d)(v.c - x.c)")
    print("  A2 : v.(c - d p) - (1-d) x.c")
    print("  A3 : v.[d p - c + (1-d)(x.c) 1]      (LaPorte column 2)")
    print("  A4 : v.[d q@SWAP - c' + (1-d)(x.c') 1] (LaPorte column 3)")
    print()
    hdr = f"{'#':>3}  {'delta':>6}  {'eps':>5}  {'sum v-1':>9} {'A1':>9} {'A2':>9} {'A3':>9} {'A4':>9}"
    print(hdr)
    bad = 0
    ncase = 40
    for k in range(ncase):
        eps_zero = (k % 2 == 0)
        p, q, delta, x, eps = rand_case(rng, eps_zero)
        v, res = residuals(p, q, delta, x, eps)
        flags = [res["sum_v_minus_1"], res["A1"], res["A2"], res["A3"], res["A4"]]
        if any(f != 0 for f in flags):
            bad += 1
        print(f"{k:>3}  {float(delta):>6.2f}  {float(eps):>5.2f}  " +
              " ".join(f"{('0' if f == 0 else str(f)):>9}" for f in flags))
    print(f"\n  {ncase} cases, {bad} with a non-zero residual.")
    print("  (eps = 0 on even rows, eps > 0 on odd rows; half the x are correlated,")
    print("   i.e. NOT of the product form -- the identity does not care.)")

    # ---- also the sign/constant question, stated the other way round -------------
    print()
    print("  Sanity on the CONSTANT and the SIGN of A1, on one worked case:")
    p, q, delta, x = [F(1), F(1, 2), F(3, 4), F(0)], [F(0)] * 4, F(4, 5), x_from_openings(1, 0)
    v, res = residuals(p, q, delta, x)
    print(f"    p={p} q=ALLD delta={delta} x=CD")
    print(f"    v            = {v}")
    print(f"    v.ptilde     = {dot(v, [p[s]-C_FOCAL[s] for s in range(4)])}")
    print(f"    v.c          = {dot(v, C_FOCAL)}   x.c = {dot(x, C_FOCAL)}")
    print(f"    (1-d)/d (v.c - x.c) = "
          f"{(1-delta)/delta*(dot(v, C_FOCAL)-dot(x, C_FOCAL))}")
    print(f"    => A1 residual {res['A1']}  "
          f"(a MINUS sign or a d/(1-d) would give a non-zero here)")

    # ---- wrong-constant controls, to show the check has teeth -------------------
    print()
    print("  Controls (deliberately wrong constants) on 10 random cases -- these MUST fail:")
    rng2 = random.Random(7)
    wrongA = wrongB = 0
    for _ in range(10):
        p, q, delta, x, eps = rand_case(rng2)
        v = occupancy(p, q, delta, x, eps)
        pt = [p[s] - C_FOCAL[s] for s in range(4)]
        if dot(v, pt) + (1 - delta) / delta * (dot(v, C_FOCAL) - dot(x, C_FOCAL)) != 0:
            wrongA += 1                      # wrong SIGN
        if dot(v, pt) - delta / (1 - delta) * (dot(v, C_FOCAL) - dot(x, C_FOCAL)) != 0:
            wrongB += 1                      # reciprocal CONSTANT
    print(f"    sign flipped      : {wrongA}/10 non-zero  (10/10 expected)")
    print(f"    constant inverted : {wrongB}/10 non-zero  (10/10 expected)")

    # ---- A5, the reduced form used in part (b) ----------------------------------
    print()
    print("  A5 (p1 = 1, focal opens C):  d(1-p2) v_CD - (1-d+d p3) v_DC - (1-d+d p4) v_DD")
    rng3 = random.Random(99)
    badA5 = 0
    for _ in range(30):
        p = [F(1)] + [rq(rng3) for _ in range(3)]
        q = [rq(rng3) for _ in range(4)]
        delta = F(rng3.randint(1, 19), 20)
        x = x_from_openings(1, rq(rng3))
        if reduced_A5(p, q, delta, x) != 0:
            badA5 += 1
    print(f"    30 random cases, {badA5} non-zero residuals.")

    # ---- A6, LaPorte's discounted determinant ------------------------------------
    print()
    print("=" * 78)
    print("(d)  LaPorte discounted Press-Dyson determinant vs v.(R,S,T,P)")
    print("=" * 78)
    rng4 = random.Random(4242)
    badPD = 0
    for k in range(20):
        p, q, delta, x, eps = rand_case(rng4, eps_zero=(k % 2 == 0))
        u = F(rng4.randint(-20, 20), 10)
        vv = F(rng4.randint(-20, 20), 10)
        v = occupancy(p, q, delta, x, eps)
        pi_direct = dot(v, [F(1), u, 1 + vv, F(0)])
        num, den, pi_det = pd_determinant_payoff(p, q, delta, x, u, vv, eps)
        ok = (pi_det is not None and pi_det == pi_direct)
        if not ok:
            badPD += 1
            print(f"   MISMATCH k={k}: det={pi_det} direct={pi_direct}")
    print(f"    20 random cases, {badPD} mismatches between the determinant formula")
    print("    of Press_Dyson_Discounting_2.pdf and the direct occupancy payoff.")

    # ---------------------------------------------------------------- part (b) numerics
    print()
    print("=" * 78)
    print("(b)  UNBEATABILITY UNDER DISCOUNTING")
    print("=" * 78)
    print("Claim R1.  For T < S (wedge S / wedge E side, u - v > 1) the focal p is a")
    print("           discounted rival for EVERY q and every q-opening iff")
    print("               p0 = 1,  p1 = 1,  and  delta (1 - p2 - p3) <= 1 - delta.")
    print("           delta = 1 recovers the project's p2 + p3 >= 1 (and the excluded")
    print("           edge p2=1,p3=0 stops being excluded once delta < 1).")
    print()
    DETV = [tuple(F(b) for b in g) for g in product([0, 1], repeat=4)]

    def worst_margin(p, delta, p0):
        """min over deterministic q and q-opening of (v_CD - v_DC).  A discounted MDP for
        the co-player, so deterministic (q, q0) attain the optimum."""
        best = None
        for q in DETV:
            for q0 in (F(0), F(1)):
                v = occupancy(list(p), list(q), delta, x_from_openings(p0, q0))
                m = v[1] - v[2]
                if best is None or m < best:
                    best = m
        return best

    print(f"  {'delta':>6}  {'p':>22}  {'p0':>3}  {'pred':>5}  {'min(v_CD - v_DC)':>20}  ok")
    rows = 0
    disagree = 0
    rng5 = random.Random(31337)
    cases = []
    for delta in (F(1, 4), F(1, 2), F(3, 4), F(9, 10), F(19, 20)):
        # the interesting deterministic corners plus random p
        ps = [(F(1), F(0), F(0), F(0)), (F(1), F(1), F(0), F(0)), (F(1), F(0), F(1), F(0)),
              (F(1), F(1), F(1), F(1)), (F(1), F(0), F(0), F(1)), (F(1), F(1), F(0), F(1)),
              (F(1), F(1, 2), F(1, 2), F(0)), (F(1), F(1, 4), F(1, 4), F(1, 2)),
              (F(1), F(0), F(1, 3), F(1)), (F(1), F(1, 3), F(0), F(0))]
        ps += [(F(1), rq(rng5), rq(rng5), rq(rng5)) for _ in range(6)]
        for p in ps:
            cases.append((delta, p))
    for delta, p in cases:
        pred = (delta * (1 - p[1] - p[2]) <= 1 - delta)
        m = worst_margin(p, delta, F(1))
        ok = (pred == (m >= 0))
        rows += 1
        if not ok:
            disagree += 1
        if not ok or rows % 8 == 1:
            print(f"  {float(delta):>6.2f}  {str([str(t) for t in p]):>22}  "
                  f"{'1':>3}  {str(pred):>5}  {str(m):>20}  {'OK' if ok else 'MISMATCH'}")
    print(f"\n  {rows} (delta, p) pairs with p1 = 1 and p0 = 1: {disagree} disagreements")
    print("  between the predicate  delta(1-p2-p3) <= 1-delta  and the brute-force minimum")
    print("  of v_CD - v_DC over all 16 deterministic q and both q-openings.")

    print()
    print("Claim R2.  p0 = 1 and p1 = 1 are NECESSARY: against q = ALLC the focal is in")
    print("           CD never and in DC with positive discounted weight as soon as")
    print("           p0 < 1 or p1 < 1.")
    for p0, p in ((F(1, 2), (F(1), F(1), F(1), F(1))), (F(1), (F(9, 10), F(1), F(1), F(1)))):
        v = occupancy(list(p), [F(1)] * 4, F(9, 10), x_from_openings(p0, F(1)))
        print(f"    p0={p0} p={[str(t) for t in p]} vs ALLC(C): "
              f"v_CD={v[1]}  v_DC={v[2]}  margin={v[1]-v[2]}")

    print()
    print("Claim R3.  T > S (every Prisoner's Dilemma): efficiency forces the focal to open")
    print("           C, rivalry forces it to open D.  Shown here for the rivalry half:")
    print("           against q = ALLD the focal is never in DC, so it must never be in CD.")
    for p0 in (F(1), F(1, 2), F(0)):
        v = occupancy([F(1), F(1), F(1), F(1)], [F(0)] * 4, F(9, 10), x_from_openings(p0, F(0)))
        print(f"    p0={p0} p=ALLC vs ALLD(D): v_CD={v[1]} v_DC={v[2]} "
              f"-> v_DC - v_CD = {v[2]-v[1]}")

    print()
    print("Claim R4.  No memory-one strategy can ENFORCE pi_x - pi_y = const when delta < 1.")
    print("           delta p = c + (S_x - S_y) + gamma 1 needs delta p1 = 1 + gamma <= delta")
    print("           (so gamma <= delta - 1 < 0) and delta p4 = gamma >= 0.  Contradiction.")
    print("           At delta = 1 both hold with gamma = 0: the generalised-TFT line")
    print("           p = (1, 1-(T-S), T-S, 0).  Numerical illustration -- the undiscounted")
    print("           equaliser stops equalising the moment delta < 1:")
    u, vv = F(3, 2), F(0)                 # T - S = 1 + v - u = -1/2  (wedge S, T < S)
    TmS = 1 + vv - u
    peq = [F(1), 1 - TmS, TmS, F(0)]      # = (1, 3/2, -1/2, 0) infeasible; use a feasible one
    u, vv = F(6, 5), F(0)
    TmS = 1 + vv - u                      # = -1/5
    peq = [F(1), 1 - TmS, TmS, F(0)]
    print(f"    (u,v)=({u},{vv}) -> T-S={TmS}; equaliser p={[str(t) for t in peq]} "
          f"is outside the cube (p3<0), as T<S forces.  Use T>S instead:")
    u, vv = F(1, 2), F(0)
    TmS = 1 + vv - u                      # = 1/2
    peq = [F(1), 1 - TmS, TmS, F(0)]
    for delta in (F(1), F(99, 100), F(9, 10), F(1, 2)):
        vals = []
        for q in DETV[:6]:
            if delta == 1:
                continue
            v = occupancy(peq, list(q), delta, x_from_openings(F(1), F(1)))
            vals.append(v[2] - v[1])
        if vals:
            print(f"    delta={float(delta):.2f}  v_DC - v_CD over 6 q's: "
                  f"{[str(t) for t in vals]}")
    print("    (at delta = 1 all of these are exactly 0 -- the equaliser; at delta < 1 they")
    print("     spread out, so no single constant is enforced.)")

    # ------------------------------------------------- (b) the existence region in (u,v)
    print()
    print("=" * 78)
    print("(b')  WHERE DO DISCOUNTED FRIENDLY RIVALS EXIST?   brute force over (p0, p)")
    print("=" * 78)
    print("  Efficient (discounted, eps -> 0): pi_delta(p,p) = Emax = max{1, (1+u+v)/2},")
    print("  with the round-0 state law of SELF-play x = (p0^2, p0(1-p0), (1-p0)p0, (1-p0)^2).")
    print("  Rival: min over (q, q0) of sign*(v_DC - v_CD) >= 0, sign = +1 if T > S else -1.")
    print()

    def is_efficient(p, p0, delta, u, vv):
        x = x_from_openings(p0, p0)
        v = occupancy(list(p), list(p), delta, x)
        pi = dot(v, [F(1), F(u), 1 + F(vv), F(0)])
        Emax = max(F(1), (1 + F(u) + F(vv)) / 2)
        return pi == Emax

    def is_rival(p, p0, delta, sign):
        for q in DETV:
            for q0 in (F(0), F(1)):
                v = occupancy(list(p), list(q), delta, x_from_openings(p0, q0))
                if sign * (v[2] - v[1]) < 0:
                    return False
        return True

    # candidate focal strategies: all 16 vertices x {p0=0,1}, plus a rational grid
    grid = [F(0), F(1, 3), F(1, 2), F(2, 3), F(1)]
    CAND = [(p0, p) for p0 in (F(0), F(1)) for p in DETV]
    rng6 = random.Random(2718)
    CAND += [(F(1), (F(1), rq(rng6), rq(rng6), rq(rng6))) for _ in range(12)]
    CAND += [(rq(rng6), tuple(rq(rng6) for _ in range(4))) for _ in range(12)]
    CAND += [(F(1), (F(1), a, b, F(0))) for a in grid for b in grid]

    delta = F(9, 10)
    print(f"  delta = {delta};  '.' = none, 'F' = a friendly rival exists")
    print(f"  wedge S = {{u+v <= 1, u-v >= 1}};  E = {{u+v>1,u-v>1}};  T>S <=> u-v<1")
    print()
    us = [F(k, 4) for k in range(0, 13)]        # u = 0 .. 3
    vs = [F(k, 4) for k in range(6, -7, -1)]    # v = 1.5 .. -1.5
    print("        u:  " + " ".join(f"{float(u):>4.2f}" for u in us))
    for vv in vs:
        row = []
        for u in us:
            TmS = 1 + vv - u
            if TmS == 0:
                row.append("  = ")      # T = S, rivalry vacuous
                continue
            sign = 1 if TmS > 0 else -1
            found = any(is_efficient(p, p0, delta, u, vv) and is_rival(p, p0, delta, sign)
                        for p0, p in CAND)
            row.append("  F " if found else "  . ")
        print(f"  v={float(vv):>5.2f}: " + "".join(row))
    print()
    print("  Compare the UNDISCOUNTED result: friendly rivals exist exactly on the closed")
    print("  wedge S = {u-v >= 1, u+v <= 1} (plus the T=S line, where rivalry is vacuous).")
    print("  The '=' column is T = S (u - v = 1); on it every efficient strategy counts.")

    print()
    print("  Witness in wedge S: p = ALLC with p0 = 1.")
    for (u, vv) in ((F(3, 2), F(0)), (F(5, 4), F(-1, 2)), (F(1), F(0))):
        d = F(9, 10)
        eff = is_efficient((F(1),) * 4, F(1), d, u, vv)
        riv = is_rival((F(1),) * 4, F(1), d, -1 if 1 + vv - u < 0 else 1)
        print(f"    (u,v)=({u},{vv})  T-S={1+vv-u}  efficient={eff}  rival={riv}")

    print()
    print("  Above the switch line (u+v > 1) discounted efficiency is EMPTY: Emax needs")
    print("  perfect alternation from round 0, but self-play round 0 has")
    print("  x_CC = p0^2 and x_DD = (1-p0)^2, which cannot both vanish.")
    u, vv = F(3, 2), F(1)                    # u+v = 5/2 > 1, u-v = 1/2 < 1 (wedge N, T>S)
    best = None
    for p0, p in CAND:
        v = occupancy(list(p), list(p), F(9, 10), x_from_openings(p0, p0))
        pi = dot(v, [F(1), u, 1 + vv, F(0)])
        if best is None or pi > best[0]:
            best = (pi, p0, p)
    print(f"    (u,v)=({u},{vv}): Emax = {max(F(1), (1+u+vv)/2)}, best self-play over "
          f"{len(CAND)} candidates = {best[0]} = {float(best[0]):.4f}")

    # ---------------------------------------------- the discounted friendly-rival volume
    print()
    print("=" * 78)
    print("(b'')  SIZE OF THE DISCOUNTED FRIENDLY-RIVAL SET INSIDE THE FACET p1 = 1")
    print("=" * 78)
    print("  Friendly rival (wedge S interior, delta < 1, eps -> 0):")
    print("      p0 = 1,  p1 = 1,  p2 + p3 >= (2 delta - 1)/delta,  p4 FREE.")
    print("  Volume in the (p2,p3,p4) unit cube:  1                      delta <= 1/2")
    print("                                       1 - (2delta-1)^2/(2 delta^2)  delta > 1/2")
    print()
    print(f"  {'delta':>7}  {'threshold (2d-1)/d':>20}  {'volume':>10}  {'MC check':>10}")
    rng7 = random.Random(161803)
    for delta in (F(1, 4), F(1, 2), F(3, 5), F(3, 4), F(9, 10), F(99, 100), F(1)):
        thr = (2 * delta - 1) / delta
        vol = F(1) if thr <= 0 else 1 - thr ** 2 / 2
        n = 40000
        hit = sum(1 for _ in range(n)
                  if (rng7.random() + rng7.random()) >= float(thr))
        print(f"  {float(delta):>7.3f}  {float(thr):>20.4f}  {float(vol):>10.5f}  "
              f"{hit / n:>10.5f}")
    print("  delta -> 1 gives 1/2, the project's undiscounted volume.  The set SHRINKS")
    print("  monotonically to the undiscounted one; it never disappears.")

    print()
    print("  Boundary case p = (1,1,0,0) with p0 = 1 -- EXCLUDED undiscounted by the")
    print("  'not (p2=1 and p3=0)' clause, but a rival for every delta < 1:")

    def stationary_exact(p, q, eps):
        M = transition(p, q, eps)
        A = [[M[j][i] - (F(1) if i == j else F(0)) for j in range(4)] for i in range(3)] \
            + [[F(1)] * 4]
        return solve(A, [F(0), F(0), F(0), F(1)])

    pb = [F(1), F(1), F(0), F(0)]
    for delta in (F(1, 2), F(9, 10), F(99, 100), F(999, 1000)):
        m = worst_margin(tuple(pb), delta, F(1))
        print(f"    delta={float(delta):>6.3f}  min_q (v_CD - v_DC) = {float(m):>12.8f}"
              f"   rival={m >= 0}")
    for eps in (F(1, 100), F(1, 1000), F(1, 10000)):
        worst = min(stationary_exact(pb, list(q), eps)[1] - stationary_exact(pb, list(q), eps)[2]
                    for q in DETV)
        print(f"    delta=1, eps={float(eps):<8} min_q (W_CD - W_DC) = {float(worst):>12.8f}"
              f"   rival={worst >= 0}")
    print("  with eps > 0 AND delta < 1 (p1 clamps to 1-eps < 1, so p1 = 1 fails):")
    for eps in (F(1, 100), F(1, 1000)):
        best = None
        for q in DETV:
            for q0 in (F(0), F(1)):
                v = occupancy(pb, list(q), F(9, 10), x_from_openings(F(1), q0), eps)
                m = v[1] - v[2]
                best = m if best is None or m < best else best
        print(f"    delta=0.9, eps={float(eps):<8} min_q (v_CD - v_DC) = {float(best):>12.8f}"
              f"   rival={best >= 0}")
    print("  So the undiscounted exclusion of (p2,p3)=(1,0) is an ORDER-OF-LIMITS artefact:")
    print("    lim_{eps->0} at fixed delta < 1   -> rival (margin exactly 0)")
    print("    lim_{delta->1} at fixed eps > 0   -> not a rival (margin -> -1/2)")
    print("  The project's undiscounted theory is lim_{eps->0} lim_{delta->1}; the")
    print("  discounted theory here is lim_{eps->0} at delta < 1.  They do not commute on")
    print("  this measure-zero face, and only there.")


if __name__ == "__main__":
    main()
