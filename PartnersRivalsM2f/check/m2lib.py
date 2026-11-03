"""Shared harness for the memory-two continuum (PartnersRivalsM2f).

Conventions inherited verbatim from ../../PartnersRivals and ../../PartnersRivalsM1f:

  base game (R,S,T,P) = (1, u, 1+v, 0);   line T = S is u - v = 1;  switch line u + v = 1.
  outcome code o = 2*(focal defects) + (co-player defects):  CC=0, CD=1, DC=2, DD=3
      ("own action first" -- CD means focal C, co-player D).
  memory-two state  s = 4*recent + before  in 0..15   (recent = o_t, before = o_{t-1}).
  a strategy is p in [0,1]^16, p[s] = Pr(focal plays C at round t+1 | s_t = s).
  binary genome: bit s of a 16-bit integer is 1 for C, so ALLC = 65535, ALLD = 0.
  convention A: p is fixed, error clamp p -> eps + (1-2eps) p on BOTH players, eps -> 0+.
"""
import numpy as np

NS = 16

# recent(s), before(s), own action of the focal player in the recent outcome
RECENT = np.array([s >> 2 for s in range(NS)])
BEFORE = np.array([s & 3 for s in range(NS)])
PHI    = np.array([1.0 if (s >> 2) in (0, 1) else 0.0 for s in range(NS)])   # focal played C at t
GVEC   = np.array([1.0 if (s >> 2) == 2 else (-1.0 if (s >> 2) == 1 else 0.0)
                   for s in range(NS)])                                      # 1{DC} - 1{CD}

def genome_to_p(g):
    return np.array([(g >> s) & 1 for s in range(NS)], dtype=float)

def p_to_genome(p):
    g = 0
    for s in range(NS):
        if p[s] > 0.5:
            g |= 1 << s
    return g

# ---------------------------------------------------------------- MDP transitions
# From state s the focal plays C w.p. p[s]. The co-player picks a in {C,D} (a=0 -> C).
# outcome o' = 2*(focal defects) + a ;  next state s' = 4*o' + recent(s).
def next_states(s):
    """returns ((sC_focC, sC_focD), (sD_focC, sD_focD)) : next state for co-action a, focal C/D."""
    r = s >> 2
    return ((4 * 0 + r, 4 * 2 + r),      # a = C : o' = CC or DC
            (4 * 1 + r, 4 * 3 + r))      # a = D : o' = CD or DD

NEXT = [next_states(s) for s in range(NS)]

def occupation_lp(p, sense=+1, eps=0.0):
    """Exact LP:  min over ALL co-player behaviours of  sense * (W_DC - W_CD).

    The pair process is a Markov decision process whose controller is the co-player:
    state s in 0..15, action a in {C,D}, transition given by p.  The set of achievable
    long-run (Cesaro) state-action frequencies of a finite MDP is exactly the polytope

        x >= 0,  sum x = 1,  sum_a x[s',a] = sum_{s,a} x[s,a] P(s'|s,a)  for every s',

    over ALL policies -- stationary, randomised, history-dependent, and multichain.  So the
    minimum below is attained and is the exact worst case over the whole co-player class
    Q_A; no outer bound, no witness search.  Returns (value, x).
    """
    pv = eps + (1 - 2 * eps) * np.asarray(p, dtype=float)
    nx = 2 * NS
    # A_eq: 16 balance rows + 1 normalisation
    A = np.zeros((NS + 1, nx))
    b = np.zeros(NS + 1)
    for s in range(NS):
        for a in (0, 1):
            col = 2 * s + a
            A[s, col] += 1.0                       # out-flow
    for s in range(NS):
        for a in (0, 1):
            col = 2 * s + a
            sC, sD = NEXT[s][a]
            A[sC, col] -= pv[s]
            A[sD, col] -= (1.0 - pv[s])
    A[NS, :] = 1.0
    b[NS] = 1.0
    c = np.zeros(nx)
    for s in range(NS):
        for a in (0, 1):
            c[2 * s + a] = sense * GVEC[s]
    from scipy.optimize import linprog
    r = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if not r.success:
        raise RuntimeError("LP failed: " + r.message)
    return r.fun, r.x

def is_rival(p, sense=+1, tol=1e-9, eps=0.0):
    v, _ = occupation_lp(p, sense, eps)
    return v >= -tol
