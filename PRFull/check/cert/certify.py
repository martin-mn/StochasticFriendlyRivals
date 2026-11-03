"""Certificate primitives for friendly rivals of the memory-two cube (eps -> 0 before any inequality).

THE STRUCTURAL LEMMA.  Fix a binary co-player q and a support stratum sigma (coordinates 0, 1 or
free; free ones x in (0,1)^d).  At eps > 0 the 16-state pair chain is irreducible and
pi_k = A_k / sum_j A_j with A_k the sum over spanning trees rooted at k of the product of edge
probabilities (Markov-chain tree theorem).  The edge k -> (4o+x) has probability (focal factor)
(co-player factor); the focal factor out of state k uses only x_k, and a tree uses exactly one
outgoing edge per non-root state, so each coordinate appears at most once per tree.  The eps-order
of a tree depends only on the pattern (and q), its leading coefficient is a product of x_s, 1-x_s
or 1 over distinct s.  With E = E(sigma, q) the least tree order, the limit is w_k = A*_k(x)/D*(x),
A*_k the eps^E coefficient of A_k (multilinear, nonnegative in the basis {x_s, 1-x_s}), and
D* = sum_k A*_k > 0 on the open stratum.  Every limit margin against q has the sign of the
multilinear polynomial N_q = D* * margin.
  L1  the support of the limit (hence the efficiency class, also in self-play where the lemma's
      multilinearity does not apply) is constant on an open stratum;
  L2  over a closed box inside the open stratum min N_q is attained at a vertex;
  L3  at a corner c of the closed stratum cube: N_q(c) = D*(c) * (binary margin of (c, q)) if the
      corner's own least tree order E_c equals E, and N_q(c) = D*(c) = 0 if E_c > E (E_c < E is
      impossible: pinning only raises tree orders);
  L4  N_q multilinear, <= 0 at every corner and < 0 at one interior point => N_q < 0 on the whole
      open cube (an interior zero would be an interior maximum; along each coordinate line N_q is
      affine, so it would vanish identically).  Symmetrically N_q >= 0 at every corner => N_q >= 0
      on the closed cube, so the margin is >= 0 on the whole open stratum.

Binary corners: for a binary pair every edge's leading coefficient is 1, so A*_k is the NUMBER of
minimum-order spanning trees rooted at k and D* * margin is an integer.  An interval enclosure of a
binary corner margin containing 0 with D*_hi * max|enclosure| < 1 therefore proves it is exactly 0.
(The '1/37 gap' of facts_m2.md and ms.tex:326 is false -- 1/963 occurs -- and is not used.)
Every interior-point margin used as a certificate is either an interval enclosure strictly of one
sign or an exact rational.  Competitiveness is quantified over the 65536 binary co-players; that
this covers every memory-two co-player is the paper's box-vertex lemma (it also follows from the
same structural lemma applied to q with the focal x held fixed).

Every function returns a dict whose key 'certified' (or 'excluded') is True only when the proof
went through; the remaining keys are the evidence.
"""
import math, time, itertools
import numpy as np
from fractions import Fraction
from multiprocessing import Pool
from . import engine as en
from . import exact as ex


# ------------------------------------------------------------------ helpers
def _prod(it):
    r = 1.0
    for v in it: r *= v
    return r

def _sign_eff(wedge=None, sign=None, eff=None):
    if wedge is not None:
        s, e = en.WEDGES[wedge]; return s, e
    return sign, eff

def _exact_job(args):
    p, g, sign = args
    return ex.margin(p, ex.g2p(int(g)), sign)

def exact_margins(points, qs, sign, nproc=1):
    """exact margins sign*(w_DC-w_CD) for the pairs (points[i], qs[i]); points are sequences of
    rationals.  Uses a pool of nproc processes (<= engine.MAX_THREADS) for long lists."""
    jobs = [(list(p), int(g), sign) for p, g in zip(points, qs)]
    nproc = max(1, min(nproc, en.MAX_THREADS))
    if nproc == 1 or len(jobs) < 2000:
        return [_exact_job(j) for j in jobs]
    with Pool(nproc) as pool:
        return pool.map(_exact_job, jobs, chunksize=256)

def _full_point(pat, vals):
    """a 16-vector of Fractions: pinned coordinates from pat, free ones from vals (length d or 16)"""
    pat = en.parse_pattern(pat); fr = en.free_coords(pat)
    vals = list(vals)
    if len(vals) == 16: vals = [vals[j] for j in fr]
    assert len(vals) == len(fr)
    x = [Fraction(int(v)) for v in pat]
    for j, v in zip(fr, vals): x[j] = Fraction(v)
    return x

def check_efficiency(pat, eff):
    """exact (L1): support class of the self-play limit on the open stratum and whether it meets eff"""
    cls = en.efficiency_class(pat)
    return dict(efficiency_class=cls, efficient=(True if eff is None else en.meets(cls, eff)))

def corner_cost_seconds(pat, nq=en.NQ, nthreads=None, us_per_eval=2.6):
    """estimated wall time of a corner scan (2^d interval evaluations per co-player).  The default
    rate is the measured single-thread rate of the x86_64 slice under Rosetta on Martin's Mac."""
    d = len(en.free_coords(pat)); nt = en._nt(nthreads or en.MAX_THREADS)
    return (1 << d)*nq*us_per_eval*1e-6/nt

# ------------------------------------------------------------------ single points
def certify_point(p, sign=None, eff=None, wedge=None, nthreads=None, nproc=1):
    """PROVES, for the single rational point p (its free coordinates strictly inside (0,1)):
    its exact efficiency class (L1), and for every binary co-player q the exact sign of the margin
    sign*(w_DC - w_CD)(p, q): interval enclosure strictly signed, or structural zero (w_CD = w_DC = 0
    by exponents), or -- when the enclosure contains 0 -- the exact rational value.
    'certified' is True iff p is efficient (if eff given) and no co-player has a negative margin.
    Also returns exact self-play frequencies, min margin and the counts
    (positive / structural zero / exact nonstructural zero ('identities' at p) / negative)."""
    sign, eff = _sign_eff(wedge, sign, eff)
    t0 = time.time(); nt = nthreads or en.MAX_THREADS
    P = [Fraction(x) for x in p]
    plo, phi = en.point_box(P)
    S = en.i_scan(plo, phi, sign, nthreads=nt)
    cls = S['cls']
    un = np.where(cls == en.CLS_UNRES)[0]
    vals = exact_margins([P]*len(un), un, sign, nproc=nproc)
    zero_un = sum(1 for v in vals if v == 0); neg_un = sum(1 for v in vals if v < 0); pos_un = len(vals) - zero_un - neg_un
    npos = int((cls == 1).sum()) + pos_un; nneg = int((cls == -1).sum()) + neg_un
    nstruct = int((cls == 0).sum())
    self_W = ex.selfplay(P)
    ecls = en.support_class(ex.lead_exact(P, P)['supp'])
    # min margin: exact.  0 if any zero; otherwise the exact minimum over every q whose lower
    # enclosure end is within 1e-9 of the smallest (enclosure widths are ~1e-12)
    if nneg:
        negq = [int(q) for q in np.where(cls == -1)[0]] + [int(un[i]) for i, v in enumerate(vals) if v < 0]
        m = min((ex.margin(P, ex.g2p(q), sign), q) for q in negq[:64])
        minmargin, argmin = m
    elif nstruct or zero_un:
        minmargin, argmin = Fraction(0), None
    else:
        m0 = S['mlo'].min(); near = np.where(S['mlo'] <= m0 + 1e-9)[0]
        minmargin, argmin = min((ex.margin(P, ex.g2p(int(q)), sign), int(q)) for q in near)
    out = dict(certified=(nneg == 0 and (eff is None or en.meets(ecls, eff))),
               efficiency_class=ecls, selfplay=self_W, n_pos=npos, n_struct_zero=nstruct,
               n_exact_zero=zero_un, n_zero=nstruct + zero_un, n_neg=nneg, n_unresolved_by_interval=len(un),
               min_margin=minmargin, argmin=argmin, seconds=time.time() - t0)
    return out

# ------------------------------------------------------------------ whole open strata (L3 + L4)
def certify_open_stratum_all_q(pat, sign=None, eff=None, wedge=None, nthreads=None, block=4096,
                               stop_on_fail=True, max_seconds=840):
    """PROVES (if out['certified']): every point of the OPEN stratum pat is competitive -- for every
    binary co-player q and every x in the open stratum, sign*(w_DC - w_CD)(x, q) >= 0 -- and, if eff
    is given, efficient of class eff (L1).
    Method (L3 + L4, positive direction), per q: E(q) of the open stratum (integer, at the centre)
    and, at each of the 2^d corners c, the corner's own E_c(q); every corner with E_c = E must have
    corner margin certified >= 0 (interval > 0, structural 0, or 0 by integrality); corners with
    E_c > E have N_q(c) = 0.  So N_q >= 0 on the closed cube and margin = N_q/D* >= 0 on the open
    stratum.  Nothing is sampled.  Refuses (out['refused']) when the estimated time exceeds
    max_seconds -- run it on Cannon instead."""
    sign, eff = _sign_eff(wedge, sign, eff)
    pat = en.parse_pattern(pat); d = len(en.free_coords(pat)); nt = nthreads or en.MAX_THREADS
    est = corner_cost_seconds(pat, nthreads=nt)
    effd = check_efficiency(pat, eff)
    out = dict(pattern=en.pattern_str(pat), d=d, sign=sign, **effd, estimated_seconds=est)
    if est > max_seconds:
        out.update(certified=False, refused=True); return out
    t0 = time.time(); nfail = 0; failq = []; tot = np.zeros(7, dtype=np.int64); nchk = 0
    for a in range(0, en.NQ, block):
        ql = np.arange(a, min(a + block, en.NQ))
        C = en.corner_scan(pat, sign, ql, nthreads=nt)
        tot += C['cnt'].sum(0); nchk += len(ql)
        bad = np.where(~C['pos_ok'])[0]
        nfail += len(bad)
        for i in bad:
            if len(failq) < 20: failq.append(int(ql[i]))
        if nfail and stop_on_fail: break
    names = ['pos', 'neg', 'zero_struct', 'zero_integral', 'E_c>E', 'E_c<E', 'unresolved']
    out.update(certified=(nfail == 0 and nchk == en.NQ and effd['efficient'] and tot[5] == 0),
               competitive_certified=(nfail == 0 and nchk == en.NQ and tot[5] == 0),
               refused=False, n_q_checked=nchk, n_q_failing=nfail, failing_q=failq,
               corner_counts=dict(zip(names, map(int, tot))), seconds=time.time() - t0)
    return out

def exclude_stratum(pat, sign=None, candidates=None, wedge=None, search_all=True, nthreads=None,
                    max_seconds=840, first_batch=64):
    """PROVES (if out['excluded']): one binary co-player q strictly beats EVERY point of the open
    stratum pat, i.e. sign*(w_DC - w_CD)(x, q) < 0 for all x in it, so the stratum contains no
    competitive (hence no friendly-rival) strategy.
    Method (L3 + L4, negative direction): every corner c with E_c(q) = E(q) has corner margin
    certified <= 0 (so N_q <= 0 on the closed cube), and the margin at the centre is < 0 --
    certified by a strictly negative interval enclosure, or exactly in rationals.  Then N_q < 0 on
    the whole open cube.  (If the corners pass but the centre margin is exactly 0, N_q is
    identically 0 and q does not exclude.)
    Search: the given candidates first (a pool of co-players that worked before: only they are
    scanned, so a pool hit costs ~ms); then, if search_all, every q whose centre margin is
    certified negative or unresolved, most negative first.  Returns the first q that works."""
    sign, _ = _sign_eff(wedge, sign, None)
    pat = en.parse_pattern(pat); nt = nthreads or en.MAX_THREADS; t0 = time.time()
    c = en.centre(pat)                               # centre 1/2 is exact in binary
    cpt = _full_point(pat, [Fraction(1, 2)]*len(en.free_coords(pat)))
    out = dict(pattern=en.pattern_str(pat), d=int((pat == 2).sum()), sign=sign, excluded=False, q=None)
    names = ['pos', 'neg', 'zero_struct', 'zero_integral', 'E_c>E', 'E_c<E', 'unresolved']
    tried = [0]
    def try_batch(batch, cls, mlo, mhi):
        C = en.corner_scan(pat, sign, batch, nthreads=nt); tried[0] += len(batch)
        for i in np.where(C['neg_ok'])[0]:
            q = int(batch[i])
            if cls[i] == en.CLS_NEG:                 # strictly negative enclosure at the centre
                out.update(excluded=True, q=q, centre_margin_enclosure=(float(mlo[i]), float(mhi[i])),
                           centre_margin_exact=None)
            elif cls[i] == en.CLS_UNRES:             # decide the centre exactly
                v = ex.margin(cpt, ex.g2p(q), sign)
                if v >= 0: continue                  # v == 0 with all corners <= 0 means N_q == 0
                out.update(excluded=True, q=q, centre_margin_enclosure=None, centre_margin_exact=v)
            else:
                continue                             # centre margin >= 0: q cannot exclude
            out['corner_counts'] = dict(zip(names, map(int, C['cnt'][i])))
            return True
        return False
    seen = set()
    if candidates is not None and len(candidates):
        cand = []
        for q in candidates:
            if int(q) not in seen: seen.add(int(q)); cand.append(int(q))
        Sc = en.i_scan(c, c, sign, cand, nthreads=nt)
        if try_batch(np.array(cand), Sc['cls'], Sc['mlo'], Sc['mhi']):
            out.update(n_candidates_tried=tried[0], found_in_pool=True, seconds=time.time() - t0); return out
    if search_all:
        S = en.i_scan(c, c, sign, nthreads=nt)
        out['n_centre_negative'] = int((S['cls'] == -1).sum())
        cand = np.where((S['cls'] == -1) | (S['cls'] == 2))[0]
        cand = cand[np.argsort(S['mhi'][cand], kind='stable')]
        cand = np.array([q for q in cand if int(q) not in seen], dtype=np.int64)
        pos = 0
        while pos < len(cand):
            if time.time() - t0 > max_seconds:
                out.update(timed_out=True); break
            batch = cand[pos: pos + (first_batch if pos == 0 else 1024)]; pos += len(batch)
            if try_batch(batch, S['cls'][batch], S['mlo'][batch], S['mhi'][batch]): break
    out.update(n_candidates_tried=tried[0], found_in_pool=False, seconds=time.time() - t0)
    return out

def is_identity(pat, sign, q):
    """PROVES (True) or refutes (False): N_q is identically 0 on the stratum, i.e. the margin
    against q vanishes on the whole open stratum.  A multilinear polynomial is determined by its
    corner values, and by L3 N_q(c) = 0 iff E_c > E or the binary corner margin is 0."""
    Dd = en.corner_detail(pat, sign, q)
    ok = (Dd['Ec'] >= Dd['E']).all()
    same = Dd['Ec'] == Dd['E']
    zero = (Dd['cls'] == 0) | (Dd['cls'] == 3)
    return bool(ok and (zero | ~same).all())

# ------------------------------------------------------------------ boxes inside a stratum (L2)
def _box_bounds(pat, lo, hi):
    pat = en.parse_pattern(pat); fr = en.free_coords(pat)
    lo = list(lo); hi = list(hi)
    if len(lo) == 16: lo = [lo[j] for j in fr]
    if len(hi) == 16: hi = [hi[j] for j in fr]
    lo = [Fraction(x) for x in lo]; hi = [Fraction(x) for x in hi]
    for a, b in zip(lo, hi):
        if not (0 < a <= b < 1): raise ValueError("box must lie inside the open stratum: 0 < lo <= hi < 1")
    return pat, fr, lo, hi

def certify_box(pat, lo, hi, sign=None, eff=None, wedge=None, nthreads=None, use_corners=None,
                nproc=1, max_seconds=840):
    """PROVES (if out['certified']): every point x of the closed box B = prod_s [lo_s, hi_s] over the
    free coordinates of pat (0 < lo_s <= hi_s < 1; pinned coordinates as in pat) is competitive
    against every binary co-player, and (with eff) efficient -- so B, of dimension d and volume
    prod (hi_s - lo_s), lies in the friendly-rival set.
    Per co-player q one of the following holds (each alone implies margin >= 0 on B):
      (a) corner test, L3 + L4: every stratum corner has N_q(c) >= 0, so N_q >= 0 on the closed cube;
      (b) interval evaluation with the whole box as input: enclosure > 0, or structural zero;
      (c) L2: N_q >= 0 at all 2^d vertices of B -- each vertex certified by a positive interval
          enclosure, a structural zero, or an exact rational margin >= 0.
    A certified-negative vertex margin (or a negative exact one) is a counterexample and ends the
    run with certified False.  use_corners: None = automatic (only when the corner scan is cheap)."""
    sign, eff = _sign_eff(wedge, sign, eff)
    pat, fr, lo, hi = _box_bounds(pat, lo, hi); d = len(fr); nt = nthreads or en.MAX_THREADS
    t0 = time.time(); effd = check_efficiency(pat, eff)
    out = dict(pattern=en.pattern_str(pat), d=d, sign=sign, lo=lo, hi=hi, **effd,
               volume=_prod(float(b - a) for a, b in zip(lo, hi)))
    # stage (a)
    if use_corners is None: use_corners = corner_cost_seconds(pat, nthreads=nt) < min(120, max_seconds/3)
    rest = np.arange(en.NQ)
    if use_corners:
        C = en.corner_scan(pat, sign, None, nthreads=nt)
        out['n_by_corners'] = int(C['pos_ok'].sum()); rest = np.where(~C['pos_ok'])[0]
        if (C['cnt'][:, en.C_LOW] > 0).any(): raise RuntimeError("E_c < E encountered: lemma violated")
    else:
        out['n_by_corners'] = 0
    # stage (b)
    blo = en.centre(pat).copy(); bhi = blo.copy()
    for j, a, b in zip(fr, lo, hi): blo[j] = en.enclose(a)[0]; bhi[j] = en.enclose(b)[1]
    B = en.i_scan(blo, bhi, sign, rest, nthreads=nt)
    okb = (B['cls'] == 1) | (B['cls'] == 0)
    out['n_by_box_interval'] = int(okb.sum()); rest = rest[~okb]
    # stage (c)
    nv = 1 << d
    if nv*len(rest)*2.6e-6/nt > max_seconds:
        out.update(certified=False, refused=True, n_left=len(rest)); return out
    V = np.array(list(itertools.product(*[(a, b) for a, b in zip(lo, hi)])), dtype=object)
    vlo = np.repeat(blo[None, :], nv, 0); vhi = vlo.copy()
    for k in range(nv):
        for t, j in enumerate(fr): vlo[k, j], vhi[k, j] = en.enclose(V[k][t])
    out['n_vertex_q'] = len(rest)
    if len(rest):
        VS = en.vertex_scan(vlo, vhi, sign, rest, nthreads=nt, cap=5_000_000)
        if VS['nunres'] > 5_000_000: raise RuntimeError("too many unresolved vertex pairs")
        negv = np.where(VS['cnt'][:, 1] > 0)[0]
        out['n_vertex_neg_interval'] = len(negv)
        if len(negv):
            i = negv[0]; out.update(certified=False, counterexample_q=int(rest[i]), seconds=time.time() - t0)
            return out
        un = VS['unres']
        out['n_vertex_exact'] = len(un)
        vals = exact_margins([_full_point(pat, V[k]) for k in un[:, 1]], rest[un[:, 0]], sign, nproc=nproc)
        neg = [i for i, v in enumerate(vals) if v < 0]
        out['n_vertex_exact_zero'] = sum(1 for v in vals if v == 0)
        if neg:
            i = neg[0]; out.update(certified=False, counterexample_q=int(rest[un[i, 0]]),
                                   counterexample_vertex=[str(x) for x in V[un[i, 1]]], seconds=time.time() - t0)
            return out
        out['min_vertex_margin_lower_bound'] = float(VS['mmin'].min())
    out.update(certified=effd['efficient'], competitive_certified=True, refused=False, seconds=time.time() - t0)
    return out

def certify_ratio_box(c, h=None, sign=None, eff=None, wedge=None, lo=None, hi=None, nthreads=None):
    """PROVES (if out['certified']): every x in the box B around the rational point c -- the free
    coordinates of c's stratum within c_s +- h (or within [lo_s, hi_s] if given), inside (0,1) -- is
    competitive against every binary co-player and (with eff) efficient.
    Ratio route: on the stratum A*_k(x) is a nonnegative combination of monomials prod_s f_s(x_s),
    f_s in {x_s, 1-x_s, 1}, so on B each A*_k(x)/A*_k(c) lies in [L, U] with
    L = prod_s min(lo_s/c_s, (1-hi_s)/(1-c_s)) and U = prod_s max(hi_s/c_s, (1-lo_s)/(1-c_s))
    (= prod (1 -+ h/m_s), m_s = min(c_s, 1-c_s), for a symmetric box).  With fav/unf the favourable
    and unfavourable outcome (DC/CD if sign > 0, CD/DC if sign < 0), margin >= 0 on B for q when
      * w_unf = 0 by exponents (structural, so on the whole stratum), or
      * w_fav(c)/w_unf(c) >= U/L  (interval lower bound first, exact rational if needed), or
      * w_fav(c) = w_unf(c) exactly and N_q is identically 0 (is_identity: all corners zero).
    One scan at the centre, no vertices: this is the route for high-dimensional strata (S, d = 15)."""
    sign, eff = _sign_eff(wedge, sign, eff)
    C = [Fraction(x) for x in c]; pat = en.pattern_of(C); fr = en.free_coords(pat); nt = nthreads or en.MAX_THREADS
    t0 = time.time()
    if lo is None:
        h = Fraction(h); lo = [C[j] - h for j in fr]; hi = [C[j] + h for j in fr]
    else:
        lo = [Fraction(x) for x in (lo if len(lo) == len(fr) else [lo[j] for j in fr])]
        hi = [Fraction(x) for x in (hi if len(hi) == len(fr) else [hi[j] for j in fr])]
    L = Fraction(1); U = Fraction(1)
    for j, a, b in zip(fr, lo, hi):
        if not (0 < a <= C[j] <= b < 1): raise ValueError("box must contain c and lie inside the open stratum")
        L *= min(a/C[j], (1 - b)/(1 - C[j])); U *= max(b/C[j], (1 - a)/(1 - C[j]))
    UL = U/L
    effd = check_efficiency(pat, eff)
    fav, unf = (2, 1) if sign > 0 else (1, 2)
    plo, phi = en.point_box(C)
    S = en.i_scanW(plo, phi, nthreads=nt)
    unf_on = en.outcome_present(S['supp'], unf); fav_on = en.outcome_present(S['supp'], fav)
    n_unf_struct = int((~unf_on).sum()); n_struct_both = int((~unf_on & ~fav_on).sum())
    live = np.where(unf_on)[0]
    fail_fav0 = live[~fav_on[live]]
    rho_lo = np.full(en.NQ, np.inf)
    a = S['Wlo'][live, fav]; b = S['Whi'][live, unf]
    rho_lo[live] = np.nextafter(a/b, -np.inf)
    # compare the float lower bound with the exact U/L: Fraction(rho_lo) >= UL
    ULf = float(UL); ULhi = float(np.nextafter(ULf, np.inf)) if Fraction(ULf) < UL else ULf
    passed_int = live[rho_lo[live] >= ULhi]
    todo = np.setdiff1d(live, passed_int); todo = np.setdiff1d(todo, fail_fav0)
    n_exact = 0; n_ident = 0; fails = [int(q) for q in fail_fav0[:10]]; nfail = len(fail_fav0)
    rho_exact_min = None; rho_arg = None
    for q in todo:
        r = ex.lead_exact(C, ex.g2p(int(q)))['W']; n_exact += 1
        rho = r[fav]/r[unf]
        if rho_exact_min is None or rho < rho_exact_min: rho_exact_min, rho_arg = rho, int(q)
        if rho >= UL: continue
        if rho == 1 and is_identity(pat, sign, int(q)): n_ident += 1; continue
        nfail += 1; fails.append(int(q))
    # the exact minimum of rho: evaluate exactly every q whose interval lower bound is within 1e-9
    # (relative) of the smallest lower bound (enclosure widths are ~1e-12)
    if len(live) and np.isfinite(rho_lo[live]).any():
        m0 = rho_lo[live].min()
        for q in live[rho_lo[live] <= m0*(1 + 1e-9) + 1e-12]:
            r = ex.lead_exact(C, ex.g2p(int(q)))['W']
            if r[unf] != 0:
                rq = r[fav]/r[unf]
                if rho_exact_min is None or rq < rho_exact_min: rho_exact_min, rho_arg = rq, int(q)
    out = dict(pattern=en.pattern_str(pat), d=len(fr), sign=sign, **effd, U_over_L=UL,
               n_struct_zero=n_struct_both, n_unfav_struct_zero=n_unf_struct - n_struct_both,
               n_ratio_interval=len(passed_int), n_ratio_exact=n_exact, n_identities=n_ident,
               n_fail=nfail, failing_q=fails[:10], rho_min=rho_exact_min, rho_argmin=rho_arg,
               volume=_prod(float(b - a) for a, b in zip(lo, hi)), seconds=time.time() - t0)
    out['competitive_certified'] = (nfail == 0)
    out['certified'] = (nfail == 0) and effd['efficient']
    return out
