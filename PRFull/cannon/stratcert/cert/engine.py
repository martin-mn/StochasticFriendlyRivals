"""ctypes wrapper for libcert.dylib (certeng.c): float and interval leading-order engines.

Conventions (as in PartnersRivalsM2f): p_s is the probability of cooperating at state s = 4x + y,
x the last outcome, y the one before, CC=0 CD=1 DC=2 DD=3, focal action first.  A binary
co-player is a genome g in 0..65535 with bit s = q_s (1 = C): ALLC = 65535, ALLD = 0.
A pattern (support stratum) is a length-16 int8 array with 0, 1, or 2 = free, or a string such as
'1*1*001*1*1**0*0' (character s is p_s).  The competitive margin is sign*(w_DC - w_CD), sign = +1
where T > S and -1 where T < S.

What is exact here and what is not:
  * exponents (E, Ek, the support, hence the efficiency class) are integer computations that never
    read a coefficient value -- exact in both engines, for any input;
  * float coefficients (f_*) carry rounding (relative ~1e-15 on typical inputs);
  * interval coefficients (i_*) are outward-rounded enclosures: the exact limit value lies inside
    [lo, hi] for every focal point of the input box and the given co-player.
"""
import ctypes, os, sys, math
import numpy as np
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))
_L = ctypes.CDLL(os.path.join(_HERE, "libcert.dylib" if sys.platform == "darwin" else "libcert.so"))
_d = ctypes.POINTER(ctypes.c_double); _i = ctypes.POINTER(ctypes.c_int)
_l = ctypes.POINTER(ctypes.c_long); _b = ctypes.POINTER(ctypes.c_byte)
_L.f_eval.argtypes = [_d, _d, _d, _i, _i, _d, _i]
_L.i_eval.argtypes = [_d, _d, _d, _d, _d, _d, _d, _i, _i, _d, _d, _i]; _L.i_eval.restype = ctypes.c_int
_L.f_scan.argtypes = [_d, ctypes.c_long, _l, _d, _i, _i, _d, ctypes.c_int]
_L.i_scan.argtypes = [_d, _d, ctypes.c_int, ctypes.c_long, _l, _i, _d, _d, _i, _i, ctypes.c_int]
_L.i_scan.restype = ctypes.c_int
_L.corner_scan.argtypes = [_b, ctypes.c_int, ctypes.c_long, _l, _i, _i, _d, _d, ctypes.c_int]
_L.corner_detail.argtypes = [_b, ctypes.c_int, ctypes.c_long, _i, _i, _d, _d, _d]; _L.corner_detail.restype = ctypes.c_int
_L.vertex_scan.argtypes = [_d, _d, ctypes.c_long, ctypes.c_int, ctypes.c_long, _l, _i, _d, _l, ctypes.c_long, _l, ctypes.c_int]
_L.vertex_scan.restype = ctypes.c_int
_L.f_self.argtypes = [_d, _d, _i, _i]

NQ = 65536
# NOTE: the '1/37 gap' of facts_m2.md / ms.tex:326 (smallest nonzero |w_DC - w_CD| over binary
# pairs) is FALSE: p=15204, q=21477 gives w_DC - w_CD = 1/963 exactly (tests.py (d)).  The
# rigorous replacement used here is integrality: for a binary pair every edge's leading
# coefficient is 1, so A*_k = #minimum-order spanning trees rooted at k and D* * (w_DC - w_CD) is
# an integer; a margin enclosure with D*_hi * max|enclosure| < 1 proves the margin is exactly 0.
# Local-compute rule: at most 3 cores on Martin's Mac.  On Cannon set CERT_MAX_THREADS.
MAX_THREADS = int(os.environ.get("CERT_MAX_THREADS", "3"))
WEDGES = {'W': (+1, 'C'), 'S': (-1, 'C'), 'N': (+1, 'A'), 'E': (-1, 'A'),
          'ray+': (+1, 'Z'), 'ray-': (-1, 'Z')}      # (sign of T-S, efficiency class required)

def _nt(n):
    return max(1, min(int(n), MAX_THREADS))
def _dp(a): return a.ctypes.data_as(_d)
def _ip(a): return a.ctypes.data_as(_i)
def _lp(a): return a.ctypes.data_as(_l)
def _f64(a): return np.ascontiguousarray(a, dtype=np.float64)

# ------------------------------------------------------------------ patterns and points
def parse_pattern(s):
    """'1*1*001*1*1**0*0' | sequence of 0/1/2 or 0/1/'*' | census code (int) -> int8[16] of 0,1,2(free)"""
    if isinstance(s, (int, np.integer)):
        g = int(s); d = []
        for _ in range(16): d.append(g % 3); g //= 3
        return np.array([{0: 0, 1: 2, 2: 1}[x] for x in d], dtype=np.int8)   # census digit 1 = free, 2 = one
    if isinstance(s, str):
        s = s.strip(); assert len(s) == 16, s
        return np.array([2 if ch == '*' else int(ch) for ch in s], dtype=np.int8)
    a = np.array([2 if (x == '*' or x is None) else int(x) for x in s], dtype=np.int8)
    assert a.shape == (16,) and set(a.tolist()) <= {0, 1, 2}
    return a
def pattern_str(pat): return ''.join('01*'[int(x)] for x in parse_pattern(pat))
def census_code(pat):
    """inverse of parse_pattern(int): census digits 0 -> 0, 1 -> free, 2 -> one"""
    pat = parse_pattern(pat); return int(sum({0: 0, 2: 1, 1: 2}[int(x)]*3**j for j, x in enumerate(pat)))
def pattern_of(p):
    """support pattern of a point"""
    return np.array([0 if x == 0 else (1 if x == 1 else 2) for x in p], dtype=np.int8)
def free_coords(pat): return [j for j, x in enumerate(parse_pattern(pat)) if x == 2]
def centre(pat):
    """the point with every free coordinate 1/2 (exact in binary)"""
    pat = parse_pattern(pat); return np.where(pat == 2, 0.5, pat.astype(float))
def corners(pat):
    """2^d binary corners of the closed stratum cube; row m has free coordinate t equal to bit t of m"""
    pat = parse_pattern(pat); fr = free_coords(pat); base = np.where(pat == 1, 1.0, 0.0)
    C = np.repeat(base[None, :], 1 << len(fr), 0)
    m = np.arange(1 << len(fr))
    for t, j in enumerate(fr): C[:, j] = (m >> t) & 1
    return C
def g2p(g): return np.array([(int(g) >> j) & 1 for j in range(16)], dtype=float)
def p2g(p):
    assert all(x in (0, 1) for x in p); return int(sum(int(x) << j for j, x in enumerate(p)))
def _qvec(q): return g2p(q) if isinstance(q, (int, np.integer)) else _f64(q)

# efficiency class from the support (exact)
REC = [k >> 2 for k in range(16)]
def support_class(supp):
    """'C' (w_CC = 1), 'A' (w_CD + w_DC = 1), 'Z' (w_DD = 0 only), '-' (none); most specific first"""
    r = [0, 0, 0, 0]
    for k in range(16):
        if supp >> k & 1: r[k >> 2] = 1
    if not (r[1] or r[2] or r[3]): return 'C'
    if not (r[0] or r[3]): return 'A'
    if not r[3]: return 'Z'
    return '-'
def meets(cls, eff):
    """does support class cls satisfy the efficiency requirement eff ('C', 'A' or 'Z')?"""
    return cls == eff or (eff == 'Z' and cls in ('C', 'A'))

# ------------------------------------------------------------------ enclosures of rationals
def enclose(x):
    """floats (lo, hi) with lo <= x <= hi for a rational (Fraction/int/str/float) x"""
    fx = Fraction(x); f = float(fx); r = Fraction(f)
    if r == fx: return f, f
    return (f, float(np.nextafter(f, np.inf))) if r < fx else (float(np.nextafter(f, -np.inf)), f)
def point_box(p):
    """interval enclosure (plo, phi) of a rational point p (pinned 0/1 stay exact)"""
    lo = np.empty(16); hi = np.empty(16)
    for j, x in enumerate(p): lo[j], hi[j] = enclose(x)
    return lo, hi

# ------------------------------------------------------------------ float engine
def f_eval(p, q):
    """float leading-order limit of focal p (point in [0,1]^16) against q (genome or point).
    Returns dict: w16, W (CC,CD,DC,DD), Ek (16 ints), E (int: min spanning-tree exponent), D (float D*),
    supp (int bitmask of the limit support).  ~5-10 us per call from python."""
    p = _f64(p); q = _qvec(q); w = np.zeros(16); Ek = np.zeros(16, dtype=np.int32)
    E = ctypes.c_int(); D = ctypes.c_double(); s = ctypes.c_int()
    _L.f_eval(_dp(p), _dp(q), _dp(w), _ip(Ek), ctypes.byref(E), ctypes.byref(D), ctypes.byref(s))
    W = w.reshape(4, 4).sum(1)
    return dict(w16=w, W=W, Ek=Ek, E=E.value, D=D.value, supp=s.value)

def f_scan(p, qlist=None, nthreads=1):
    """float limits of p against every co-player in qlist (default all 65536 genomes).
    Returns dict: q (int64), W (n x 4), E (n), supp (n), D (n).  ~70 ms per 65536 (1 thread, native)."""
    p = _f64(p)
    ql = None if qlist is None else np.ascontiguousarray(qlist, dtype=np.int64)
    n = NQ if ql is None else len(ql)
    W = np.zeros((n, 4)); E = np.zeros(n, dtype=np.int32); s = np.zeros(n, dtype=np.int32); D = np.zeros(n)
    _L.f_scan(_dp(p), n, None if ql is None else _lp(ql), _dp(W), _ip(E), _ip(s), _dp(D), _nt(nthreads))
    return dict(q=np.arange(NQ) if ql is None else ql, W=W, E=E, supp=s, D=D)

def f_self(p):
    """self-play: float (w_CC,..,w_DD), exact support bitmask, E"""
    p = _f64(p); W = np.zeros(4); s = ctypes.c_int(); E = ctypes.c_int()
    _L.f_self(_dp(p), _dp(W), ctypes.byref(s), ctypes.byref(E))
    return W, s.value, E.value

def efficiency_class(pat):
    """exact efficiency class ('C','A','Z','-') of the open stratum pat (pattern string/array or
    census code).  By lemma L1 the support of the self-play limit is decided by exponents alone,
    which depend only on the pattern, so one evaluation at the centre decides it for every point
    of the open stratum -- no tolerance involved."""
    return support_class(f_self(centre(parse_pattern(pat)))[1])

def efficiency_class_of_point(p):
    """efficiency class of the stratum containing the point p (same as at p itself)"""
    return efficiency_class(pattern_of(p))

# ------------------------------------------------------------------ interval engine
def _box(plo, phi):
    plo = _f64(plo); phi = plo if phi is None else _f64(phi)
    return plo, phi

def i_eval(plo, phi, q):
    """interval leading-order limit, valid for every focal point x with plo <= x <= phi (free
    coordinates strictly inside (0,1); pinned ones exactly 0 or 1) against co-player q.
    Returns dict: wlo, whi (16), Wlo, Whi (4), Ek, E, Dlo, Dhi, supp.  Raises on an invalid box."""
    plo, phi = _box(plo, phi); q = _qvec(q)
    wlo = np.zeros(16); whi = np.zeros(16); Wlo = np.zeros(4); Whi = np.zeros(4); Ek = np.zeros(16, dtype=np.int32)
    E = ctypes.c_int(); Dlo = ctypes.c_double(); Dhi = ctypes.c_double(); s = ctypes.c_int()
    r = _L.i_eval(_dp(plo), _dp(phi), _dp(q), _dp(wlo), _dp(whi), _dp(Wlo), _dp(Whi), _ip(Ek),
                  ctypes.byref(E), ctypes.byref(Dlo), ctypes.byref(Dhi), ctypes.byref(s))
    if r: raise ValueError("invalid interval input (code %d): free coordinates must lie strictly in (0,1)" % r)
    return dict(wlo=wlo, whi=whi, Wlo=Wlo, Whi=Whi, Ek=Ek, E=E.value, Dlo=Dlo.value, Dhi=Dhi.value, supp=s.value)

CLS_STRUCT, CLS_POS, CLS_NEG, CLS_UNRES = 0, 1, -1, 2
def i_scan(plo, phi=None, sign=+1, qlist=None, nthreads=1):
    """certified-sign scan of the margin sign*(w_DC - w_CD) over co-players (default all 65536).
    For each q: cls = +1 certified > 0 on the whole input box, -1 certified < 0, 0 structural zero
    (w_CD = w_DC = 0 by exponents alone, hence on the whole open stratum), 2 unresolved (enclosure
    contains 0: an exact tie or too close to call).  Also mlo, mhi (enclosure), E, supp.
    ~0.1-0.2 s per 65536 (1 thread)."""
    plo, phi = _box(plo, phi)
    ql = None if qlist is None else np.ascontiguousarray(qlist, dtype=np.int64)
    n = NQ if ql is None else len(ql)
    cls = np.zeros(n, dtype=np.int32); mlo = np.zeros(n); mhi = np.zeros(n)
    E = np.zeros(n, dtype=np.int32); s = np.zeros(n, dtype=np.int32)
    r = _L.i_scan(_dp(plo), _dp(phi), int(sign), n, None if ql is None else _lp(ql), _ip(cls), _dp(mlo), _dp(mhi),
                  _ip(E), _ip(s), _nt(nthreads))
    if r: raise ValueError("invalid interval input")
    return dict(q=np.arange(NQ) if ql is None else ql, cls=cls, mlo=mlo, mhi=mhi, E=E, supp=s)

# columns of corner_scan counts
C_POS, C_NEG, C_ZS, C_ZI, C_HIGH, C_LOW, C_UNR = range(7)
def corner_scan(pat, sign=+1, qlist=None, nthreads=1):
    """Lemma L3 corner table.  For each co-player q (default all 65536) and each of the 2^d corners
    c of the closed stratum cube: compare the corner's integer minimum tree exponent E_c(q) with the
    open stratum's E(q), and certify the sign of the binary corner margin (interval engine; an
    enclosure straddling 0 with D*_hi * max|enclosure| < 1 is exactly 0 by integrality, since
    D* * margin is an integer for a binary pair).
    Returns dict: q, E (n), cnt (n x 7: npos, nneg, nzero_struct, nzero_int, nhigh, nlow, nunres),
    mmin, mmax, and
      pos_ok[i]: every corner has N_q(c) >= 0  (no certified-negative corner with E_c = E, nothing
                 unresolved, no E_c < E)  => N_q >= 0 on the closed cube (L2), i.e. margin >= 0 on
                 the whole open stratum against this q;
      neg_ok[i]: every corner has N_q(c) <= 0  (needs an interior point with margin < 0 to exclude).
    Cost: 2^d interval evaluations per co-player (~2-4 us each)."""
    pat = np.ascontiguousarray(parse_pattern(pat), dtype=np.int8)
    ql = None if qlist is None else np.ascontiguousarray(qlist, dtype=np.int64)
    n = NQ if ql is None else len(ql)
    E = np.zeros(n, dtype=np.int32); cnt = np.zeros((n, 7), dtype=np.int32); mn = np.zeros(n); mx = np.zeros(n)
    _L.corner_scan(pat.ctypes.data_as(_b), int(sign), n, None if ql is None else _lp(ql), _ip(E), _ip(cnt),
                   _dp(mn), _dp(mx), _nt(nthreads))
    bad = (cnt[:, C_LOW] > 0) | (cnt[:, C_UNR] > 0)
    return dict(q=np.arange(NQ) if ql is None else ql, E=E, cnt=cnt, mmin=mn, mmax=mx,
                pos_ok=(cnt[:, C_NEG] == 0) & ~bad, neg_ok=(cnt[:, C_POS] == 0) & ~bad)

def corner_detail(pat, sign, g):
    """per-corner data for ONE co-player: dict(E (open stratum), Ec[2^d], cls[2^d] (0 structural
    zero, 1 > 0, -1 < 0, 3 zero by integrality, 2 unresolved, 9 if E_c != E), mlo, mhi, Dhi);
    corner m has free coordinate t = bit t of m."""
    pat = np.ascontiguousarray(parse_pattern(pat), dtype=np.int8); d = int((pat == 2).sum()); n = 1 << d
    Ec = np.zeros(n, dtype=np.int32); cls = np.zeros(n, dtype=np.int32); lo = np.zeros(n); hi = np.zeros(n); Dh = np.zeros(n)
    E = _L.corner_detail(pat.ctypes.data_as(_b), int(sign), int(g), _ip(Ec), _ip(cls), _dp(lo), _dp(hi), _dp(Dh))
    return dict(E=E, Ec=Ec, cls=cls, mlo=lo, mhi=hi, Dhi=Dh)

def binary_pair_exact(pg, qg):
    """EXACT limit of a binary pair from the interval engine alone: returns (A, D) with A[k] the
    integer number of minimum-order spanning trees rooted at k and D = sum(A) = D*, so that
    w_k = A[k]/D exactly.  Each integer is the unique one inside its enclosure (checked)."""
    p = g2p(pg); r = i_eval(p, p, qg)
    D = math.ceil(r['Dlo'])
    if not (D <= r['Dhi'] < D + 1): raise ArithmeticError("D* enclosure does not isolate an integer")
    A = []
    for k in range(16):
        if r['whi'][k] == 0: A.append(0); continue
        lo = float(np.nextafter(r['wlo'][k]*r['Dlo'], -np.inf)); hi = float(np.nextafter(r['whi'][k]*r['Dhi'], np.inf))
        a = math.ceil(lo)
        if not (a <= hi < a + 1): raise ArithmeticError("A*_k enclosure does not isolate an integer")
        A.append(a)
    if sum(A) != D: raise ArithmeticError("tree counts do not add up")
    return A, D

_L.i_scanW.argtypes = [_d, _d, ctypes.c_long, _l, _d, _d, _i, _i, ctypes.c_int]; _L.i_scanW.restype = ctypes.c_int
def i_scanW(plo, phi=None, qlist=None, nthreads=1):
    """interval enclosures of all four outcome frequencies against each co-player (default all
    65536): dict(q, Wlo (n x 4), Whi (n x 4), E, supp).  Valid for every focal point in the box."""
    plo, phi = _box(plo, phi)
    ql = None if qlist is None else np.ascontiguousarray(qlist, dtype=np.int64)
    n = NQ if ql is None else len(ql)
    Wlo = np.zeros((n, 4)); Whi = np.zeros((n, 4)); E = np.zeros(n, dtype=np.int32); s = np.zeros(n, dtype=np.int32)
    r = _L.i_scanW(_dp(plo), _dp(phi), n, None if ql is None else _lp(ql), _dp(Wlo), _dp(Whi), _ip(E), _ip(s), _nt(nthreads))
    if r: raise ValueError("invalid interval input")
    return dict(q=np.arange(NQ) if ql is None else ql, Wlo=Wlo, Whi=Whi, E=E, supp=s)

def outcome_present(supp, o):
    """exact: does the limit support contain a state whose recent outcome is o (0 CC,1 CD,2 DC,3 DD)?
    Vectorised over an int array of support masks."""
    m = 0xF << (4*o)
    return (np.asarray(supp) & m) != 0

def vertex_scan(vlo, vhi, sign=+1, qlist=None, nthreads=1, cap=1_000_000):
    """interval margins at nv box vertices (vlo, vhi: nv x 16 enclosures) against co-players.
    Returns dict: q, cnt (n x 4: #vertices certified >0, <0, structural zero, unresolved), mmin,
    unres (k x 2 array of (index into q, vertex index) with an enclosure containing 0), nunres."""
    vlo = _f64(vlo); vhi = _f64(vhi); nv = vlo.shape[0]
    ql = None if qlist is None else np.ascontiguousarray(qlist, dtype=np.int64)
    n = NQ if ql is None else len(ql)
    cnt = np.zeros((n, 4), dtype=np.int32); mn = np.zeros(n); un = np.zeros(2*cap, dtype=np.int64); nun = ctypes.c_long()
    r = _L.vertex_scan(_dp(vlo), _dp(vhi), nv, int(sign), n, None if ql is None else _lp(ql), _ip(cnt), _dp(mn),
                       _lp(un), cap, ctypes.byref(nun), _nt(nthreads))
    if r: raise ValueError("invalid vertex enclosure")
    k = min(nun.value, cap)
    return dict(q=np.arange(NQ) if ql is None else ql, cnt=cnt, mmin=mn, unres=un[:2*k].reshape(k, 2), nunres=nun.value)
