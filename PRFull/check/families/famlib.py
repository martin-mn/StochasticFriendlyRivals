"""Shared helpers for the exact family certificates (cert_pt.py, cert_ab.py) and the witness archive.

Everything here is exact rational arithmetic (cert/exact.py, gmpy2.mpq inside, Fractions out) except
where a function says it calls the float or interval engine.

The objects computed per (focal point x, binary co-player q):
    E     the least spanning-tree eps-order of the pair chain (integer),
    D     D*(x), the eps^E coefficient of the total tree sum  sum_k A_k(eps),
    W     the exact limit frequencies (w_CC, w_CD, w_DC, w_DD), focal view,
    m     the competitive margin sign*(w_DC - w_CD),
    N     D* * m  = sign*(A*_DC - A*_CD), which by the structural lemma is a MULTILINEAR polynomial in
          the free coordinates on the open stratum (for fixed q), with E constant there,
    supp  bitmask of the limit support.
Local-compute rule: pools are capped at 3 processes (cert.engine.MAX_THREADS).
"""
import os, sys, time
sys.dont_write_bytecode = True
_CHECK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _CHECK not in sys.path: sys.path.insert(0, _CHECK)
from fractions import Fraction as F
from multiprocessing import Pool
from cert import engine as en, exact as ex

NQ = 65536

def _fr(x): return F(int(x.numerator), int(x.denominator))

def point_data(p, g, sign):
    """exact (E, D*, N, W, supp) of focal point p against binary co-player g"""
    r = ex.lead_exact(p, ex.g2p(int(g)), as_fraction=False)
    W = r['W']; D = r['D']
    m = sign*(W[2] - W[1])
    return dict(E=int(r['E']), D=_fr(D), N=_fr(D*m), m=_fr(m), W=tuple(_fr(w) for w in W), supp=int(r['supp']))

def _chunk(args):
    points, sign, qs = args
    out = []
    for g in qs:
        out.append((int(g), [point_data(p, g, sign) for p in points]))
    return out

def exact_scan(points, sign, qlist=None, nproc=None, chunk=512):
    """exact point_data for every co-player in qlist (default all 65536) at every point of `points`
    (a list of 16-vectors of rationals).  Returns dict g -> [point_data per point]."""
    qlist = list(range(NQ)) if qlist is None else [int(q) for q in qlist]
    pts = [[F(x) for x in p] for p in points]
    nproc = max(1, min(nproc or en.MAX_THREADS, en.MAX_THREADS))
    jobs = [(pts, sign, qlist[i:i + chunk]) for i in range(0, len(qlist), chunk)]
    res = {}
    if nproc == 1 or len(jobs) == 1:
        for j in jobs:
            for g, v in _chunk(j): res[g] = v
    else:
        with Pool(nproc) as pool:
            for part in pool.imap_unordered(_chunk, jobs):
                for g, v in part: res[g] = v
    return res

def fstr(x):
    """a rational as a string 'n/d' (or 'n')"""
    x = F(x); return str(x)

def vec_str(v): return [fstr(x) for x in v]

def unf_fav(sign):
    """(favourable, unfavourable) outcome indices for the margin sign*(w_DC - w_CD)"""
    return (2, 1) if sign > 0 else (1, 2)

def structural_zero(supp):
    """w_CD = w_DC = 0 by exponents alone (neither outcome CD nor DC in the limit support):
    by L1 this holds on the whole open stratum"""
    return not (en.outcome_present(supp, 1) or en.outcome_present(supp, 2))

def timer():
    t0 = time.time()
    return lambda: time.time() - t0
