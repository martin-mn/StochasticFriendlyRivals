"""ctypes driver for ../cannon/libm2f.dylib -- the exact eps->0 continuum engine."""
import ctypes, os, numpy as np
_L = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cannon", "libm2f.dylib"))
_d = ctypes.POINTER(ctypes.c_double); _l = ctypes.POINTER(ctypes.c_long)
_L.c_freqs.argtypes=[_d,_d,_d]; _L.c_self.argtypes=[_d,_d]
_L.c_freqs16.argtypes=[_d,_d,_d]
_L.c_rival.argtypes=[_d,ctypes.c_int,_l,ctypes.c_double]; _L.c_rival.restype=ctypes.c_double
_L.c_batch.argtypes=[_d,ctypes.c_long,_d,ctypes.c_double]
_L.c_batch_wit.argtypes=[_d,ctypes.c_long,_l,ctypes.c_long,_d]

def _p(a): return np.ascontiguousarray(a,dtype=np.float64).ctypes.data_as(_d)

def freqs(p,q):
    W=np.zeros(4); _L.c_freqs(_p(p),_p(q),_p(W)); return W
def freqs16(p,q):
    w=np.zeros(16); _L.c_freqs16(_p(p),_p(q),_p(w)); return w
def selfplay(p):
    W=np.zeros(4); _L.c_self(_p(p),_p(W)); return W
def rival_margin(p,sign,earlyquit=-1e18):
    a=ctypes.c_long(0)
    v=_L.c_rival(_p(p),sign,ctypes.byref(a),earlyquit); return v,a.value
def batch(P,earlyquit=-1e-9):
    P=np.ascontiguousarray(P,dtype=np.float64); n=P.shape[0]
    out=np.zeros(8*n)
    _L.c_batch(P.ctypes.data_as(_d),n,out.ctypes.data_as(_d),earlyquit)
    return out.reshape(n,8)
def batch_wit(P,wit):
    P=np.ascontiguousarray(P,dtype=np.float64); n=P.shape[0]
    w=np.ascontiguousarray(wit,dtype=np.int64); out=np.zeros(2*n)
    _L.c_batch_wit(P.ctypes.data_as(_d),n,w.ctypes.data_as(_l),len(w),out.ctypes.data_as(_d))
    return out.reshape(n,2)
def g2p(g): return np.array([(g>>j)&1 for j in range(16)],dtype=float)
def p2g(p):
    g=0
    for j in range(16):
        if p[j]>0.5: g|=1<<j
    return g
