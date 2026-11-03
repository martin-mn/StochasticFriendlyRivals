"""Loaders for the production strategy dumps of PartnersRivalsM2f/cannon/wfc (read-only copies in ../data)."""
import os, numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); KIT = os.path.dirname(HERE)
REC = np.dtype([('irep', '<i4'), ('count', '<i4'), ('atom', '<i4'), ('p', '<f8', (16,))])
assert REC.itemsize == 140
def games(path=os.path.join(KIT, 'games.dat')):
    """list over tasks 1..512 (index task-1) of dict(gid, ipt, cs, ct, tag, wedge, below, sign)"""
    G = []
    for line in open(path):
        if line.startswith('#'): continue
        t = line.split()
        if len(t) < 8: continue
        cs, ct = float(t[3]), float(t[4]); assert float(t[2]) == 1.0
        below = Fraction(cs) + Fraction(ct) < 2; sign = +1 if Fraction(ct) > Fraction(cs) else -1
        assert Fraction(cs) + Fraction(ct) != 2 and Fraction(cs) != Fraction(ct)
        wedge = {(True, 1): 'W', (True, -1): 'S', (False, 1): 'N', (False, -1): 'E'}[(bool(below), sign)]
        # production's own double tests (wfc.c classify): below = UU+VV<1, sign = UU-VV<1 ? +1 : -1
        UU = cs; VV = ct - 1.0
        assert (UU + VV < 1.0) == below and ((UU - VV < 1.0) == (sign == 1))
        G.append(dict(gid=int(t[0]), ipt=int(t[7][3:6]), cs=cs, ct=ct, tag=t[7], wedge=wedge, below=bool(below), sign=sign))
    return G
def dump(regime, task):
    return np.fromfile(os.path.join(KIT, 'data', regime, 'S_1o3_%d_1.bin' % task), dtype=REC)
