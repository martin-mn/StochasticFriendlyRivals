#!/usr/bin/env python3
"""Decisive check: does the exact MDP occupation-measure LP reproduce the binary memory-two
rival sets of ../../PartnersRivals/check (rivP = rivM = 2640, computed there by an entirely
different route -- eps->0 GTH leading-order arithmetic against all 65536 memory-two co-players)?"""
import os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m2lib

MASKS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "..", "PartnersRivals", "Figures", "m2_masks.npz")

def job(args):
    lo, hi = args
    outP = np.zeros(hi - lo, dtype=bool); outM = np.zeros(hi - lo, dtype=bool)
    valP = np.zeros(hi - lo); valM = np.zeros(hi - lo)
    for i, g in enumerate(range(lo, hi)):
        p = m2lib.genome_to_p(g)
        vP, _ = m2lib.occupation_lp(p, +1)
        vM, _ = m2lib.occupation_lp(p, -1)
        valP[i] = vP; valM[i] = vM
        outP[i] = vP >= -1e-9; outM[i] = vM >= -1e-9
    return lo, outP, outM, valP, valM

if __name__ == "__main__":
    N = 65536
    chunks = [(i, min(i + 512, N)) for i in range(0, N, 512)]
    rP = np.zeros(N, dtype=bool); rM = np.zeros(N, dtype=bool)
    vP = np.zeros(N); vM = np.zeros(N)
    with Pool(12) as pool:
        for lo, oP, oM, aP, aM in pool.imap_unordered(job, chunks):
            rP[lo:lo+len(oP)] = oP; rM[lo:lo+len(oM)] = oM
            vP[lo:lo+len(aP)] = aP; vM[lo:lo+len(aM)] = aM
    M = np.load(MASKS)
    ref_P, ref_M = M["rivP"], M["rivM"]
    print("LP rivals T>S : %d   reference %d" % (rP.sum(), ref_P.sum()))
    print("LP rivals T<S : %d   reference %d" % (rM.sum(), ref_M.sum()))
    print("exact set match T>S:", bool(np.array_equal(rP, ref_P)))
    print("exact set match T<S:", bool(np.array_equal(rM, ref_M)))
    if not np.array_equal(rP, ref_P):
        d = np.where(rP != ref_P)[0]
        print("  %d disagreements, first 20: %s" % (len(d), d[:20].tolist()))
        print("  LP values there:", vP[d[:20]])
    # margin structure
    nz = vP[vP < -1e-12]
    print("worst-case margin: min %.6g, max negative %.6g, #zero %d"
          % (vP.min(), nz.max() if len(nz) else float('nan'), int((np.abs(vP) < 1e-12).sum())))
    np.savez_compressed(os.path.join(os.path.dirname(os.path.abspath(__file__)), "binary_lp.npz"),
                        rivP=rP, rivM=rM, valP=vP, valM=vM)
