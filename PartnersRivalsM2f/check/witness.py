#!/usr/bin/env python3
"""Find a small set of binary memory-two co-players that certifies non-rivalry.

Greedy set cover: repeatedly add the co-player beating the most still-uncovered non-rivals.
Candidate pool = the observed first-beaters plus a random sample, so the cover is exact on
the 65536 binary strategies (verified at the end) and is used only as a *screen* in the
continuum: a strategy surviving the screen is still put through the full 65536 scan."""
import os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m2f

M = np.load("/Users/martin/Documents/CL1/PartnersRivals/Figures/m2_masks.npz")
rivP = M["rivP"]
D = np.load("binary_margin.npz"); argP = D["argP"]

nonriv = np.where(~rivP)[0]
pool_cand = np.unique(np.concatenate([argP[nonriv], np.arange(0, 65536, 257)]))
print("non-rivals %d, candidate co-players %d" % (len(nonriv), len(pool_cand)))

P = np.stack([m2f.g2p(g) for g in nonriv])

def beats(t):
    """which non-rivals does co-player t beat (W_DC - W_CD < 0)?"""
    o = m2f.batch_wit(P, np.array([t], dtype=np.int64))
    return o[:, 0] < -1e-9

def job(ts):
    return [(t, beats(t)) for t in ts]

if __name__ == "__main__":
    chunks = [pool_cand[i::12] for i in range(12)]
    cover = {}
    with Pool(12) as pool:
        for res in pool.imap_unordered(job, chunks):
            for t, b in res:
                cover[int(t)] = b
    unc = np.ones(len(nonriv), dtype=bool)
    chosen = []
    while unc.any():
        best, bn = None, -1
        for t, b in cover.items():
            n = int((b & unc).sum())
            if n > bn: bn, best = n, t
        if bn <= 0:
            print("STUCK: %d uncovered" % unc.sum()); break
        chosen.append(best); unc &= ~cover[best]
        print("  + co-player %5d covers %6d, remaining %6d" % (best, bn, unc.sum()))
    print("witness set (%d):" % len(chosen), chosen)
    np.save("witnesses.npy", np.array(chosen, dtype=np.int64))
