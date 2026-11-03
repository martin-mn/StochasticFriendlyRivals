"""Attack, take 2 (vectorised).  Is there ANY wedge-W friendly-rival stratum of dim >= 8 ?

Stage 1 (already done, dim8_C.npy): every dim-8 stratum classified from 2 interior draws;
        491127 are class C (self-play W_CC = 1).  ALL of them have p0 pinned to 1.
Stage 2 screen every class-C dim-8 stratum with many random interior draws against a witness
        set = the project's 19 refuting co-players + ALLD(0) + ALLC(65535).  The screen has NO
        false negatives: a genuine rival has margin >= 0 against EVERY co-player.
        Strata with p15 pinned to 0 -- the only ones Riv+ subset {p15=0} allows -- get NDEEP
        draws; the rest get NSHALLOW as a belt-and-braces check of that lemma.
Stage 3 every survivor gets the exact full 65536-co-player scan.
"""
import sys, time, itertools, numpy as np, multiprocessing as mp
sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivalsM2f/check")
B  = "/Users/martin/Documents/CL1/PartnersRivalsM2f/check/working/verify_task1/"
NF = int(sys.argv[1]) if len(sys.argv) > 1 else 8
NDEEP    = int(sys.argv[2]) if len(sys.argv) > 2 else 512
NSHALLOW = int(sys.argv[3]) if len(sys.argv) > 3 else 64
COMB  = np.array(list(itertools.combinations(range(16), NF)), dtype=np.int64)
NMASK = 1 << (16-NF)
WIT   = np.concatenate([np.load("/Users/martin/Documents/CL1/PartnersRivalsM2f/check/witnesses.npy"),
                        np.array([0, 65535], dtype=np.int64)])   # + ALLD + ALLC

def digits(ids):
    """(n,16) array of base-3 digits (0 pinned to 0, 1 free, 2 pinned to 1) for stratum ids."""
    ids = np.asarray(ids, dtype=np.int64)
    ci, mask = ids // NMASK, ids % NMASK
    D = np.zeros((len(ids), 16), dtype=np.int8)
    free = COMB[ci]                                    # (n,NF)
    np.put_along_axis(D, free, 1, axis=1)
    for i in range(len(ids)):
        fixed = [j for j in range(16) if D[i, j] != 1]
        for k, j in enumerate(fixed):
            D[i, j] = 2 * ((mask[i] >> k) & 1)
    return D

def run(a):
    lo, hi, ids, nd, sign, seed = a
    import m2f
    rng = np.random.default_rng(seed)
    G = 256
    surv = []; npts = 0; nsurv = 0
    for s in range(lo, hi, G):
        blk = ids[s:min(s+G, hi)]
        D = digits(blk); n = len(blk)
        P = np.repeat(D.astype(np.float64), nd, axis=0)      # (n*nd,16)
        fr = (P == 1)
        P[P == 2] = 1.0
        P[fr] = rng.random(int(fr.sum())) * 0.998 + 0.001
        w = m2f.batch_wit(P, WIT)
        npts += len(P)
        col = 0 if sign > 0 else 1
        for i in np.flatnonzero(w[:, col] >= -1e-12):
            nsurv += 1
            m, arg = m2f.rival_margin(P[i], sign, -1e-9)      # exact 65536 scan
            if m >= -1e-12:
                surv.append((int(blk[i // nd]), float(m), int(arg)))
    return surv, npts, nsurv

if __name__ == "__main__":
    NPROC = 12
    ids = np.load(B + "dim%d_C.npy" % NF)
    D   = digits(ids)
    print("class-C dim-%d strata: %d ; p0 digit histogram %s ; p15 digit histogram %s"
          % (NF, len(ids), np.bincount(D[:, 0], minlength=3).tolist(),
             np.bincount(D[:, 15], minlength=3).tolist()), flush=True)
    deep    = ids[D[:, 15] == 0]
    shallow = ids[D[:, 15] != 0]
    print("  p15 pinned to 0 : %d  (deep, %d draws)   others : %d  (shallow, %d draws)"
          % (len(deep), NDEEP, len(shallow), NSHALLOW), flush=True)
    for tag, sub, nd in [("DEEP", deep, NDEEP), ("SHALLOW", shallow, NSHALLOW)]:
        t0 = time.time(); m = len(sub)
        bnds = np.linspace(0, m, NPROC+1).astype(int)
        with mp.Pool(NPROC) as pool:
            res = pool.map(run, [(bnds[k], bnds[k+1], sub, nd, +1, 424242+k) for k in range(NPROC)])
        sv = [x for r in res for x in r[0]]
        print("  %-7s %6d strata x %4d draws = %9d points ; %d survived the %d-witness screen ;"
              " %d CONFIRMED rivals (T>S) against all 65536   [%.1f s]"
              % (tag, m, nd, sum(r[1] for r in res), sum(r[2] for r in res), len(WIT),
                 len(sv), time.time()-t0), flush=True)
        if sv:
            for sid in sorted(set(x[0] for x in sv))[:20]:
                d = digits([sid])[0]
                print("      FOUND dim-%d wedge-W stratum: %s" % (NF, "".join("0*1"[x] for x in d)))
    print("  ==> no wedge-W friendly rival of dimension %d found." % NF)
