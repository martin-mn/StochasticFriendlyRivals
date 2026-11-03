"""Independent attack: is there ANY wedge-W friendly-rival stratum of dimension >= 8 ?

Stage 1  enumerate EVERY dimension-NF support stratum -- C(16,NF)*2^(16-NF) of them,
         3 294 720 at NF=8 -- and classify it by two independent interior draws.
         Keep the class-C ones (self-play W_CC = 1, i.e. efficient below the switch line).
Stage 2  give each class-C stratum many random interior draws and screen them against the
         19 refuting co-players.  The screen has NO false negatives: a genuine rival has
         margin >= 0 against every co-player, those 19 included.
Stage 3  every screen survivor gets the exact full 65536-co-player scan.
"""
import sys, time, itertools, numpy as np, multiprocessing as mp
sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivalsM2f/check")
B = "/Users/martin/Documents/CL1/PartnersRivalsM2f/check/working/verify_task1/"
NF    = int(sys.argv[1]) if len(sys.argv) > 1 else 8
DRAWS = int(sys.argv[2]) if len(sys.argv) > 2 else 512
COMB  = list(itertools.combinations(range(16), NF))
NMASK = 1 << (16-NF)

def stage1(a):
    lo, hi = a
    import m2f
    rng = np.random.default_rng(20260909 + lo)
    keep = []; nonconst = 0
    for ci in range(lo, hi):
        free  = list(COMB[ci])
        fixed = [j for j in range(16) if j not in free]
        R = rng.random((NMASK, 2, NF)) * 0.998 + 0.001
        for mask in range(NMASK):
            p = np.empty(16)
            for k, j in enumerate(fixed):
                p[j] = float((mask >> k) & 1)
            p[free] = R[mask, 0]; W0 = m2f.selfplay(p)
            p[free] = R[mask, 1]; W1 = m2f.selfplay(p)
            if abs(W0[0]-W1[0]) > 1e-9: nonconst += 1
            if abs(W0[0]-1.0) < 1e-9 and abs(W1[0]-1.0) < 1e-9:
                keep.append(ci*NMASK + mask)
    return keep, nonconst

def stage23(a):
    lo, hi, ids, sign = a
    import m2f
    rng = np.random.default_rng(777000 + lo)
    wit = np.load(B + "../../witnesses.npy") if False else np.load(
        "/Users/martin/Documents/CL1/PartnersRivalsM2f/check/witnesses.npy")
    surv = []; nsurv = 0; nconf = 0
    for t in range(lo, hi):
        sid  = int(ids[t]); ci, mask = divmod(sid, NMASK)
        free = list(COMB[ci]); fixed = [j for j in range(16) if j not in free]
        P = np.empty((DRAWS, 16))
        for k, j in enumerate(fixed): P[:, j] = float((mask >> k) & 1)
        P[:, free] = rng.random((DRAWS, NF)) * 0.998 + 0.001
        w = m2f.batch_wit(P, wit)              # [:,0] = min over witnesses of (W_DC-W_CD)
        col = 0 if sign > 0 else 1
        cand = np.flatnonzero(w[:, col] >= -1e-12)
        nsurv += len(cand)
        for i in cand:
            m, arg = m2f.rival_margin(P[i], sign, -1e-9)
            if m >= -1e-12:
                nconf += 1
                surv.append((sid, float(m), int(arg)))
    return surv, nsurv, nconf

if __name__ == "__main__":
    NPROC = 12
    n = len(COMB)
    print("dimension %d : %d free-sets x %d pin-patterns = %d strata"
          % (NF, n, NMASK, n*NMASK), flush=True)
    t0 = time.time()
    bnds = np.linspace(0, n, NPROC+1).astype(int)
    with mp.Pool(NPROC) as pool:
        res = pool.map(stage1, [(bnds[k], bnds[k+1]) for k in range(NPROC)])
    ids = np.array(sorted(x for r in res for x in r[0]), dtype=np.int64)
    nonconst = sum(r[1] for r in res)
    print("  stage 1  %.1f s :  class-C (W_CC=1) dim-%d strata = %d ;"
          " strata whose class differed between the two draws = %d"
          % (time.time()-t0, NF, len(ids), nonconst), flush=True)
    np.save(B + "dim%d_C.npy" % NF, ids)

    t0 = time.time()
    m = len(ids); bnds = np.linspace(0, m, NPROC+1).astype(int)
    with mp.Pool(NPROC) as pool:
        res = pool.map(stage23, [(bnds[k], bnds[k+1], ids, +1) for k in range(NPROC)])
    surv  = [x for r in res for x in r[0]]
    nsurv = sum(r[1] for r in res); nconf = sum(r[2] for r in res)
    print("  stage 2  %.1f s :  %d class-C strata x %d draws = %d points;"
          " %d survived the 19-witness screen" % (time.time()-t0, m, DRAWS, m*DRAWS, nsurv), flush=True)
    print("  stage 3  :  %d of those confirmed rivals (T>S) against ALL 65536 co-players" % nconf)
    if surv:
        st = sorted(set(s[0] for s in surv))
        print("  ==> WEDGE-W FRIENDLY RIVALS OF DIMENSION %d FOUND: %d strata" % (NF, len(st)))
        for sid in st[:20]:
            ci, mask = divmod(sid, NMASK)
            free = list(COMB[ci]); fixed = [j for j in range(16) if j not in free]
            pat = ["*"]*16
            for k, j in enumerate(fixed): pat[j] = "01"[(mask >> k) & 1]
            print("      ", "".join(pat))
    else:
        print("  ==> NO wedge-W friendly rival of dimension %d exists in %d confirmed draws."
              % (NF, m*DRAWS))
