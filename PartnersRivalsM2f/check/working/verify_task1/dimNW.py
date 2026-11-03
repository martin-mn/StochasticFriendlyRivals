"""Same attack at dimension NF (stage 1 enumeration + stage 2 screen + stage 3 exact scan)."""
import sys, time, itertools, numpy as np, multiprocessing as mp
sys.path.insert(0, "/Users/martin/Documents/CL1/PartnersRivalsM2f/check")
B  = "/Users/martin/Documents/CL1/PartnersRivalsM2f/check/working/verify_task1/"
NF = int(sys.argv[1]); NDEEP = int(sys.argv[2]); NSHALLOW = int(sys.argv[3])
COMB  = np.array(list(itertools.combinations(range(16), NF)), dtype=np.int64)
NMASK = 1 << (16-NF)
WIT   = np.concatenate([np.load("/Users/martin/Documents/CL1/PartnersRivalsM2f/check/witnesses.npy"),
                        np.array([0, 65535], dtype=np.int64)])

def digits(ids):
    ids = np.asarray(ids, dtype=np.int64)
    ci, mask = ids // NMASK, ids % NMASK
    D = np.zeros((len(ids), 16), dtype=np.int8)
    np.put_along_axis(D, COMB[ci], 1, axis=1)
    for i in range(len(ids)):
        fixed = [j for j in range(16) if D[i, j] != 1]
        for k, j in enumerate(fixed):
            D[i, j] = 2 * ((mask[i] >> k) & 1)
    return D

def stage1(a):
    lo, hi = a
    import m2f
    rng = np.random.default_rng(555 + lo)
    keep = []; flip = 0
    for ci in range(lo, hi):
        free = list(COMB[ci]); fixed = [j for j in range(16) if j not in free]
        R = rng.random((NMASK, 2, NF)) * 0.998 + 0.001
        for mask in range(NMASK):
            p = np.empty(16)
            for k, j in enumerate(fixed): p[j] = float((mask >> k) & 1)
            p[free] = R[mask,0]; a0 = abs(m2f.selfplay(p)[0]-1.0) < 1e-9
            p[free] = R[mask,1]; a1 = abs(m2f.selfplay(p)[0]-1.0) < 1e-9
            if a0 != a1: flip += 1          # true test of stratum-constancy of the EFFICIENCY INDICATOR
            if a0 and a1: keep.append(ci*NMASK + mask)
    return keep, flip

def run(a):
    lo, hi, ids, nd, seed = a
    import m2f
    rng = np.random.default_rng(seed); G = 256
    surv = []; npts = 0; nsurv = 0
    for s in range(lo, hi, G):
        blk = ids[s:min(s+G, hi)]
        D = digits(blk)
        P = np.repeat(D.astype(np.float64), nd, axis=0)
        fr = (P == 1); P[P == 2] = 1.0
        P[fr] = rng.random(int(fr.sum())) * 0.998 + 0.001
        w = m2f.batch_wit(P, WIT); npts += len(P)
        for i in np.flatnonzero(w[:, 0] >= -1e-12):
            nsurv += 1
            m, arg = m2f.rival_margin(P[i], +1, -1e-9)
            if m >= -1e-12: surv.append((int(blk[i//nd]), float(m), int(arg)))
    return surv, npts, nsurv

if __name__ == "__main__":
    NPROC = 12; n = len(COMB)
    print("=== dimension %d : %d free-sets x %d pin-patterns = %d strata ==="
          % (NF, n, NMASK, n*NMASK), flush=True)
    t0 = time.time(); bnds = np.linspace(0, n, NPROC+1).astype(int)
    with mp.Pool(NPROC) as pool:
        res = pool.map(stage1, [(bnds[k], bnds[k+1]) for k in range(NPROC)])
    ids = np.array(sorted(x for r in res for x in r[0]), dtype=np.int64)
    print("  stage 1 %.1f s : class-C (W_CC=1) strata = %d ; efficiency INDICATOR disagreed "
          "between the two draws on %d strata" % (time.time()-t0, len(ids), sum(r[1] for r in res)), flush=True)
    if len(ids) == 0:
        print("  ==> no class-C stratum of dimension %d at all." % NF); sys.exit()
    D = digits(ids)
    print("  p0 digit histogram %s (2 = pinned to 1) ; p15 digit histogram %s"
          % (np.bincount(D[:,0], minlength=3).tolist(), np.bincount(D[:,15], minlength=3).tolist()), flush=True)
    deep = ids[D[:,15]==0]; shal = ids[D[:,15]!=0]
    tot = 0
    for tag, sub, nd in [("DEEP", deep, NDEEP), ("SHALLOW", shal, NSHALLOW)]:
        if len(sub)==0:
            print("  %-7s none" % tag, flush=True); continue
        t0 = time.time(); m = len(sub); bn = np.linspace(0, m, NPROC+1).astype(int)
        with mp.Pool(NPROC) as pool:
            res = pool.map(run, [(bn[k], bn[k+1], sub, nd, 9090+k) for k in range(NPROC)])
        sv = [x for r in res for x in r[0]]; tot += len(sv)
        print("  %-7s %6d strata x %4d draws = %9d pts ; %6d survived the %d-witness screen ; "
              "%d CONFIRMED rivals (T>S) vs all 65536  [%.1f s]"
              % (tag, m, nd, sum(r[1] for r in res), sum(r[2] for r in res), len(WIT), len(sv),
                 time.time()-t0), flush=True)
        for sid in sorted(set(x[0] for x in sv))[:10]:
            print("      FOUND: %s" % "".join("0*1"[x] for x in digits([sid])[0]))
    print("  ==> dimension %d : %d wedge-W friendly rivals found." % (NF, tot), flush=True)
