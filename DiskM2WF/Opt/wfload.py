"""
Load a dw1 (WF) run and present it in exactly the shape DiskM2IN's analysis
expects, so that maps.py, squaremaps.py and squarereport.py port with the
loader swapped and nothing else.

Two sources, same output:

  * a packed pair `dw1_<e>.{idx,bin}` written on Cannon by
    Cannon/dw1/pack.py  -- the normal path;
  * a raw `run/` directory of `<isl>`, `w<isl>`, `h<isl>` files -- for the
    probe, or for anything that was not packed.

`load_run(path)` returns {gid: game dict} with the DiskM2IN field names:

    gid tag u v cR cS cT emax pi eff neff maxpi top toppi ps

and the WF-specific extras

    pay    mean raw per-round payoff over the cell's replicates
    efsd   s.d. of pay/Emax over the replicates -- the run-to-run spread,
           which is the honest error bar on `eff`
    nrep   replicate lines found (10 if the cell is complete)
    ndist  strategies with a non-zero summed count
    ptot   the grand total count, nrep*(itend/2)*N

`dead`, `absorbing`, `it` and `L1` are set so that the ported code's branches
behave: the WF process has no dead states and no iterative solve, so
dead=False, absorbing=False, it=0, L1=0.0.

WHAT `pi` IS, AND WHAT IT IS NOT.  Here `pi` is the ABUNDANCE vector: the
(time x individual) frequency of each strategy over the sampled second half of
the runs, so it includes the polymorphic mutant cloud.  DiskM2IN's `pi` is the
fraction of time the embedded chain sits at the monomorphic state.  They are
not the same quantity; any figure that puts them side by side has to say so.
"""
import os, glob
import numpy as np

NC = 65536


def ne_uv(cR, cS, cT, udisp=None, vdisp=None, where=""):
    """NE's coordinates of the game (cR, cS, cT, 0).

    NE fixes the representative by R - P = 1, so (cR,cS,cT,0) is the same game
    as (1, cS/cR, cT/cR, 0) and therefore

        u = cS/cR        v = cT/cR - 1

    NOT u = cS, v = cT - 1.  Those coincide only when cR = 1, which is true of
    every disk and square row but FALSE of the 19 DONwf3 rows, where
    cR = 1-c: at c = 0.65 the game is (0.35, -0.65, 1, 0) and its NE
    coordinates are (-1.857, +1.857), not (-0.65, 0).  Getting this wrong puts
    the whole donation ray in the wrong place on the disk.

    Cross-checked against the u/v display columns the kit copies out of
    games.dat, which is what makes it a check rather than an assertion.
    """
    if not cR > 0.0:
        raise ValueError("%s: cR = %r, not positive" % (where, cR))
    u = cS / cR
    v = cT / cR - 1.0
    if udisp is not None:
        d = max(abs(u - udisp), abs(v - vdisp))
        if d > 1e-6:
            raise ValueError("%s: NE coords from (cR,cS,cT) are (%.9f,%.9f) "
                             "but games.dat says (%.9f,%.9f), off by %.2e"
                             % (where, u, v, udisp, vdisp, d))
    return u, v


def _game(gid, isl, tag, u, v, cR, cS, cT, emax, pay, ef, efsd,
          ndist, nrep, ptot, ps, p):
    p = np.asarray(p, dtype=np.float64)
    s = p.sum()
    if s > 0:
        p = p / s
    top = np.argsort(p)[::-1][:10]
    return dict(gid=gid, isl=isl, tag=tag, u=u, v=v, cR=cR, cS=cS, cT=cT,
                emax=emax, pi=p, eff=ef, pay=pay, efsd=efsd,
                neff=float(1.0 / np.square(p).sum()), maxpi=float(p.max()),
                top=[int(t) for t in top],
                toppi=[float(p[t]) for t in top],
                ps=list(ps), ndist=ndist, nrep=nrep, ptot=ptot,
                dead=False, absorbing=False, it=0, L1=0.0)


def _attach_win(g, win):
    """What the 10 independent replicates each led with.

    The abundance vector is summed over the replicates, so its argmax can be a
    strategy that led in only a few of them.  `nlead` is how many replicates
    made the cell's overall leader their own rank-1, and `runlead` is the full
    multiset.  This is a direct measurement of leader reproducibility, and it
    is the WF analogue of DiskM2IN's residual audit -- with the advantage of
    being a measurement rather than a bound.  (DiskM2IN §3g: its residual
    audit missed its two worst cells entirely.)
    """
    codes = [c for _, c, _, _ in win]
    g["runlead"] = codes
    g["runab"] = [a for _, _, a, _ in win]
    if codes:
        lead = g["top"][0]
        g["nlead"] = sum(1 for c in codes if c == lead)
        g["nlead_distinct"] = len(set(codes))
    else:
        g["nlead"] = -1
        g["nlead_distinct"] = -1
    return g


# ------------------------------------------------------------- packed pair
def load_packed(idxpath):
    binpath = idxpath[:-4] + ".bin"
    winpath = idxpath[:-4] + ".win"
    rows = [l.split() for l in open(idxpath) if not l.startswith("#")]
    blob = np.fromfile(binpath, dtype="<f4")
    if blob.size != len(rows) * NC:
        raise ValueError("%s holds %d floats, index has %d cells x %d"
                         % (binpath, blob.size, len(rows), NC))
    blob = blob.reshape(len(rows), NC)
    G = {}
    for k, r in enumerate(rows):
        gid, isl, tag = int(r[0]), int(r[1]), r[2]
        udisp, vdisp, cR, cS, cT, emax, pay, ef, efsd = map(float, r[3:12])
        u, v = ne_uv(cR, cS, cT, udisp, vdisp, "%s gid %d" % (idxpath, gid))
        ndist, nrep, ptot = int(r[12]), int(r[13]), int(r[14])
        ps = [float(x) for x in r[15:31]]
        G[gid] = _game(gid, isl, tag, u, v, cR, cS, cT, emax, pay, ef, efsd,
                       ndist, nrep, ptot, ps, blob[k])
    # The .win is NOT optional.  maps.audit reads nlead, and a cell with no
    # nlead gets frac = NaN, which its `frac == frac and frac < minfrac` test
    # short-circuits to False -- so a missing sidecar makes the leader
    # reproducibility instrument report every cell safe, from no data at all.
    # That is exactly the failure mode §8a introduced it to avoid.  pack.py
    # always writes it, so its absence means the copy down was incomplete.
    if not os.path.exists(winpath):
        raise IOError("%s is missing.  It carries the per-replicate rank-1 "
                      "strategies, which are the leader-reproducibility "
                      "instrument; without it maps.audit reports every cell "
                      "safe from no data.  Copy it down with the .idx/.bin."
                      % winpath)
    W = {}
    for ln in open(winpath):
        if ln.startswith("#"):
            continue
        f = ln.split()
        W.setdefault(int(f[0]), []).append(
            (int(f[1]), int(f[2]), float(f[3]), float(f[4])))
    for gid, g in G.items():
        _attach_win(g, W.get(gid, []))
    return G


# ---------------------------------------------------------------- raw run/
def _read_h(path):
    meta, code, cnt = {}, [], []
    for ln in open(path):
        if ln[0] == "#":
            t = ln.split()
            if ln.startswith("# cell:"):
                meta["isl"] = int(t[3]); meta["gid"] = int(t[5])
                meta["tag"] = t[6]
            elif ln.startswith("# game:"):
                meta["cR"], meta["cS"], meta["cT"] = map(float, t[5:8])
                meta["udisp"], meta["vdisp"] = float(t[10]), float(t[11])
                meta["emax"] = float(t[-1])
            elif ln.startswith("# distinct"):
                meta["ndist"] = int(t[2]); meta["ptot"] = int(t[4])
            continue
        a, b = ln.split()
        code.append(int(a)); cnt.append(int(b))
    # Emax > 0 always (Emax >= cR > 0 for every game in games.dat).  A header
    # of zeros is the signature of the first dw1 build, which wrote the
    # `# game:` line before assigning the payoffs; refuse that data rather
    # than silently reporting u = 0, v = -1, Emax = 0 for every cell.
    if not meta.get("emax", 0.0) > 0.0:
        raise ValueError("%s: '# game:' header has Emax = %r -- this h-file "
                         "predates the header fix, rerun it"
                         % (path, meta.get("emax")))
    # An h-file is written in one burst after the last replicate, so a task
    # still running -- or killed in that window -- leaves the header and no
    # counts, or some of them.  Both are refused here, on exactly the check
    # Cannon/dw1/pack.py makes: `# distinct ndist total ntotc` precedes the
    # counts and so survives truncation, and both invariants below hold for a
    # whole h-file whatever the replicate count.  Without this a partial
    # histogram is renormalised by its own partial sum and drawn as a normal
    # cell; an empty one gives neff = inf, maxpi = 0 and leader 0.
    if "ndist" not in meta:
        raise ValueError("%s: no '# distinct' line and so no histogram.  "
                         "sf.f writes the h-file's headers when the task "
                         "STARTS and the counts only after the last "
                         "replicate, so this is what a task that is still "
                         "running looks like -- wait for it, or rerun it."
                         % path)
    tot = sum(cnt)
    if len(cnt) != meta["ndist"] or tot != meta["ptot"]:
        raise ValueError("%s: TRUNCATED histogram -- %d count lines summing "
                         "to %d, but the header says ndist %d total %d.  The "
                         "task was still running, or died before its final "
                         "write; rerun it."
                         % (path, len(cnt), tot, meta["ndist"], meta["ptot"]))
    p = np.zeros(NC)
    p[np.array(code, dtype=np.int64)] = np.array(cnt, dtype=np.float64)
    return meta, p


def _read_w(path):
    """the w-file: the rank-1 strategy of each replicate.
       columns: isl eps gid irun rank code abund selfpay selfef w0..w3 bits"""
    out = []
    for ln in open(path):
        if ln[0] == "#":
            continue
        f = ln.split()
        if int(f[4]) != 1:
            continue
        out.append((int(f[3]), int(f[5]), float(f[6]), float(f[7])))
    return out


def load_raw(rundir, ie=None):
    """A raw run/ directory, optionally restricted to one error rate.

    `ie` is the error-rate index the task decode assigns, 1/2/3 for
    eps = 1e-2/1e-3/1e-4, and it is column 1 of every summary file.  It has to
    be part of the key, not just an option: run/ holds every array that has
    ever written into it, and the SAME GAME at two error rates has the same
    gid.  Keying on gid alone let the later isl silently overwrite the earlier
    one -- which is not hypothetical, it is the state of the probe's run/,
    where isl 37 (gid 13 at eps=1e-2) and isl 39 (gid 13 at eps=1e-4) are both
    present and only one survived the load, chosen by sort order.

    `load_raw(rundir, ie=3)` returns the eps=1e-4 layer.  With no `ie` the
    directory must hold exactly one error rate, or this raises rather than
    picking one.
    """
    G, seen = {}, {}
    for hp in sorted(glob.glob(os.path.join(rundir, "h[0-9]*")),
                     key=lambda q: int(os.path.basename(q)[1:])):
        isl = int(os.path.basename(hp)[1:])
        sp = os.path.join(rundir, str(isl))
        if not os.path.exists(sp):
            continue
        a = np.array([[float(x) for x in l.split()]
                      for l in open(sp) if l.strip()])
        if a.size == 0:
            continue
        cie = int(a[0, 1])
        if ie is not None and cie != ie:
            continue
        meta, p = _read_h(hp)
        if meta["gid"] in seen:
            raise ValueError(
                "%s: gid %d appears at isl %d (ie %d) and isl %d (ie %d).  "
                "This directory holds more than one error rate; pass "
                "ie=1/2/3 (eps = 1e-2/1e-3/1e-4) to choose one."
                % (rundir, meta["gid"], seen[meta["gid"]][0],
                   seen[meta["gid"]][1], isl, cie))
        seen[meta["gid"]] = (isl, cie)
        u, v = ne_uv(meta["cR"], meta["cS"], meta["cT"],
                     meta["udisp"], meta["vdisp"], hp)
        G[meta["gid"]] = _game(
            meta["gid"], isl, meta["tag"], u, v,
            meta["cR"], meta["cS"], meta["cT"], meta["emax"],
            float(a[:, 5].mean()), float(a[:, 6].mean()),
            float(a[:, 6].std(ddof=1)) if a.shape[0] > 1 else 0.0,
            meta.get("ndist", -1), a.shape[0], meta.get("ptot", -1),
            a[:, 7:23].mean(axis=0), p)
        wp = os.path.join(rundir, "w%d" % isl)
        if os.path.exists(wp):
            _attach_win(G[meta["gid"]], _read_w(wp))
    short = sorted(g["isl"] for g in G.values() if g["nrep"] < 10)
    if short:
        print("!! %s: %d cell(s) with fewer than 10 replicates, isl %s -- "
              "their eff, efsd and ps are means over what finished"
              % (rundir, len(short), short))
    return G


EPSIE = {"e2": 1, "e3": 2, "e4": 3}


def load_run(path, ie=None):
    """packed .idx, a directory holding packed pairs, or a raw run/

    `ie` = 1/2/3 for eps = 1e-2/1e-3/1e-4 picks the layer when the directory
    holds more than one.  The narrowed u=1e-2, eps=1e-4 campaign leaves
    exactly that situation behind: run/ carries the 455 eps=1e-4 cells plus
    the single eps=1e-2 cell of probe task 37, so packing it writes both
    dw1_e4.* and a one-cell dw1_e2.*.
    """
    if path.endswith(".idx"):
        return load_packed(path)
    idx = sorted(glob.glob(os.path.join(path, "dw1_*.idx")))
    if idx and ie is not None:
        want = [i for i in idx
                if EPSIE.get(os.path.basename(i)[4:-4]) == ie]
        if len(want) != 1:
            raise ValueError("%s holds %s; none of them is ie = %d"
                             % (path, [os.path.basename(i) for i in idx], ie))
        return load_packed(want[0])
    if len(idx) == 1:
        return load_packed(idx[0])
    if len(idx) > 1:
        raise ValueError("%s holds %d packed error rates (%s); pass ie=1/2/3 "
                         "for eps=1e-2/1e-3/1e-4, or name the .idx"
                         % (path, len(idx), ", ".join(os.path.basename(i)
                                                      for i in idx)))
    return load_raw(path, ie=ie)


if __name__ == "__main__":
    import sys
    ie = int(sys.argv[2]) if len(sys.argv) > 2 else None
    G = load_run(sys.argv[1], ie=ie)
    print("%d cells" % len(G))
    for gid in sorted(G)[:40]:
        g = G[gid]
        print("  gid %5d %-22s u %8.3f v %8.3f  Emax %8.4f  ef %8.5f "
              "+-%7.5f  lead %5d %8.5f  Neff %9.2f  nrep %2d"
              % (gid, g["tag"], g["u"], g["v"], g["emax"], g["eff"],
                 g["efsd"], g["top"][0], g["toppi"][0], g["neff"],
                 g["nrep"]))
