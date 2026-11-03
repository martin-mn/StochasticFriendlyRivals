#!/usr/bin/env python3
"""Process the Cannon stratum census (array 45751662) into the friendly-rival tables.

Input: the concatenated out/st_*.txt of the census, one line per stratum with any confirmed
rival:   g  class  nP  nM  K
where class is C (self-play W_CC = 1), A (perfect alternation), Z (W_DD = 0 only) or '-' (not
efficient anywhere), and nP/nM count the draws confirmed as rivals for T>S / T<S out of K.

    python3 census.py st_all.txt

Reports, per locus of the plane: the number of support strata carrying friendly rivals, the
dimension histogram, the binary vertex count, and the maximal strata.  Also applies the two exact
consistency checks:
  (1) the C/D mirror mbar(s) = 15 - s is a bijection fixing the alternating class and swapping the
      rival signs, so in a COMPLETE census wedge N and wedge E must have equal counts and their
      stratum sets must be exact mirror images;
  (2) every binary friendly rival must appear as a zero-dimensional stratum, and there must be no
      others -- 8 / 1519 / 80 / 80 for W / S / N / E.
"""
import sys, collections
import numpy as np

COMP = {0: 2, 1: 1, 2: 0}
SYM = {0: '0', 1: '*', 2: '1'}

def dig(g):
    d = []
    for _ in range(16):
        d.append(g % 3); g //= 3
    return d

def undig(d):
    g = 0
    for j in range(15, -1, -1): g = g*3 + d[j]
    return g

def mirror(g):
    d = dig(g)
    return undig([COMP[d[15-s]] for s in range(16)])

def nfree(g): return sum(1 for x in dig(g) if x == 1)

def pat(g): return "".join(SYM[x] for x in dig(g))

def cell_mask_value(g):
    """A stratum's closed cell contains binary genome b iff (b & M) == V, where M marks the
    coordinates the stratum pins and V their values.  Free coordinates impose nothing."""
    d = dig(g); M = 0; V = 0
    for j in range(16):
        if d[j] != 1:
            M |= 1 << j
            if d[j] == 2: V |= 1 << j
    return M, V

def visible_flags(strata, ref):
    """For each stratum, does its closed cell contain any genome of `ref`?
    Grouped by the pin-mask M so the cost is O(#distinct M x |ref| + #strata) rather than
    O(#strata x 2^free)."""
    ref = np.asarray(sorted(ref), dtype=np.int64)
    mv = [cell_mask_value(g) for g in strata]
    bym = {}
    for M, V in mv: bym.setdefault(M, set()).add(V)
    ok = {}
    for M in bym:
        ok[M] = set(np.unique(ref & M).tolist())
    return [V in ok[M] for (M, V) in mv]

def main(path):
    C_P, C_M, A_P, A_M, Z_P, Z_M = (set() for _ in range(6))
    nline = 0
    for line in open(path):
        a = line.split()
        if len(a) < 4: continue
        g = int(a[0]); cls = a[1]; nP = int(a[2]); nM = int(a[3]); nline += 1
        if cls == 'C':
            if nP: C_P.add(g)
            if nM: C_M.add(g)
        elif cls == 'A':
            if nP: A_P.add(g)
            if nM: A_M.add(g)
        elif cls == 'Z':
            if nP: Z_P.add(g)
            if nM: Z_M.add(g)
    print("census lines read: %d" % nline)

    loci = [("wedge W  (u+v<1, T>S)", C_P),
            ("wedge S  (u+v<1, T<S)", C_M),
            ("wedge N  (u+v>1, T>S)", A_P),
            ("wedge E  (u+v>1, T<S)", A_M),
            ("ray u+v=1, u<1 (T>S)", C_P | A_P | Z_P),
            ("ray u+v=1, u>1 (T<S)", C_M | A_M | Z_M)]

    M = np.load('/Users/martin/Documents/CL1/PartnersRivals/Figures/m2_masks.npz')
    effCC, effALT, rivP, rivM = (M[k] for k in ("effCC", "effALT", "rivP", "rivM"))
    binloci = {"wedge W  (u+v<1, T>S)": set(np.where(effCC & rivP)[0]),
               "wedge S  (u+v<1, T<S)": set(np.where(effCC & rivM)[0]),
               "wedge N  (u+v>1, T>S)": set(np.where(effALT & rivP)[0]),
               "wedge E  (u+v>1, T<S)": set(np.where(effALT & rivM)[0])}

    print("\n%-24s %10s  %6s  %s" % ("locus", "strata", "maxdim", "dimension histogram"))
    for name, S in loci:
        if not S:
            print("%-24s %10d" % (name, 0)); continue
        h = collections.Counter(nfree(g) for g in S)
        print("%-24s %10d  %6d  %s" % (name, len(S), max(h),
              " ".join("%d:%d" % (k, h[k]) for k in sorted(h))))

    print("\n--- exact check 1: the C/D mirror must map wedge N onto wedge E ---")
    mN = {mirror(g) for g in A_P}
    print("  |N| = %d, |E| = %d, equal: %s" % (len(A_P), len(A_M), len(A_P) == len(A_M)))
    print("  mirror(N) == E as SETS: %s   (|mirror(N) \\ E| = %d, |E \\ mirror(N)| = %d)"
          % (mN == A_M, len(mN - A_M), len(A_M - mN)))
    print("  union N u mirror(E), the tightened lower bound: %d" % len(A_P | {mirror(g) for g in A_M}))

    print("\n--- exact check 2: binary vertices are exactly the 0-dimensional strata ---")
    for name, S in loci[:4]:
        zero = {g for g in S if nfree(g) == 0}
        # a 0-dimensional stratum IS a binary strategy: its genome is the bits with digit 2
        binfound = {sum(1 << j for j, x in enumerate(dig(g)) if x == 2) for g in zero}
        ref = binloci[name]
        print("  %-24s 0-dim strata %5d   reference binary count %5d   set-equal: %s"
              % (name, len(zero), len(ref), binfound == ref))

    print("\n--- binary-invisible components (closed cell contains no binary friendly rival) ---")
    for name, S in loci[:4]:
        ref = binloci[name]
        Sl = sorted(S)
        vis = visible_flags(Sl, ref)
        inv = [g for g, v in zip(Sl, vis) if not v]
        if inv:
            h = collections.Counter(nfree(g) for g in inv)
            print("  %-24s invisible %6d of %6d   maxdim %d" % (name, len(inv), len(S), max(h)))
        else:
            print("  %-24s invisible      0 of %6d" % (name, len(S)))

    print("\n--- maximal strata per wedge ---")
    for name, S in loci[:4]:
        if not S: continue
        mx = max(nfree(g) for g in S)
        ex = [g for g in S if nfree(g) == mx][:4]
        print("  %-24s dim %2d : %s" % (name, mx, "  ".join(pat(g) for g in ex)))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "st_all.txt")
