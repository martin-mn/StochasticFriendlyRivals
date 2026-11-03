// Exact eps->0 limit of the stationary distribution for a CONTINUOUS memory-one resident
// p = (p1,p2,p3,p4) in [0,1]^4 against (a) the 16 binary memory-one co-players and (b) all
// 65536 binary memory-two co-players, by GTH state reduction in leading-order
// (coefficient, exponent) arithmetic.  Every operation is on positive quantities, so the
// leading term of each entry is exact and there is no threshold and no extrapolation.
//
// This is ../../PartnersRivals/check/pairs.c generalised from a binary resident to a
// continuous one.  The only change is the alphabet of leading terms:
//
//     p_i = 0  ->  C branch (coefficient 1, exponent 1),  D branch (1, 0)
//     p_i = 1  ->  C branch (1, 0),                       D branch (1, 1)
//     0<p_i<1  ->  C branch (p_i, 0),                     D branch (1-p_i, 0)
//
// which is exactly the eps-clamp p^eps = eps + (1-2 eps) p read to leading order.  Products of
// two such terms multiply coefficients and add exponents, so the pair chain's entries are
// leading terms again and GTH goes through verbatim.
//
// Convention (../../PartnersRivalsM1): outcomes CC=0 CD=1 DC=2 DD=3 with the FOCAL player's
// action first; a memory-one strategy answers the last outcome; a memory-two strategy answers
// state j = 4*recent + before and bit j of the 16-bit genome is 1 for C.
//
// Build:   cc -O3 -march=native -o m1f_pairs m1f_pairs.c -lm
// Usage:   m1f_pairs m1 <n> <task> <ntasks>        grid resolution 1/n, vs the 16 vertices
//          m1f_pairs m2 <n> <task> <ntasks>        grid resolution 1/n, vs all 65536 memory-two

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EINF 1000000000

// ---------------------------------------------------------------- leading-order GTH, size NS
// c[i][j] is the coefficient and e[i][j] the eps-exponent of the transition i->j; e = EINF
// means the entry is identically zero.  On return w[] holds the limiting distribution.
static void gth(int NS, double *c, int *e, double *w)
{
#define C(i,j) c[(i)*NS+(j)]
#define E(i,j) e[(i)*NS+(j)]
    for (int k = NS - 1; k >= 1; k--) {
        int es = EINF; double cs = 0;
        for (int j = 0; j < k; j++) {
            int ej = E(k,j);
            if (ej < es) { es = ej; cs = C(k,j); }
            else if (ej == es && es < EINF) cs += C(k,j);
        }
        for (int i = 0; i < k; i++)
            if (E(i,k) < EINF) { C(i,k) /= cs; E(i,k) -= es; }
        for (int i = 0; i < k; i++) {
            if (E(i,k) >= EINF) continue;
            int eik = E(i,k); double cik = C(i,k);
            for (int j = 0; j < k; j++) {
                if (E(k,j) >= EINF) continue;
                int en = eik + E(k,j); double cn = cik * C(k,j);
                if (en < E(i,j)) { E(i,j) = en; C(i,j) = cn; }
                else if (en == E(i,j)) C(i,j) += cn;
            }
        }
    }
    double pc[16]; int pe[16];                       // NS is 4 or 16, never more
    pc[0] = 1; pe[0] = 0;
    for (int k = 1; k < NS; k++) {
        int ek = EINF; double ck = 0;
        for (int i = 0; i < k; i++) {
            if (E(i,k) >= EINF || pe[i] >= EINF) continue;
            int en = pe[i] + E(i,k); double cn = pc[i] * C(i,k);
            if (en < ek) { ek = en; ck = cn; }
            else if (en == ek) ck += cn;
        }
        pc[k] = ck; pe[k] = ek;
    }
    int emin = EINF; double ct = 0;
    for (int k = 0; k < NS; k++) {
        if (pe[k] < emin) { emin = pe[k]; ct = pc[k]; }
        else if (pe[k] == emin) ct += pc[k];
    }
    for (int k = 0; k < NS; k++) w[k] = (pe[k] == emin) ? pc[k] / ct : 0.0;
#undef C
#undef E
}

// ------------------------------------------------------- leading term of one C-probability
// x is a probability in [0,1]; EXACTLY 0 and EXACTLY 1 are the degenerate cases that pick up
// a factor of eps.  cC/eC is the C branch, cD/eD the D branch.
static inline void branch(double x, double *cC, int *eC, double *cD, int *eD)
{
    if (x == 0.0)      { *cC = 1.0; *eC = 1; *cD = 1.0;     *eD = 0; }
    else if (x == 1.0) { *cC = 1.0; *eC = 0; *cD = 1.0;     *eD = 1; }
    else               { *cC = x;   *eC = 0; *cD = 1.0 - x; *eD = 0; }
}

static inline int swapo(int o) { return ((o & 1) << 1) | (o >> 1); }   // CD <-> DC

// ------------------------------------- four-state chain: continuous p vs binary memory-one q
// W[0..3] = limiting frequencies of CC, CD, DC, DD from the focal player's side.
static void freqs_cont(const double p[4], const double q[4], double W[4])
{
    double c[16]; int e[16];
    for (int i = 0; i < 16; i++) { c[i] = 0; e[i] = EINF; }
    for (int s = 0; s < 4; s++) {
        double pC, pD; int epC, epD;
        branch(p[s], &pC, &epC, &pD, &epD);
        double qC, qD; int eqC, eqD;
        branch(q[swapo(s)], &qC, &eqC, &qD, &eqD);
        // next state = 2*(focal plays D) + (co plays D)
        c[s*4 + 0] = pC * qC; e[s*4 + 0] = epC + eqC;
        c[s*4 + 1] = pC * qD; e[s*4 + 1] = epC + eqD;
        c[s*4 + 2] = pD * qC; e[s*4 + 2] = epD + eqC;
        c[s*4 + 3] = pD * qD; e[s*4 + 3] = epD + eqD;
    }
    double w[4];
    gth(4, c, e, w);
    for (int i = 0; i < 4; i++) W[i] = w[i];
}

static void freqs_m1(const double p[4], unsigned q, double W[4])
{
    double qq[4];
    for (int i = 0; i < 4; i++) qq[i] = (double)((q >> i) & 1u);
    freqs_cont(p, qq, W);
}

// ------------------------------- sixteen-state chain: continuous p vs binary memory-two tau
// State j = 4*recent + before, outcomes from the focal player's side.  The memory-one resident
// answers `recent` alone; the memory-two co-player answers bit jt of tau.
static void freqs_m2(const double p[4], unsigned tau, double W[4])
{
    double c[256]; int e[256];
    for (int i = 0; i < 256; i++) { c[i] = 0; e[i] = EINF; }
    for (int j = 0; j < 16; j++) {
        int o1 = j >> 2, o0 = j & 3;
        int jt = 4 * swapo(o1) + swapo(o0);
        double pC, pD; int epC, epD;
        branch(p[o1], &pC, &epC, &pD, &epD);
        double tC, tD; int etC, etD;
        branch((double)((tau >> jt) & 1u), &tC, &etC, &tD, &etD);
        // the four successors 4*onew + o1 are distinct, so these are assignments, not sums
        int n;
        n = 4 * 0 + o1; c[j*16 + n] = pC * tC; e[j*16 + n] = epC + etC;
        n = 4 * 1 + o1; c[j*16 + n] = pC * tD; e[j*16 + n] = epC + etD;
        n = 4 * 2 + o1; c[j*16 + n] = pD * tC; e[j*16 + n] = epD + etC;
        n = 4 * 3 + o1; c[j*16 + n] = pD * tD; e[j*16 + n] = epD + etD;
    }
    double w[16];
    gth(16, c, e, w);
    W[0] = W[1] = W[2] = W[3] = 0;
    for (int j = 0; j < 16; j++) W[j >> 2] += w[j];
}

// ------------------------------------------------------------------ the closed form under test
// Riv+ (T > S) : p4 = 0, p2 + p3 <= 1, not (p2 = 1 and p3 = 0)
// Riv- (T < S) : p1 = 1, p2 + p3 >= 1, not (p2 = 1 and p3 = 0)
static inline int pred_plus(const double p[4])
{
    return p[3] == 0.0 && p[1] + p[2] <= 1.0 && !(p[1] == 1.0 && p[2] == 0.0);
}
static inline int pred_minus(const double p[4])
{
    return p[0] == 1.0 && p[1] + p[2] >= 1.0 && !(p[1] == 1.0 && p[2] == 0.0);
}

// GTH returns exact leading coefficients, so a limit is either 0 or bounded away from it.
// 1e-12 is far below any coefficient the arithmetic can produce and far above rounding.
#define ZTOL 1e-12

// splitmix64: a deterministic PRNG, so a task's random sample is reproducible from its index
static inline unsigned long long sm64(unsigned long long *s)
{
    unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline double sm64u(unsigned long long *s) { return (double)(sm64(s) >> 11) / 9007199254740992.0; }


// ---------------------------------------------------- item 2: efficiency and friendly rivals
// Predicate-FREE search: for each p we compute the self-play limit and the two rival margins and
// count friendly rivals wedge by wedge.  Nothing here depends on any conjectured closed form, so
// a hit in W, N or E would refute the impossibility outright.
//   W = {u+v<1, u-v<1}: efficient = all-C  and rival for T>S
//   S = {u+v<1, u-v>1}: efficient = all-C  and rival for T<S
//   N = {u+v>1, u-v<1}: efficient = alternation and rival for T>S
//   E = {u+v>1, u-v>1}: efficient = alternation and rival for T<S
static void classify(const double p[4], int *effCC, int *effALT, int *effLINE,
                     int *rivP, int *rivM)
{
    double W[4];
    freqs_cont(p, p, W);                       // self-play
    *effCC   = (W[0] > 1.0 - ZTOL);
    *effALT  = (fabs(W[1] - 0.5) < ZTOL) && (fabs(W[2] - 0.5) < ZTOL);
    *effLINE = (W[3] < ZTOL);
    double mp = 1e300, mm = 1e300;
    for (unsigned q = 0; q < 16u; q++) {
        double V[4];
        freqs_m1(p, q, V);
        double d = V[2] - V[1];
        if (d < mp) mp = d;
        if (-d < mm) mm = -d;
    }
    *rivP = (mp >= -ZTOL);
    *rivM = (mm >= -ZTOL);
}

static int run_fr(long n, long task, long ntasks, int random_mode)
{
    long side = n + 1, total = random_mode ? n : side * side * side * side;
    long nW = 0, nS = 0, nN = 0, nE = 0, nEffCC = 0, nEffALT = 0, nEffLINE = 0, tested = 0;
    long listed = 0;
    unsigned long long rng = 0x9E3779B97F4A7C15ULL + 0x100000001B3ULL * (unsigned long long)task;
    for (long idx = task; idx < total; idx += ntasks) {
        double p[4];
        if (random_mode) {
            for (int i = 0; i < 4; i++) {          // draw a stratum, then a value inside it
                unsigned long long r = sm64(&rng) % 5ULL;
                p[i] = (r == 0) ? 0.0 : (r == 1) ? 1.0 : sm64u(&rng);
            }
        } else {
            long t = idx;
            long i1 = t % side; t /= side;
            long i2 = t % side; t /= side;
            long i3 = t % side; t /= side;
            p[0] = (double)i1 / n; p[1] = (double)i2 / n;
            p[2] = (double)i3 / n; p[3] = (double)t / n;
        }
        int eCC, eALT, eLINE, rP, rM;
        classify(p, &eCC, &eALT, &eLINE, &rP, &rM);
        nEffCC += eCC; nEffALT += eALT; nEffLINE += eLINE;
        int inW = eCC && rP, inS = eCC && rM, inN = eALT && rP, inE = eALT && rM;
        nW += inW; nS += inS; nN += inN; nE += inE;
        if ((inW || inN || inE) && listed < 200) {   // any of these refutes the impossibility
            printf("IMPOSSIBILITY-HIT p=(%.17g,%.17g,%.17g,%.17g) W=%d N=%d E=%d\n",
                   p[0], p[1], p[2], p[3], inW, inN, inE);
            listed++;
        }
        tested++;
    }
    printf("FR-SUMMARY task=%ld ntasks=%ld %s n=%ld tested=%ld "
           "effCC=%ld effALT=%ld effLINE=%ld  FR_W=%ld FR_S=%ld FR_N=%ld FR_E=%ld\n",
           task, ntasks, random_mode ? "frand" : "fr", n, tested,
           nEffCC, nEffALT, nEffLINE, nW, nS, nN, nE);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 5) {
        fprintf(stderr,
            "usage: %s m1|m2 <n> <task> <ntasks>       grid of spacing 1/n\n"
            "       %s rand|rand2 <N> <task> <ntasks>  N random points per task, all 81 strata\n"
            "       %s near <N> <task> <ntasks>        N points at controlled distance from the faces\n"
            "       %s fr|frand <n|N> <task> <ntasks>   item 2: friendly rivals per wedge, predicate-free\n",
            argv[0], argv[0], argv[0], argv[0]);
        return 1;
    }
    const char *mode = argv[1];
    long n = atol(argv[2]), task = atol(argv[3]), ntasks = atol(argv[4]);
    if (!strcmp(mode, "fr"))    return run_fr(n, task, ntasks, 0);
    if (!strcmp(mode, "frand")) return run_fr(n, task, ntasks, 1);
    long side = n + 1, total = side * side * side * side;
    int use_m2 = !strcmp(mode, "m2") || !strcmp(mode, "rand2");
    int is_rand = !strcmp(mode, "rand") || !strcmp(mode, "rand2");
    int is_near = !strcmp(mode, "near");
    long ncop = use_m2 ? 65536 : 16;
    if (is_rand || is_near) total = n;                       // n is then the sample size

    long tested = 0, bad_plus = 0, bad_minus = 0, listed = 0;
    long nriv_plus = 0, nriv_minus = 0;
    double smallest_nonzero = 1e300, largest_violation = 0.0;
    unsigned long long rng = 0x243F6A8885A308D3ULL + 0x100000001B3ULL * (unsigned long long)task;

    for (long idx = task; idx < total; idx += ntasks) {
        double p[4];
        if (is_rand) {
            // draw a stratum first, so all 81 of them are sampled and not only the interior one
            for (int i = 0; i < 4; i++) {
                unsigned long long r = sm64(&rng) % 5ULL;
                p[i] = (r == 0) ? 0.0 : (r == 1) ? 1.0 : sm64u(&rng);
            }
        } else if (is_near) {
            // a random point pushed to within delta of one of the structural sets, delta
            // sweeping 1e-1 down to 1e-12, to document how the margin degrades at the boundary
            double delta = pow(10.0, -(double)(1 + (sm64(&rng) % 12ULL)));
            for (int i = 0; i < 4; i++) p[i] = sm64u(&rng);
            switch (sm64(&rng) % 3ULL) {
                case 0: p[3] = delta; break;                                  // near the face p4 = 0
                case 1: p[0] = 1.0 - delta; break;                            // near the face p1 = 1
                default: p[2] = 1.0 - p[1] + delta;                           // near p2 + p3 = 1
                         if (p[2] < 0) p[2] = 0; if (p[2] > 1) p[2] = 1;
            }
        } else {
            long t = idx;
            long i1 = t % side; t /= side;
            long i2 = t % side; t /= side;
            long i3 = t % side; t /= side;
            long i4 = t;
            p[0] = (double)i1 / n; p[1] = (double)i2 / n;
            p[2] = (double)i3 / n; p[3] = (double)i4 / n;
        }

        double mplus = 1e300, mminus = 1e300;
        for (long q = 0; q < ncop; q++) {
            double W[4];
            if (use_m2) freqs_m2(p, (unsigned)q, W); else freqs_m1(p, (unsigned)q, W);
            double d = W[2] - W[1];
            double ad = fabs(d);
            if (ad > ZTOL && ad < smallest_nonzero) smallest_nonzero = ad;
            if (d < mplus) mplus = d;
            if (-d < mminus) mminus = -d;
        }
        int comp_plus  = (mplus  >= -ZTOL);
        int comp_minus = (mminus >= -ZTOL);
        int pp = pred_plus(p), pm = pred_minus(p);
        nriv_plus += comp_plus; nriv_minus += comp_minus;
        if (comp_plus != pp || comp_minus != pm) {
            if (comp_plus != pp)  bad_plus++;
            if (comp_minus != pm) bad_minus++;
            double viol = 0;
            if (comp_plus != pp)  viol = fabs(mplus);
            if (comp_minus != pm && fabs(mminus) > viol) viol = fabs(mminus);
            if (viol > largest_violation) largest_violation = viol;
            if (listed < 200) {
                printf("DISAGREE p=(%.6f,%.6f,%.6f,%.6f) comp=(%d,%d) pred=(%d,%d) "
                       "margins=(%.6g,%.6g)\n",
                       p[0], p[1], p[2], p[3], comp_plus, comp_minus, pp, pm, mplus, mminus);
                listed++;
            }
        }
        tested++;
    }
    printf("SUMMARY task=%ld ntasks=%ld mode=%s n=%ld tested=%ld "
           "bad_plus=%ld bad_minus=%ld rivals_plus=%ld rivals_minus=%ld "
           "smallest_nonzero=%.6g largest_violation=%.6g\n",
           task, ntasks, mode, n, tested, bad_plus, bad_minus, nriv_plus, nriv_minus,
           smallest_nonzero == 1e300 ? 0.0 : smallest_nonzero, largest_violation);
    return 0;
}
