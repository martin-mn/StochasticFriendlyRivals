/* kappa.c -- the KERNEL MEASURE of the eight atoms and 81 strata at one task's game: the neutral
 * reference for enrichment (share / share-of-strategy-space).
 *
 *   kappa.x <isl> [M]        games.dat in the cwd; M kernel draws (default 10^7)
 *
 * Task decode is sf1.f's:  ie = mod(isl-1,3)+1, eps = (1e-2,1e-3,1e-4)(ie);  ig = (isl-1)/3+1 is the
 * data row of games.dat;  u = cS, v = cT - 1 (both exact: every u, v is a multiple of 1/1024).
 *
 * The kernel is the one the simulation uses for mutants and the initial population: m1m2.F's umut4
 * with the EXACT alpha construction (-DALPHA path),
 *     pmin = eps, pmax = 1 - pmin, uw = (1-2 pmin)/(1-2 alpha), ua = pmin - alpha uw,
 *     ran = ua + uw U  clipped to [pmin, pmax],      one U per coordinate,
 * so each coordinate is exactly pmin w.p. alpha, exactly pmax w.p. alpha, Uniform(pmin,pmax) otherwise.
 * U has sf1's buffer convention (mtb.f mtfill = genrand_real3: (r + 0.5)/2^32).
 *
 * The RNG is this program's OWN MT19937, init_genrand(41000000 + isl), independent of the dynamics
 * (whose seeds are 40000000 + ... and which never run in this process).
 *
 * Output: ONE line
 *   isl gid ie eps u v M  f_atom[0..7]  f_stratum[0..80]  f_neartie1  f_exacttie2
 * (98 fields; the last two are the fractions with m1cls neartie == 1 and == 2, see m1cls.h)
 * with atom index = 4E+2N+R (binary ENR) and stratum = sum_k code_k 3^k (code 0: p=0, 1: p=1,
 * 2: interior).  Fractions are counts/M printed with 12 significant digits (counts exactly recoverable).
 *
 * Build:  cc -O2 -o kappa.x kappa.c m1cls.c -lm      (no -ffast-math: m1cls relies on exact TwoSum)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "m1cls.h"

#ifndef ALPHA
#define ALPHA (1.0/3.0)
#endif

/* ------------------------------------------------------------------ MT19937 (reference algorithm) */
#define MTN 624
#define MTM 397
static uint32_t mt[MTN];
static int mti = MTN + 1;
static void init_genrand(uint32_t s)
{
    mt[0] = s;
    for (mti = 1; mti < MTN; mti++)
        mt[mti] = 1812433253u * (mt[mti-1] ^ (mt[mti-1] >> 30)) + (uint32_t)mti;
}
static uint32_t genrand_int32(void)
{
    static const uint32_t mag01[2] = {0u, 0x9908b0dfu};
    uint32_t y;
    if (mti >= MTN) {
        int kk;
        for (kk = 0; kk < MTN - MTM; kk++) {
            y = (mt[kk] & 0x80000000u) | (mt[kk+1] & 0x7fffffffu);
            mt[kk] = mt[kk+MTM] ^ (y >> 1) ^ mag01[y & 1u];
        }
        for (; kk < MTN - 1; kk++) {
            y = (mt[kk] & 0x80000000u) | (mt[kk+1] & 0x7fffffffu);
            mt[kk] = mt[kk+(MTM-MTN)] ^ (y >> 1) ^ mag01[y & 1u];
        }
        y = (mt[MTN-1] & 0x80000000u) | (mt[0] & 0x7fffffffu);
        mt[MTN-1] = mt[MTM-1] ^ (y >> 1) ^ mag01[y & 1u];
        mti = 0;
    }
    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680u;
    y ^= (y << 15) & 0xefc60000u;
    y ^= (y >> 18);
    return y;
}
static inline double genrand_real3(void) { return ((double)genrand_int32() + 0.5) * (1.0 / 4294967296.0); }

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <isl> [M]\n", argv[0]); return 1; }
    long isl = atol(argv[1]);
    long long M = (argc > 2) ? (long long)atof(argv[2]) : 10000000LL;
    if (isl < 1 || M < 1) { fprintf(stderr, "bad isl or M\n"); return 1; }
    static const double epsv[3] = {1e-2, 1e-3, 1e-4};
    int ie = (int)((isl - 1) % 3) + 1;
    long ig = (isl - 1) / 3 + 1;
    double eps = epsv[ie - 1];

    FILE *fg = fopen("games.dat", "r");
    if (!fg) { fprintf(stderr, "games.dat not found\n"); return 2; }
    char line[1024]; long row = 0; int gid = -1; double cR = 0, cS = 0, cT = 0;
    while (fgets(line, sizeof line, fg)) {
        char *q = line; while (*q == ' ' || *q == '\t') q++;
        if (*q == '#' || *q == '\n' || *q == 0) continue;
        if (++row == ig) {
            int ic; double uu, vv;
            if (sscanf(q, "%d %d %lf %lf %lf %lf %lf", &gid, &ic, &cR, &cS, &cT, &uu, &vv) != 7) {
                fprintf(stderr, "games.dat row %ld unreadable\n", row); return 2;
            }
            break;
        }
    }
    fclose(fg);
    if (gid < 0) { fprintf(stderr, "games.dat has no row %ld\n", ig); return 2; }
    if (cR != 1.0) { fprintf(stderr, "row %ld: cR != 1\n", ig); return 2; }
    double u = cS, v = cT - 1.0;

    const double pmin = eps, pmax = 1.0 - pmin;
    const double uw = (1.0 - 2.0 * pmin) / (1.0 - 2.0 * ALPHA);
    const double ua = pmin - ALPHA * uw;

    init_genrand((uint32_t)(41000000L + isl));
    long long na[8] = {0}, ns[81] = {0}, nnear = 0, ntie = 0;
    for (long long it = 0; it < M; it++) {
        double s[4];
        for (int k = 0; k < 4; k++) {
            double ran = ua + genrand_real3() * uw;
            s[k] = (ran < pmin) ? pmin : (ran > pmax) ? pmax : ran;
        }
        int a, st, nt;
        m1cls(s, &eps, &u, &v, &a, &st, &nt);
        na[a]++; ns[st]++; nnear += (nt == 1); ntie += (nt == 2);
    }
    long long ta = 0, ts = 0;
    for (int i = 0; i < 8; i++) ta += na[i];
    for (int i = 0; i < 81; i++) ts += ns[i];
    if (ta != M || ts != M) { fprintf(stderr, "count invariant failed\n"); return 3; }

    printf("%ld %d %d %.17g %.17g %.17g %lld", isl, gid, ie, eps, u, v, M);
    for (int i = 0; i < 8; i++)  printf(" %.12g", (double)na[i] / (double)M);
    for (int i = 0; i < 81; i++) printf(" %.12g", (double)ns[i] / (double)M);
    printf(" %.12g %.12g\n", (double)nnear / (double)M, (double)ntie / (double)M);
    return 0;
}
