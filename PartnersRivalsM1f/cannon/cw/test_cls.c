/* test_cls.c -- C test driver and CLI for m1cls (used by test_cls.py).
 *
 *   test_cls.x corners [games.dat]    V1 input: for each of the 512 games and each of the 16 corners
 *                                     (sf1 bit convention, s_k = 1-eps if bit k else eps, eps = 1e-4)
 *                                     one line "ipt k atom stratum neartie"; asserts early == full.
 *   test_cls.x points                 stdin: lines "u v eps s1 s2 s3 s4" (%la hex or decimal);
 *                                     stdout per line: "atom stratum neartie E N R maxgap_full
 *                                     atom_e stratum_e neartie_e ndev_e p1 p2 p3 p4" (%a for doubles)
 *   test_cls.x bench u v M [eps]      mean ns per m1cls call over M exact-alpha kernel draws, with
 *                                     and without early exit, plus the atom shares of the draws
 *   test_cls.x selftest               internal consistency checks (exit status 0 = pass)
 *
 * Build:  cc -O2 -o test_cls.x test_cls.c m1cls.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include "m1cls.h"

static int read_games(const char *fn, double *u, double *v, int max)
{
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", fn); exit(2); }
    char line[1024]; int n = 0;
    while (fgets(line, sizeof line, f)) {
        char *q = line; while (*q == ' ' || *q == '\t') q++;
        if (*q == '#' || *q == '\n' || *q == 0) continue;
        int gid, ic; double cR, cS, cT, uu, vv;
        if (sscanf(q, "%d %d %lf %lf %lf %lf %lf", &gid, &ic, &cR, &cS, &cT, &uu, &vv) != 7) exit(2);
        if (n >= max) exit(2);
        u[n] = cS; v[n] = cT - 1.0; n++;
    }
    fclose(f);
    return n;
}

static void cls_both(const double s[4], double eps, double u, double v, int *out /* 12 ints */,
                     double *mg, double p[4])
{
    int a, st, nt, E, N, R, nd, a2, st2, nt2, E2, N2, R2, nd2; double mg2, p2[4];
    m1cls_detail(s, eps, u, v, 0, &a, &st, &nt, mg, p, &E, &N, &R, &nd);
    m1cls_detail(s, eps, u, v, 1, &a2, &st2, &nt2, &mg2, p2, &E2, &N2, &R2, &nd2);
    out[0] = a; out[1] = st; out[2] = nt; out[3] = E; out[4] = N; out[5] = R;
    out[6] = a2; out[7] = st2; out[8] = nt2; out[9] = nd2;
    /* the production entry point must agree with the detail call */
    int a3, st3, nt3;
    m1cls(s, &eps, &u, &v, &a3, &st3, &nt3);
    out[10] = (a3 == a2 && st3 == st2 && nt3 == nt2);
    out[11] = (a == a2 && st == st2 && nt == nt2);
}

static int do_corners(const char *fn)
{
    static double u[2048], v[2048];
    int ng = read_games(fn, u, v, 2048), bad = 0;
    const double eps = 1e-4;
    for (int g = 0; g < ng; g++)
        for (int k = 0; k < 16; k++) {
            double s[4], mg, p[4]; int o[12];
            for (int j = 0; j < 4; j++) s[j] = ((k >> j) & 1) ? 1.0 - eps : eps;
            cls_both(s, eps, u[g], v[g], o, &mg, p);
            if (!o[10] || !o[11]) bad++;
            printf("%d %d %d %d %d\n", g, k, o[0], o[1], o[2]);
        }
    if (bad) { fprintf(stderr, "early-exit mismatches at corners: %d\n", bad); return 1; }
    return 0;
}

static int do_points(void)
{
    char line[4096]; int bad = 0;
    while (fgets(line, sizeof line, stdin)) {
        double u, v, eps, s[4];
        char *q = line, *e;
        double *dst[7] = {&u, &v, &eps, &s[0], &s[1], &s[2], &s[3]};
        int ok = 1;
        for (int i = 0; i < 7; i++) { *dst[i] = strtod(q, &e); if (e == q) { ok = 0; break; } q = e; }
        if (!ok) continue;
        double mg, p[4]; int o[12];
        cls_both(s, eps, u, v, o, &mg, p);
        if (!o[10] || !o[11]) bad++;
        printf("%d %d %d %d %d %d %a %d %d %d %d %a %a %a %a\n", o[0], o[1], o[2], o[3], o[4], o[5],
               mg, o[6], o[7], o[8], o[9], p[0], p[1], p[2], p[3]);
    }
    fflush(stdout);
    if (bad) { fprintf(stderr, "early-exit mismatches: %d\n", bad); return 1; }
    return 0;
}

/* splitmix64 for the benchmark's kernel draws (independent of everything else) */
static uint64_t sm = 0x243F6A8885A308D3ULL;
static inline double smu(void)
{
    uint64_t z = (sm += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return ((double)(uint32_t)(z >> 32) + 0.5) * (1.0 / 4294967296.0);
}

static double now(void)
{
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static int do_bench(double u, double v, long M, double eps)
{
    const double a = 1.0 / 3.0, pmin = eps, pmax = 1.0 - pmin;
    const double uw = (1.0 - 2.0 * pmin) / (1.0 - 2.0 * a), ua = pmin - a * uw;
    double *S = malloc(sizeof(double) * 4 * M);
    for (long i = 0; i < 4 * M; i++) {
        double r = ua + smu() * uw;
        S[i] = r < pmin ? pmin : r > pmax ? pmax : r;
    }
    long na[8] = {0}, nnear = 0, ntie = 0, ndsum = 0, mism = 0;
    int *A1 = malloc(sizeof(int) * M);
    /* warm-up */
    for (long i = 0; i < (M < 10000 ? M : 10000); i++) {
        int at, st, nt; m1cls(S + 4 * i, &eps, &u, &v, &at, &st, &nt);
    }
    double t0 = now();
    for (long i = 0; i < M; i++) {
        int at, st, nt; m1cls(S + 4 * i, &eps, &u, &v, &at, &st, &nt);
        A1[i] = at; na[at]++; nnear += (nt == 1); ntie += (nt == 2);
    }
    double t1 = now();
    for (long i = 0; i < M; i++) {
        int at, st, nt, E, N, R, nd; double mg, p[4];
        m1cls_detail(S + 4 * i, eps, u, v, 0, &at, &st, &nt, &mg, p, &E, &N, &R, &nd);
        if (at != A1[i]) mism++;
    }
    double t2 = now();
    for (long i = 0; i < M; i++) {
        int at, st, nt, E, N, R, nd; double mg, p[4];
        m1cls_detail(S + 4 * i, eps, u, v, 1, &at, &st, &nt, &mg, p, &E, &N, &R, &nd);
        ndsum += nd;
    }
    printf("bench u=%.10g v=%.10g M=%ld eps=%g  early: %.1f ns/call (mean deviators %.3f)  "
           "full: %.1f ns/call  early-vs-full mismatches %ld  neartie(1) %ld  exact-tie(2) %ld\n",
           u, v, M, eps, 1e9 * (t1 - t0) / M, (double)ndsum / M, 1e9 * (t2 - t1) / M, mism, nnear, ntie);
    static const char *nm[8] = {"000", "001", "010", "011", "100", "101", "110", "111"};
    printf("  kernel atom shares (ENR):");
    for (int i = 0; i < 8; i++) printf(" %s=%.5f", nm[i], (double)na[i] / M);
    printf("\n");
    free(S); free(A1);
    return mism != 0;
}

/* a few hand checks with known answers */
static int do_selftest(void)
{
    int fail = 0;
    const double eps = 1e-4, om = 1.0 - eps;
    /* (1) self-play of WSLS: W = (1,0,0,0); TFT vs TFT: (1/4,1/4,1/4,1/4)... convention A
     *     TFT self-play at eps->0 is the uniform distribution over the 4 states. */
    double tft[4] = {1, 0, 1, 0}, wsls[4] = {1, 0, 0, 1}, W[4];
    m1cls_freqs(wsls, wsls, W);
    if (!(W[0] > 1 - 1e-15)) { printf("FAIL WSLS self-play %g %g %g %g\n", W[0], W[1], W[2], W[3]); fail++; }
    m1cls_freqs(tft, tft, W);
    for (int k = 0; k < 4; k++) if (!(W[k] > 0.25 - 1e-15 && W[k] < 0.25 + 1e-15)) { printf("FAIL TFT self-play\n"); fail++; break; }
    /* (2) friendly.md: (1,0,1,1/2) has W = (2/5, 3/10, 3/10, 0) in self-play */
    double x[4] = {1, 0, 1, 0.5};
    m1cls_freqs(x, x, W);
    if (fabs(W[0] - 0.4) > 1e-15 || fabs(W[1] - 0.3) > 1e-15 || fabs(W[2] - 0.3) > 1e-15 || W[3] != 0)
        { printf("FAIL (1,0,1,1/2) self-play %g %g %g %g\n", W[0], W[1], W[2], W[3]); fail++; }
    /* (3) the trap sub-face p1=1,p4=0 counterexample of nash.md: W_CC=3/5, W_DD=2/5 */
    /* (4) ALLD at the PD (u=-0.5, v=0.5): atom 011 = 3 (Nash, rival, not efficient) */
    double s[4] = {eps, eps, eps, eps}, u = -0.5, v = 0.5; int a, st, nt;
    m1cls(s, &eps, &u, &v, &a, &st, &nt);
    if (a != 3 || st != 0) { printf("FAIL ALLD at PD atom %d stratum %d\n", a, st); fail++; }
    /* (5) WSLS at a weak PD u=-0.1, v=0.1: efficient, Nash (2R>T+P), not rival -> 110 = 6 */
    double sw[4] = {om, eps, eps, om}; u = -0.1; v = 0.1;
    m1cls(sw, &eps, &u, &v, &a, &st, &nt);
    if (a != 6 || st != 1 + 27) { printf("FAIL WSLS atom %d stratum %d\n", a, st); fail++; }
    /* (6) the interior preimage stays interior next to the faces */
    double sn[4] = {nextafter(eps, 1), nextafter(om, 0), 0.5, 0.25};
    double mg, p[4]; int E, N, R, nd;
    m1cls_detail(sn, eps, u, v, 0, &a, &st, &nt, &mg, p, &E, &N, &R, &nd);
    if (st != 80 || !(p[0] > 0) || !(p[1] < 1)) { printf("FAIL near-face preimage st=%d p=%a %a\n", st, p[0], p[1]); fail++; }
    printf("selftest %s (%d failures)\n", fail ? "FAILED" : "passed", fail);
    return fail != 0;
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: see header\n"); return 1; }
    if (!strcmp(argv[1], "corners")) return do_corners(argc > 2 ? argv[2] : "games.dat");
    if (!strcmp(argv[1], "points")) return do_points();
    if (!strcmp(argv[1], "bench") && argc >= 5)
        return do_bench(atof(argv[2]), atof(argv[3]), atol(argv[4]), argc > 5 ? atof(argv[5]) : 1e-4);
    if (!strcmp(argv[1], "selftest")) return do_selftest();
    fprintf(stderr, "unknown mode\n");
    return 1;
}
