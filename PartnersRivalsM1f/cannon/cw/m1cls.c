/* m1cls.c -- the exact eps->0 (convention A) atom (E,N,R) of a continuous memory-one strategy.
 * Interface and conventions: m1cls.h.
 *
 * WHERE EACH BIT COMES FROM
 *   E  the closed forms of ../../Text/friendly.md sec.1 = ../../check/m1lib.py in_eff_cc (u+v<1),
 *      in_eff_alt (u+v>1), in_eff_line (u+v=1).  They depend on p only through its STRATUM (which
 *      coordinates are exactly 0, exactly 1, interior), so E is decided from the codes alone.
 *   R  m1lib in_riv_plus (T>S, u-v<1), in_riv_minus (T<S, u-v>1), everything on T=S.  The one
 *      non-stratum clause, p2+p3 <= 1 (>= 1), is decided EXACTLY: for two interior coordinates
 *      p2+p3-1 = (s2+s3-1)/(1-2eps) in real arithmetic, and sign(s2+s3-1) is computed exactly by
 *      TwoSum.  So R has no tolerance and no near-tie.
 *   N  Nash iff max_sigma [pi_lim(sigma,p) - pi_lim(p,p)] <= 0 over the 16 deterministic memory-one
 *      deviators (the quantifier collapse, Text/partners.md sec.1, Text/nash.md), each limit taken by
 *      the leading-order (coefficient, exponent) GTH of ../m1f_pairs.c.  Decided with
 *      tau = 1e-12*scale; near-tie band 1e-9*scale; scale = max(1,|u|,|1+v|).
 *      Early exit: the deviators are taken in the order ALLD, ALLC, then 1..14, and the loop stops at
 *      the first gap > 1e-9*scale.  Such a gap is a violation (> tau) AND puts the max outside the
 *      near-tie band, so all three outputs are exactly those of the full 16-deviator test.
 *      At a corner the deviator equal to p (gap identically 0) is skipped; it cannot change N.
 *      neartie = 2 if the max over the (other) deviators is exactly 0.0 -- a structural tie such as
 *      pi(ALLC,p) = pi(p,p) = 1 on the slab p1 = 1, decided by the weak inequality; = 1 if
 *      0 < |max| <= 1e-9*scale, the rounding-sensitive case; 0 otherwise.
 *
 * THE PORT FROM ../m1f_pairs.c (read-only), verbatim in logic:
 *   gth()        <- m1f_pairs.c lines 36-80 (specialised to NS = 4; the loop structure, the
 *                   leading-order min/sum rules and the final normalisation are unchanged)
 *   branch()     <- lines 85-90, unchanged, and evaluated once per strategy (struct lt4) instead of
 *                   once per chain entry; for an interior coordinate the D-branch coefficient 1-p is
 *                   formed from s directly (see m1cls_detail) rather than as 1.0 - (rounded p)
 *   swapo()      <- line 92, unchanged
 *   freqs()      <- freqs_cont, lines 96-114, unchanged except that the co-player q is any vector of
 *                   preimage values (freqs_m1, lines 116-121, is the special case of a binary q; here
 *                   the binary vector is built once per deviator in DEV[][])
 *   EINF         <- line 31
 * No threshold enters the GTH: every operation is on positive leading terms.
 */
#include <math.h>
#include "m1cls.h"

#define EINF 1000000000

/* ---------------------------------------------- leading-order GTH, 4 states (m1f_pairs.c 36-80) */
static void gth4(double *c, int *e, double *w)
{
#define NS 4
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
    double pc[NS]; int pe[NS];
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
#undef NS
}

/* ------------------------------------------ leading term of one C-probability (m1f_pairs.c 85-90) */
static inline void branch(double x, double *cC, int *eC, double *cD, int *eD)
{
    if (x == 0.0)      { *cC = 1.0; *eC = 1; *cD = 1.0;     *eD = 0; }
    else if (x == 1.0) { *cC = 1.0; *eC = 0; *cD = 1.0;     *eD = 1; }
    else               { *cC = x;   *eC = 0; *cD = 1.0 - x; *eD = 0; }
}

static inline int swapo(int o) { return ((o & 1) << 1) | (o >> 1); }   /* CD <-> DC (line 92) */

/* the four branch() results of one strategy, computed once and reused against every co-player */
typedef struct { double cC[4], cD[4]; int eC[4], eD[4]; } lt4;

static void lt_from_p(const double x[4], lt4 *L)
{
    for (int k = 0; k < 4; k++) branch(x[k], &L->cC[k], &L->eC[k], &L->cD[k], &L->eD[k]);
}

/* ----------------------------------------- four-state pair chain, focal x vs co-player y (96-114) */
static void freqs_lt(const lt4 *x, const lt4 *y, double W[4])
{
    double c[16]; int e[16];
    for (int s = 0; s < 4; s++) {
        const double pC = x->cC[s], pD = x->cD[s];
        const int epC = x->eC[s], epD = x->eD[s];
        const int so = swapo(s);
        const double qC = y->cC[so], qD = y->cD[so];
        const int eqC = y->eC[so], eqD = y->eD[so];
        /* next state = 2*(focal plays D) + (co plays D) */
        c[s*4 + 0] = pC * qC; e[s*4 + 0] = epC + eqC;
        c[s*4 + 1] = pC * qD; e[s*4 + 1] = epC + eqD;
        c[s*4 + 2] = pD * qC; e[s*4 + 2] = epD + eqC;
        c[s*4 + 3] = pD * qD; e[s*4 + 3] = epD + eqD;
    }
    gth4(c, e, W);
}

void m1cls_freqs(const double x[4], const double y[4], double W[4])
{
    lt4 X, Y;
    lt_from_p(x, &X); lt_from_p(y, &Y);
    freqs_lt(&X, &Y, W);
}

/* ------------------------------------------------------------------- exact sign of a + b - 1 */
/* TwoSum: a + b = x + t exactly; then x - 1 = y + t2 exactly.  y != 0 dominates t + t2 (both are at
 * most an ulp of y), and y == 0 means x == 1 exactly, when the sign is that of t. */
static inline int sgn_sum_minus1(double a, double b)
{
    double x = a + b;
    double bb = x - a;
    double t = (a - (x - bb)) + (b - bb);
    double y = x - 1.0;
    if (y > 0) return 1;
    if (y < 0) return -1;
    return (t > 0) - (t < 0);
}

/* the 16 deterministic deviators in the order ALLD, ALLC, 1..14 (sf1's bit convention: bit k -> C
 * after state k) as preimage vectors */
static const int DORD[16] = {0, 15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
static lt4 DEVL[16];
static int dev_ready = 0;
static void dev_init(void)
{
    for (int i = 0; i < 16; i++) {
        double d[4];
        for (int k = 0; k < 4; k++) d[k] = (double)((DORD[i] >> k) & 1);
        lt_from_p(d, &DEVL[i]);
    }
    dev_ready = 1;
}

void m1cls_detail(const double s[4], double eps, double u, double v, int early,
                  int *atom, int *stratum, int *neartie,
                  double *maxgap, double p[4], int *Eo, int *No, int *Ro, int *ndev)
{
    if (!dev_ready) dev_init();
    const double om = 1.0 - eps, den = 1.0 - 2.0 * eps;
    const double omerr = (1.0 - om) - eps;          /* (1 - eps) - om exactly: both steps are exact */
    int c[4], st = 0, pw = 1;
    lt4 P;
    for (int k = 0; k < 4; k++) {
        if (s[k] == eps)     { c[k] = 0; p[k] = 0.0; }
        else if (s[k] == om) { c[k] = 1; p[k] = 1.0; }
        else {
            c[k] = 2;
            /* p = (s - eps)/(1 - 2eps) and ITS COMPLEMENT 1 - p = ((1 - eps) - s)/(1 - 2eps), each
             * from s with one rounding: s - eps is exact for s <= 2eps and om - s for s >= om/2
             * (Sterbenz), so both are relatively accurate next to their own face.  (1 - p computed
             * from a rounded p, as m1f_pairs.c's branch() does, loses all relative accuracy when p is
             * within ~1e-12 of 1.) */
            double x = (s[k] - eps) / den;
            double y = ((om - s[k]) + omerr) / den;
            /* keep an interior coordinate interior: rounding must not land it on a face */
            if (!(x > 0.0)) x = 2.2250738585072014e-308;   /* DBL_MIN: normal, FTZ-safe */
            if (!(x < 1.0)) x = 0.99999999999999989;       /* nextafter(1, 0) */
            if (!(y > 0.0)) y = 2.2250738585072014e-308;
            p[k] = x;
            P.cC[k] = x; P.eC[k] = 0; P.cD[k] = y; P.eD[k] = 0;   /* branch(), interior case */
        }
        if (c[k] < 2) branch(p[k], &P.cC[k], &P.eC[k], &P.cD[k], &P.eD[k]);
        st += c[k] * pw; pw *= 3;
    }
    *stratum = st;

    /* ---- E: from the stratum alone (m1lib in_eff_cc / in_eff_alt / in_eff_line) */
    const int alt = (c[1] == 0 && c[2] == 1) || (c[1] == 1 && c[2] == 0);
    const int effalt = alt && c[0] != 1 && c[3] != 0 && !(c[0] == 0 && c[3] == 1);
    const int sw = sgn_sum_minus1(u, v);                 /* sign(u + v - 1) */
    int E;
    if (sw < 0)
        E = (c[0] == 1) && !alt && (c[3] != 0 || c[1] == 1 || c[2] == 1);
    else if (sw > 0)
        E = effalt;
    else
        E = (c[0] == 1) ? (c[3] != 0 || ((c[1] == 1 || c[2] == 1) && c[1] != 0 && c[2] != 0))
                        : effalt;

    /* ---- R: m1lib in_riv_plus / in_riv_minus, the p2+p3 clause exact */
    const int ts = sgn_sum_minus1(u, -v);                /* sign(u - v - 1): <0 T>S, >0 T<S */
    int R;
    if (ts == 0) R = 1;
    else {
        int sg;                                          /* sign(p2 + p3 - 1), exact */
        if (c[1] == 2 && c[2] == 2) sg = sgn_sum_minus1(s[1], s[2]);
        else if (c[1] == 2 || c[2] == 2) {               /* one pinned, one interior */
            int pin = (c[1] == 2) ? c[2] : c[1];
            sg = pin ? 1 : -1;
        } else sg = c[1] + c[2] - 1;                     /* both pinned */
        const int trap = (c[1] == 1 && c[2] == 0);
        if (ts < 0) R = (c[3] == 0) && sg <= 0 && !trap;
        else        R = (c[0] == 1) && sg >= 0 && !trap;
    }

    /* ---- N: eps->0 Nash against the 16 deterministic memory-one deviators */
    double scale = 1.0;
    if (fabs(u) > scale) scale = fabs(u);
    if (fabs(1.0 + v) > scale) scale = fabs(1.0 + v);
    const double tau = 1e-12 * scale, band = 1e-9 * scale, T = 1.0 + v;
    double W[4];
    freqs_lt(&P, &P, W);
    const double own = W[0] + W[1] * u + W[2] * T;
    /* at a corner the deviator equal to p itself has gap exactly 0 <= tau: it cannot decide N, so it
     * is skipped (otherwise every Nash corner would read as a near-tie) */
    int self = -1;
    if (c[0] < 2 && c[1] < 2 && c[2] < 2 && c[3] < 2) self = c[0] | c[1] << 1 | c[2] << 2 | c[3] << 3;
    double mx = -1e300;
    int nd = 0;
    for (int i = 0; i < 16; i++) {
        if (DORD[i] == self) continue;
        double V[4];
        freqs_lt(&DEVL[i], &P, V);
        double g = (V[0] + V[1] * u + V[2] * T) - own;
        nd++;
        if (g > mx) mx = g;
        if (early && g > band) break;
    }
    const int N = (mx <= tau);
    *neartie = (mx == 0.0) ? 2 : (fabs(mx) <= band) ? 1 : 0;
    *maxgap = mx;
    *ndev = nd;
    *Eo = E; *No = N; *Ro = R;
    *atom = 4 * E + 2 * N + R;
}

void m1cls(const double s[4], const double *eps, const double *u, const double *v,
           int *atom, int *stratum, int *neartie)
{
    double mg, p[4]; int E, N, R, nd;
    m1cls_detail(s, *eps, *u, *v, 1, atom, stratum, neartie, &mg, p, &E, &N, &R, &nd);
}

void m1cls_(const double s[4], const double *eps, const double *u, const double *v,
            int *atom, int *stratum, int *neartie)
{
    m1cls(s, eps, u, v, atom, stratum, neartie);
}
