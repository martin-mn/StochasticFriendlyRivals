/* m1cls.h -- the exact eps->0 (convention A) atom of a continuous memory-one strategy.
 *
 * A strategy is handed over as the EFFECTIVE vector s in [eps,1-eps]^4 that the simulation uses
 * (P(C | CC,CD,DC,DD), own action first).  Its preimage in the closed cube is
 *     s_k == eps        -> p_k = 0            code 0
 *     s_k == 1.0 - eps  -> p_k = 1            code 1   (the double the kernel writes as pmax)
 *     otherwise         -> p_k = (s_k-eps)/(1-2eps), interior   code 2
 * and the classification is that of p in the limit eps->0, weak inequalities taken after the limit.
 *
 *     atom    = 4*E + 2*N + R          (binary ENR; '110' = 6 is a partner, '111' = 7 a friendly rival)
 *     stratum = sum_k code_k 3^k       (0..80; the 16 corners are the strata with no code 2)
 *     neartie = Nash margin m = max_sigma[pi(sigma,p)-pi(p,p)] (sigma != p at a corner):
 *               2 if m == 0.0 exactly (structural tie), 1 if 0 < |m| <= 1e-9*scale, 0 otherwise;
 *               scale = max(1,|u|,|1+v|).  Test "neartie != 0" for the spec's flag.
 *
 * Game (R,S,T,P) = (1, u, 1+v, 0).
 *
 * Fortran:  call m1cls(s, eps, u, v, iatom, istrat, inear)   with real*8 s(4), eps, u, v and default
 * integers.  Both link conventions are exported: the plain symbol `m1cls` (for an interface block with
 * bind(C, name='m1cls')) and the lowercase trailing-underscore symbol `m1cls_` that gfortran and ifx
 * (Linux) emit for an implicit external call from fixed-form code.  Arguments are by reference.
 */
#ifndef M1CLS_H
#define M1CLS_H

#ifdef __cplusplus
extern "C" {
#endif

/* the production entry point: Nash test with early exit (identical results to the full test) */
void m1cls (const double s[4], const double *eps, const double *u, const double *v,
            int *atom, int *stratum, int *neartie);
void m1cls_(const double s[4], const double *eps, const double *u, const double *v,
            int *atom, int *stratum, int *neartie);

/* the same, and also the pieces, for tests.  early = 0 evaluates all 16 deviators.
 * maxgap = max_sigma [pi_lim(sigma,p) - pi_lim(p,p)] over the deviators evaluated (with early exit
 * this is only a lower bound on the true max, and is > 1e-9*scale whenever the exit was taken);
 * p[4] = the preimage actually used; E, N, R the three bits; ndev = deviators evaluated. */
void m1cls_detail(const double s[4], double eps, double u, double v, int early,
                  int *atom, int *stratum, int *neartie,
                  double *maxgap, double p[4], int *E, int *N, int *R, int *ndev);

/* leading-order eps->0 limit of the pair chain, focal x (as p-preimage values in [0,1]) against
 * co-player y: W[0..3] = limiting frequencies of CC,CD,DC,DD from x's side. */
void m1cls_freqs(const double x[4], const double y[4], double W[4]);

#ifdef __cplusplus
}
#endif
#endif
