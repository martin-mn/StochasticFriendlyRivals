/* m2f.c -- exact eps->0 limit of the stationary distribution of a pair of memory-two
 * strategies drawn from the FULL continuum [0,1]^16, with independent implementation
 * errors p -> eps + (1-2eps)p on both players (convention A).
 *
 * Continuum generalisation of ../../PartnersRivals/check/pairs.c.  Leading-order
 * (coefficient, exponent-of-eps) arithmetic:
 *     Pr(focal C at state j) = eps + (1-2eps) p[j]  ->  (p[j], 0)   if p[j] > 0
 *                                                       (1,     1)   if p[j] = 0
 *     Pr(focal D at state j) = 1 - that              ->  (1-p[j], 0) if p[j] < 1
 *                                                       (1,      1)  if p[j] = 1
 * Every intermediate of GTH state reduction is a sum/product/quotient of positive
 * quantities, so the leading term is exact: no threshold, no extrapolation.
 *
 * Convention (NE.pdf / PartnersRivals): state j = 4*recent + before, outcome
 * o = 2*(focal defects) + (co-player defects), CC=0 CD=1 DC=2 DD=3, own action first.
 * Bit j of a 16-bit genome is 1 for C, so ALLC = 65535, ALLD = 0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NS 16
#define EINF 1000000000

static inline int swapo(int o){ return ((o&1)<<1)|(o>>1); }          /* CD <-> DC */
static inline int mirror_state(int j){ return 4*swapo(j>>2) + swapo(j&3); }

/* leading order of Pr(C) and Pr(D) for a probability x in [0,1] */
static inline void lead(double x,double*cC,int*eC,double*cD,int*eD){
  if(x>0.0){ *cC=x;   *eC=0; } else { *cC=1.0; *eC=1; }
  if(x<1.0){ *cD=1-x; *eD=0; } else { *cD=1.0; *eD=1; }
}

static void build(const double*p,const double*q,double c[NS][NS],int e[NS][NS]){
  for(int i=0;i<NS;i++) for(int j=0;j<NS;j++){ c[i][j]=0.0; e[i][j]=EINF; }
  for(int j=0;j<NS;j++){
    int jt=mirror_state(j);
    double sC,sD,tC,tD; int esC,esD,etC,etD;
    lead(p[j],&sC,&esC,&sD,&esD);
    lead(q[jt],&tC,&etC,&tD,&etD);
    double cs[2]={sC,sD}; int es[2]={esC,esD};      /* index 0 = focal C, 1 = focal D */
    double ct[2]={tC,tD}; int et[2]={etC,etD};
    for(int dF=0;dF<2;dF++) for(int dA=0;dA<2;dA++){
      int onew=2*dF+dA, nxt=4*onew+(j>>2);
      double cc=cs[dF]*ct[dA]; int ee=es[dF]+et[dA];
      if(ee<e[j][nxt]){ e[j][nxt]=ee; c[j][nxt]=cc; }
      else if(ee==e[j][nxt]) c[j][nxt]+=cc;
    }
  }
}

static void gth(double c[NS][NS],int e[NS][NS],double w[NS]){
  for(int k=NS-1;k>=1;k--){
    int es=EINF; double cs=0;
    for(int j=0;j<k;j++){ int ej=e[k][j];
      if(ej<es){es=ej;cs=c[k][j];} else if(ej==es&&es<EINF){cs+=c[k][j];} }
    for(int i=0;i<k;i++){ if(e[i][k]<EINF){ c[i][k]/=cs; e[i][k]-=es; } }
    for(int i=0;i<k;i++){ if(e[i][k]>=EINF) continue; int eik=e[i][k]; double cik=c[i][k];
      for(int j=0;j<k;j++){ if(e[k][j]>=EINF) continue; int en=eik+e[k][j]; double cn=cik*c[k][j];
        if(en<e[i][j]){e[i][j]=en;c[i][j]=cn;} else if(en==e[i][j]){c[i][j]+=cn;} } }
  }
  double pc[NS]; int pe[NS]; pc[0]=1; pe[0]=0;
  for(int k=1;k<NS;k++){ int ek=EINF; double ck=0;
    for(int i=0;i<k;i++){ if(e[i][k]>=EINF||pe[i]>=EINF) continue; int en=pe[i]+e[i][k]; double cn=pc[i]*c[i][k];
      if(en<ek){ek=en;ck=cn;} else if(en==ek){ck+=cn;} }
    pc[k]=ck; pe[k]=ek; }
  int emin=EINF; double ct=0;
  for(int k=0;k<NS;k++){ if(pe[k]<emin){emin=pe[k];ct=pc[k];} else if(pe[k]==emin){ct+=pc[k];} }
  for(int k=0;k<NS;k++) w[k]=(pe[k]==emin)?pc[k]/ct:0.0;
}

/* limit outcome frequencies from the focal player's view: W[0]=CC W[1]=CD W[2]=DC W[3]=DD */
void freqs(const double*p,const double*q,double W[4]){
  double c[NS][NS]; int e[NS][NS]; double w[NS];
  build(p,q,c,e); gth(c,e,w);
  W[0]=W[1]=W[2]=W[3]=0.0;
  for(int j=0;j<NS;j++) W[j>>2]+=w[j];
}
void freqs16(const double*p,const double*q,double w[NS]){
  double c[NS][NS]; int e[NS][NS];
  build(p,q,c,e); gth(c,e,w);
}

static void genome_to_p(unsigned g,double*p){ for(int j=0;j<NS;j++) p[j]=(double)((g>>j)&1); }

/* ------------------------------------------------------------------ rival margin
 * sign=+1 : T>S reading, p is a rival iff  W_DC - W_CD >= 0  against every co-player.
 * sign=-1 : T<S reading, iff W_CD - W_DC >= 0.
 * By the box-vertex lemma (see Text/rivals.md) the worst co-player over the whole
 * memory-two continuum is attained at a binary co-player, so the min over the 65536
 * binary genomes is the exact worst case over q in [0,1]^16.
 * Returns the minimum of sign*(W_DC - W_CD); *arg receives the arg-min genome. */
double rival_margin(const double*p,int sign,long*arg,double earlyquit){
  double best=1e18; long ba=-1; double q[NS],W[4];
  for(long t=0;t<65536;t++){
    genome_to_p((unsigned)t,q);
    freqs(p,q,W);
    double d=sign*(W[2]-W[1]);
    if(d<best){ best=d; ba=t; if(best<earlyquit) break; }
  }
  if(arg)*arg=ba;
  return best;
}
