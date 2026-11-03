/* rcheck.c -- probe for the rarity-sweep efficiency audit (Ref4 item 5).
 *
 * Replays the kernel stream of PRFull/cannon/sweep13 (census.c psweep, alpha = 1/3, seed 3000000+i), and for every
 * draw reads the self-play class twice:
 *   read  : from the weights, exactly as census.c selfclass() (thresholds 1e-12);
 *   exact : from the support Sigma = states carrying the least exponent (no tolerance; L1 of SI 'Memory two').
 * A draw is flagged when the two differ.  Flagged draws then get the psweep competitiveness test
 * (screen + rival_margin, tolerance -1e-9), to see which friendly-rival count they entered.
 *
 * usage: rcheck <task> <ndraws> [classify]
 *   prints: RC task n effC effA effZ | exactC exactA exactZ | flagged | contaminated W S N E
 *   and one FLAG line per flagged draw.
 * kdraw() and the seeding are copied verbatim from sweep13/census.c; freqs()/gth() from sweep13/m2f.c,
 * with gth extended to return the exponent of each state.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "witnesses.h"
#define NS 16
#define EINF 1000000000

static inline int swapo(int o){ return ((o&1)<<1)|(o>>1); }
static inline int mirror_state(int j){ return 4*swapo(j>>2) + swapo(j&3); }
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
    double cs[2]={sC,sD}; int es[2]={esC,esD};
    double ct[2]={tC,tD}; int et[2]={etC,etD};
    for(int dF=0;dF<2;dF++) for(int dA=0;dA<2;dA++){
      int onew=2*dF+dA, nxt=4*onew+(j>>2);
      double cc=cs[dF]*ct[dA]; int ee=es[dF]+et[dA];
      if(ee<e[j][nxt]){ e[j][nxt]=ee; c[j][nxt]=cc; }
      else if(ee==e[j][nxt]) c[j][nxt]+=cc;
    }
  }
}
/* as m2f.c gth, plus insig[k] = 1 iff state k carries the least exponent */
static void gth(double c[NS][NS],int e[NS][NS],double w[NS],int insig[NS]){
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
  for(int k=0;k<NS;k++){ insig[k]=(pe[k]==emin); w[k]=insig[k]?pc[k]/ct:0.0; }
}
static void freqs(const double*p,const double*q,double W[4]){
  double c[NS][NS]; int e[NS][NS]; double w[NS]; int s[NS];
  build(p,q,c,e); gth(c,e,w,s);
  W[0]=W[1]=W[2]=W[3]=0.0; for(int j=0;j<NS;j++) W[j>>2]+=w[j];
}
static void genome_to_p(unsigned g,double*p){ for(int j=0;j<NS;j++) p[j]=(double)((g>>j)&1); }
static double rival_margin(const double*p,int sign,double earlyquit){
  double best=1e18; double q[NS],W[4];
  for(long t=0;t<65536;t++){ genome_to_p((unsigned)t,q); freqs(p,q,W);
    double d=sign*(W[2]-W[1]); if(d<best){ best=d; if(best<earlyquit) break; } }
  return best;
}
static double WP[NWIT][NS];
static int screen(const double*p,int sign){ double W[4];
  for(int k=0;k<NWIT;k++){ freqs(p,WP[k],W); if(sign*(W[2]-W[1]) < -1e-9) return 0; } return 1; }

/* ---- kernel, verbatim from sweep13/census.c (rational form only) */
static uint64_t rs;
static inline double urnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17;
  return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }
static inline uint64_t k53(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (rs>>11)&0x1FFFFFFFFFFFFFull; }
static uint64_t AN=1, AD=3, AQ; static long NREDRAW;
static inline void kdraw(double*p,long*npin){
  for(int j=0;j<NS;j++){
    uint64_t k; while((k=k53())>=AD*AQ) NREDRAW++;
    if(k<AN*AQ){ p[j]=0.0; (*npin)++; } else if(k<2*AN*AQ){ p[j]=1.0; (*npin)++; }
    else { while((p[j]=urnd())==0.0) NREDRAW++; } } }

int main(int argc,char**argv){
  if(argc>1&&!strcmp(argv[1],"margin")){ /* read 16 coordinates per line from stdin; full scans, no early quit */
    for(int k=0;k<NWIT;k++) for(int j=0;j<NS;j++) WP[k][j]=(double)((WITG[k]>>j)&1);
    double p[NS]; while(1){ for(int j=0;j<NS;j++) if(scanf("%lf",&p[j])!=1) return 0;
      printf("MARGIN\tP %.6e\tM %.6e\tscreenP %d\tscreenM %d\n",rival_margin(p,+1,-1e300),rival_margin(p,-1,-1e300),screen(p,+1),screen(p,-1)); } }
  if(argc<3){ fprintf(stderr,"usage: rcheck <task> <ndraws> [classify] | rcheck margin < coords\n"); return 2; }
  long task=atol(argv[1]), n=atol(argv[2]); int classify=argc>3;
  for(int k=0;k<NWIT;k++) for(int j=0;j<NS;j++) WP[k][j]=(double)((WITG[k]>>j)&1);
  AQ=(1ull<<53)/AD;
  rs=(uint64_t)(3000000+task)*2862933555777941757ULL+3037000493ULL; if(!rs)rs=7;
  long npin=0, rC=0,rA=0,rZ=0, xC=0,xA=0,xZ=0, nflag=0, cW=0,cS=0,cN=0,cE=0;
  double p[NS];
  for(long i=0;i<n;i++){
    kdraw(p,&npin);
    double c[NS][NS]; int e[NS][NS]; double w[NS]; int sg[NS];
    build(p,p,c,e); gth(c,e,w,sg);
    double W[4]={0,0,0,0}; for(int j=0;j<NS;j++) W[j>>2]+=w[j];
    char rd = W[0]>1-1e-12?'C':(W[1]+W[2]>1-1e-12?'A':(W[3]<1e-12?'Z':'-'));
    int hasO[4]={0,0,0,0}; double minw=2; for(int j=0;j<NS;j++) if(sg[j]){ hasO[j>>2]=1; if(w[j]<minw)minw=w[j]; }
    char ex = (!hasO[1]&&!hasO[2]&&!hasO[3])?'C':((!hasO[0]&&!hasO[3])?'A':(!hasO[3]?'Z':'-'));
    if(rd=='C')rC++; else if(rd=='A')rA++; else if(rd=='Z')rZ++;
    if(ex=='C')xC++; else if(ex=='A')xA++; else if(ex=='Z')xZ++;
    if(rd!=ex){ nflag++;
      int isP=-1,isM=-1;
      if(classify){ isP=screen(p,+1)&&rival_margin(p,+1,-1e-9)>=-1e-9; isM=screen(p,-1)&&rival_margin(p,-1,-1e-9)>=-1e-9;
        if(rd=='C'&&isP)cW++; if(rd=='C'&&isM)cS++; if(rd=='A'&&isP)cN++; if(rd=='A'&&isM)cE++; }
      printf("FLAG\t%ld\t%ld\t%c\t%c\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%d\t%d",task,i,rd,ex,minw,W[0],W[1],W[2],W[3],isP,isM);
      for(int j=0;j<NS;j++) printf("\t%.17g",p[j]); printf("\n"); }
  }
  printf("RC\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",task,n,rC,rA,rZ,xC,xA,xZ,nflag,cW,cS,cN,cE);
  return 0;
}
