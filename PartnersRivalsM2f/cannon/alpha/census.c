/* census.c -- the definitive census of FRIENDLY RIVALS over the full memory-two continuum.
 *
 * Two products, selected by argv[1]:
 *
 *  strata <stripe> <nstripes> <Keff> <out> <pts> [Kother]
 *      For every one of the 3^16 = 43 046 721 support strata (each coordinate p[j] pinned to 0,
 *      pinned to 1, or free in (0,1)) draw K random interior points and decide, EXACTLY, for
 *      each: the eps->0 self-play class, and the two rival margins
 *          m+ = min over all 65536 binary memory-two co-players of (W_DC - W_CD)
 *          m- = min over the same of (W_CD - W_DC).
 *      Screened first against WITG (no false negatives), then confirmed against all 65536.
 *      Emits one line per stratum with any confirmed rival:  g cls nP nM
 *      and, for strata that are ALSO efficient, one exemplar point per confirmed draw.
 *
 *  psweep <n> <seed> <alpha_permille> <out>
 *      Predicate-free, PURELY STRATIFIED at MUTATIONAL BIAS alpha = alpha_permille/1000: each of
 *      the sixteen coordinates is set to 0 with probability alpha, to 1 with probability alpha,
 *      and drawn uniformly on (0,1) with probability 1-2alpha.  n such strategies, fully
 *      classified.
 *
 *      alpha is the parameter of Nowak et al. (Literature/Nowak/WF.pdf, SI 3.1), where the same
 *      three-outcome draw defines the mutational bias of the Wright-Fisher process: alpha = 0 is
 *      unbiased/uniform, alpha = 1/2 gives binary strategies, and alpha = 1/3 is that paper's
 *      working value.  One difference: there p_min = eps and p_max = 1-eps at a fixed noise floor
 *      eps = 1e-4, whereas here the pinned values are exactly 0 and 1 because this engine takes
 *      the eps -> 0 limit analytically rather than simulating at a floor.
 *
 *      This replaces the earlier `sweep`, whose mixture wasted half its draws on the uniform law
 *      -- the friendly-rival set is Lebesgue-null (efficiency needs a pinned coordinate, rivalry
 *      needs p15=0), so a uniform draw is a friendly rival with probability exactly zero -- and
 *      which tallied both arms into ONE set of counters, so the uniform arm's null could not be
 *      read off the output.  Here there is no uniform arm and pi is explicit.
 *
 *      alpha = 0.500 pins every coordinate, so the draw is uniform on the 65536 binary vertices and
 *      the eleven rates must converge to the EXACT deterministic census.  That is a proved
 *      endpoint and hence an end-to-end validation:
 *          FR_W 8/65536   FR_S 1519/65536   FR_N 80/65536   FR_E 80/65536
 *          effCC 7639/65536   effALT 3072/65536   rivP 2640/65536   rivM 2640/65536
 *      alpha -> 0 gives all-interior draws, where every one of the eleven counts must be 0.
 *
 *      Purpose: the absolute rates are a property of alpha, not of the set.  What the manuscript
 *      claims to be alpha-independent is the COMPARISON -- wedge against wedge, and the joint
 *      against the product of the marginals.  This sweep measures that claim instead of
 *      asserting it.
 *
 * Wedge bookkeeping (base game (R,S,T,P) = (1,u,1+v,0)):
 *      below the switch line u+v<1 efficiency is W_CC = 1;  above it is W_CD = W_DC = 1/2;
 *      on it is W_DD = 0.   T>S is u-v<1, T<S is u-v>1.
 *      W = effCC & riv+   S = effCC & riv-   N = effALT & riv+   E = effALT & riv-
 *      switch-line rays use effLINE = (W_DD = 0).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "witnesses.h"
#define NS 16
void freqs(const double*p,const double*q,double W[4]);
double rival_margin(const double*p,int sign,long*arg,double earlyquit);

static double WP[NWIT][NS];
static uint64_t rs;
static inline double urnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17;
  return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }
static void initwit(void){ for(int k=0;k<NWIT;k++) for(int j=0;j<NS;j++) WP[k][j]=(double)((WITG[k]>>j)&1); }

/* returns 1 if p survives the screen for this sign (i.e. might still be a rival) */
static int screen(const double*p,int sign){
  double W[4];
  for(int k=0;k<NWIT;k++){ freqs(p,WP[k],W); if(sign*(W[2]-W[1]) < -1e-9) return 0; }
  return 1;
}
/* 0 = not efficient anywhere, 'C','A','Z' as in the header */
static char selfclass(const double*p){
  double W[4]; freqs(p,p,W);
  if(W[0]>1-1e-12) return 'C';
  if(W[1]+W[2]>1-1e-12) return 'A';
  if(W[3]<1e-12) return 'Z';
  return '-';
}
static void draw(const int*dg,double*p){
  for(int j=0;j<NS;j++) p[j] = dg[j]==0?0.0 : (dg[j]==2?1.0 : 0.02+0.96*urnd());
}

int main(int argc,char**argv){
  initwit();
  if(argc>1 && !strcmp(argv[1],"strata")){
    long stripe=atol(argv[2]), nstr=atol(argv[3]);
    int KEFF=atoi(argv[4]);            /* draws per stratum when the stratum is efficient */
    int KOTH=argc>7?atoi(argv[7]):8;   /* draws per stratum otherwise */
    FILE*fo=fopen(argv[5],"w"); FILE*fp=fopen(argv[6],"w");
    long TOT=1; for(int i=0;i<16;i++) TOT*=3;
    long nfull=0, nemit=0; double p[NS];  int dg[16];
    for(long g=stripe;g<TOT;g+=nstr){
      long x=g; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; }
      rs=(uint64_t)g*6364136223846793005ULL+1442695040888963407ULL; if(!rs) rs=88172645463325252ULL;
      draw(dg,p); char cls=selfclass(p);
      int K = (cls=='-')?KOTH:KEFF;
      int nP=0,nM=0, dumpP=0,dumpM=0;
      double bestP=-1e18,bestM=-1e18;
      for(int r=0;r<K;r++){
        draw(dg,p);
        long a; double m;
        if(screen(p,+1)){ nfull++; m=rival_margin(p,+1,&a,-1e-9);
          if(m>bestP)bestP=m;
          if(m>=-1e-9){ nP++;
            if(cls!='-' && dumpP<4){ dumpP++;
              fprintf(fp,"P %ld %c",g,cls); for(int j=0;j<NS;j++) fprintf(fp," %.17g",p[j]); fprintf(fp,"\n"); } } }
        if(screen(p,-1)){ nfull++; m=rival_margin(p,-1,&a,-1e-9);
          if(m>bestM)bestM=m;
          if(m>=-1e-9){ nM++;
            if(cls!='-' && dumpM<4){ dumpM++;
              fprintf(fp,"M %ld %c",g,cls); for(int j=0;j<NS;j++) fprintf(fp," %.17g",p[j]); fprintf(fp,"\n"); } } }
      }
      if(nP||nM){ nemit++; fprintf(fo,"%ld %c %d %d %d\n",g,cls,nP,nM,K); }
    }
    fprintf(stderr,"stripe %ld/%ld KEFF=%d KOTH=%d: emitted %ld strata, %ld full 65536-scans\n",
            stripe,nstr,KEFF,KOTH,nemit,nfull);
    fclose(fo); fclose(fp);
  } else if(argc>1 && !strcmp(argv[1],"psweep")){
    long n=atol(argv[2]); rs=(uint64_t)atol(argv[3])*2862933555777941757ULL+3037000493ULL; if(!rs)rs=7;
    long am=atol(argv[4]);
    if(am<0||am>500){ fprintf(stderr,"alpha_permille must be in [0,500]\n"); return 2; }
    const double al=am/1000.0;
    FILE*fo=fopen(argv[5],"w");
    long npin=0;                                        /* total pinned coordinates, for a check */
    long cnt[4][4]; memset(cnt,0,sizeof cnt);           /* [class C/A/Z/-][sign none/P/M/both] */
    long frW=0,frS=0,frN=0,frE=0,frRayP=0,frRayM=0,effC=0,effA=0,effZ=0,rivP=0,rivM=0;
    double p[NS];
    for(long i=0;i<n;i++){
      /* purely stratified at alpha: 0 w.p. alpha, 1 w.p. alpha, uniform w.p. 1-2alpha */
      for(int j=0;j<NS;j++){ double u=urnd();
        if(u<al){ p[j]=0.0; npin++; } else if(u<2.0*al){ p[j]=1.0; npin++; } else p[j]=urnd(); }
      char cls=selfclass(p);
      if(cls=='C')effC++; else if(cls=='A')effA++; else if(cls=='Z')effZ++;
      int isP=0,isM=0; long a;
      if(screen(p,+1)&&rival_margin(p,+1,&a,-1e-9)>=-1e-9) isP=1;
      if(screen(p,-1)&&rival_margin(p,-1,&a,-1e-9)>=-1e-9) isM=1;
      rivP+=isP; rivM+=isM;
      if(cls=='C'&&isP)frW++;
      if(cls=='C'&&isM)frS++;
      if(cls=='A'&&isP)frN++;
      if(cls=='A'&&isM)frE++;
      if((cls=='C'||cls=='A'||cls=='Z')&&isP)frRayP++;
      if((cls=='C'||cls=='A'||cls=='Z')&&isM)frRayM++;
      int ci = cls=='C'?0:(cls=='A'?1:(cls=='Z'?2:3));
      cnt[ci][isP+2*isM]++;
      if(frW && isP && cls=='C'){ fprintf(fo,"W"); for(int j=0;j<NS;j++) fprintf(fo," %.17g",p[j]); fprintf(fo,"\n"); }
    }
    /* one machine-readable line; alpha first so the collector can group by it */
    fprintf(stderr,"psweep alpha=%.3f n=%ld: effCC %ld effALT %ld effZ %ld rivP %ld rivM %ld | "
                   "FR_W %ld FR_S %ld FR_N %ld FR_E %ld | rayP %ld rayM %ld | meanpin %.6f\n",
            al,n,effC,effA,effZ,rivP,rivM,frW,frS,frN,frE,frRayP,frRayM,
            (double)npin/((double)n*NS));
    printf("ALPHA\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
            am,n,effC,effA,effZ,rivP,rivM,frW,frS,frN,frE,frRayP,frRayM);
    fclose(fo);
  } else { fprintf(stderr,"usage: census strata|psweep ...\n"); return 1; }
  return 0;
}
