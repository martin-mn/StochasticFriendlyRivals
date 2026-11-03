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
 *  psweep <n> <seed> <alpha> <out>
 *      <alpha> is an integer permille ("333": alpha = 333/1000.0, a double threshold on the double
 *      urnd(), exactly as in every archived run, which it reproduces bit for bit) or an EXACT
 *      rational "num/den" ("1/3"; PRFull sweep13 patch, 2026-09-27).  In the rational form the pin
 *      decision uses no floating point: with k the 53-bit integer behind urnd() (urnd = k/2^53) and
 *      q = floor(2^53/den), k is redrawn while k >= den*q, and then k < num*q pins to 0 and
 *      k < 2*num*q pins to 1, each with probability num*q/(den*q) = num/den exactly.  The redraw
 *      has probability (2^53 mod den)/2^53, 2.2e-16 at den = 3, so the stream consumed is, to that
 *      probability, the permille path's; the interior value is urnd() redrawn while it is 0.0
 *      (probability 2^-53), so it lies strictly in (0,1) and a 0 is only ever a pin.
 *      The ALPHA line carries <alpha> as given ("333" or "1/3").
 *
 *  pdraw <n> <seed> <alpha>
 *      The kernel alone, no classification: n draws from the same code path as psweep with the
 *      same seeding, tallied per coordinate (0 / 1 / interior), for the law check.
 *
 *  psweep, continued.
 *      Predicate-free, PURELY STRATIFIED at MUTATIONAL BIAS alpha: each of
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

/* ---- the mutation kernel (sweep13 patch).  urnd() above is untouched; k53() takes the same xorshift
 * step and returns the integer urnd() would have divided by 2^53. */
static inline uint64_t k53(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (rs>>11)&0x1FFFFFFFFFFFFFull; }
static double AL;                      /* permille form: alpha = am/1000.0, the archived double threshold */
static uint64_t AN, AD, AQ;            /* rational form: alpha = AN/AD, AQ = floor(2^53/AD); AD = 0 means permille */
static long NREDRAW;                   /* rational form: redraws of k (k >= AD*AQ) and of an interior 0.0 */
/* parse "333" (permille, 0..500) or "num/den" (0 <= 2 num <= den <= 2^53); returns 0 on success */
static int parse_alpha(const char*s,long*am,char*tag,size_t ntag){
  const char*sl=strchr(s,'/'); char*e;
  if(!sl){ *am=strtol(s,&e,10); if(e==s||*e||*am<0||*am>500) return 1;
    AD=0; AL=*am/1000.0; snprintf(tag,ntag,"%ld",*am); return 0; }
  AN=strtoull(s,&e,10); if(e!=sl) return 1;
  AD=strtoull(sl+1,&e,10); if(e==sl+1||*e||AD==0||AD>(1ull<<53)||2*AN>AD) return 1;
  AQ=(1ull<<53)/AD; AL=(double)AN/(double)AD; *am=-1; snprintf(tag,ntag,"%llu/%llu",(unsigned long long)AN,(unsigned long long)AD);
  return 0; }
/* one strategy from the kernel: 0 w.p. alpha, 1 w.p. alpha, uniform on (0,1) otherwise */
static inline void kdraw(double*p,long*npin){
  for(int j=0;j<NS;j++){
    if(AD){ uint64_t k; while((k=k53())>=AD*AQ) NREDRAW++;         /* exact rational: integers only */
      if(k<AN*AQ){ p[j]=0.0; (*npin)++; } else if(k<2*AN*AQ){ p[j]=1.0; (*npin)++; }
      else { while((p[j]=urnd())==0.0) NREDRAW++; } }
    else { double u=urnd();                                          /* permille: exactly as archived */
      if(u<AL){ p[j]=0.0; (*npin)++; } else if(u<2.0*AL){ p[j]=1.0; (*npin)++; } else p[j]=urnd(); } } }

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
    long am; char atag[48], ahead[48];
    if(argc<6||parse_alpha(argv[4],&am,atag,sizeof atag)){
      fprintf(stderr,"alpha must be a permille in [0,500] or a rational num/den with 2 num <= den\n"); return 2; }
    const double al=AL;
    if(AD) snprintf(ahead,sizeof ahead,"%s",atag); else snprintf(ahead,sizeof ahead,"%.3f",al);
    FILE*fo=fopen(argv[5],"w");
    long npin=0;                                        /* total pinned coordinates, for a check */
    long cnt[4][4]; memset(cnt,0,sizeof cnt);           /* [class C/A/Z/-][sign none/P/M/both] */
    long frW=0,frS=0,frN=0,frE=0,frRayP=0,frRayM=0,effC=0,effA=0,effZ=0,rivP=0,rivM=0;
    double p[NS];
    for(long i=0;i<n;i++){
      /* purely stratified at alpha: 0 w.p. alpha, 1 w.p. alpha, uniform w.p. 1-2alpha */
      kdraw(p,&npin);
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
    fprintf(stderr,"psweep alpha=%s n=%ld: effCC %ld effALT %ld effZ %ld rivP %ld rivM %ld | "
                   "FR_W %ld FR_S %ld FR_N %ld FR_E %ld | rayP %ld rayM %ld | meanpin %.6f\n",
            ahead,n,effC,effA,effZ,rivP,rivM,frW,frS,frN,frE,frRayP,frRayM,
            (double)npin/((double)n*NS));
    if(AD) fprintf(stderr,"psweep exact alpha=%s: %ld redraws\n",atag,NREDRAW);
    printf("ALPHA\t%s\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",
            atag,n,effC,effA,effZ,rivP,rivM,frW,frS,frN,frE,frRayP,frRayM);
    fclose(fo);
  } else if(argc>4 && !strcmp(argv[1],"pdraw")){
    /* the kernel alone (sweep13 patch): same seeding, same kdraw() as psweep, no classification */
    long n=atol(argv[2]); rs=(uint64_t)atol(argv[3])*2862933555777941757ULL+3037000493ULL; if(!rs)rs=7;
    long am; char atag[48];
    if(parse_alpha(argv[4],&am,atag,sizeof atag)){ fprintf(stderr,"bad alpha\n"); return 2; }
    long c0[NS]={0},c1[NS]={0},cu[NS]={0},npin=0; double p[NS], umin=1.0, umax=0.0, usum=0.0;
    for(long i=0;i<n;i++){ kdraw(p,&npin);
      for(int j=0;j<NS;j++){ if(p[j]==0.0) c0[j]++; else if(p[j]==1.0) c1[j]++;
        else { cu[j]++; usum+=p[j]; if(p[j]<umin)umin=p[j]; if(p[j]>umax)umax=p[j]; } } }
    long t0=0,t1=0,tu=0; for(int j=0;j<NS;j++){ t0+=c0[j]; t1+=c1[j]; tu+=cu[j]; }
    printf("PDRAW\t%s\t%ld\t%ld\t%ld\t%ld\t%ld\t%.17g\t%.17g\t%.17g\n",atag,n,t0,t1,tu,NREDRAW,
           tu?usum/tu:0.0,umin,umax);
    for(int j=0;j<NS;j++) printf("COORD\t%d\t%ld\t%ld\t%ld\n",j,c0[j],c1[j],cu[j]);
  } else { fprintf(stderr,"usage: census strata|psweep|pdraw ...\n"); return 1; }
  return 0;
}
