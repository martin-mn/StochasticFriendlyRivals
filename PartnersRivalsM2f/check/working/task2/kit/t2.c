/* t2.c -- task-2 verification kit.  Three modes.
 *
 *  ncensus <sign> <lo> <hi> <K1> <K2> <seed>
 *     Exhaustive re-census of one alternating wedge.  For sign=+1 (wedge N) the PROVED
 *     necessary condition Riv+ subset {p15 = 0} restricts the candidates to d15 = 0, i.e.
 *     g < 3^15; for sign=-1 (wedge E) its mirror image gives d0 = 2, i.e. g % 3 == 2.
 *     For each candidate g in [lo,hi): classify the eps->0 self-play limit from one draw
 *     (support-determined); if it is the alternating class, take K1 draws with the full
 *     screen + 65536-co-player confirmation and record h1 = hits.  If h1 == 0, take K2
 *     further draws and record h2.  Emit "g nf h1 h2 K1 K2" for every g with h1+h2 > 0.
 *
 *  list <file> <shard> <nshard> <seed>
 *     Lines "g sign K"; process those whose 0-based line index % nshard == shard.
 *     Emit "g sign K nA nC nZ nX scr hit".
 *
 *  seedtest <g> <sign> <s0> <ns> <K>
 *     census.c's own seeding rs = s*6364136223846793005 + 1442695040888963407 for s in
 *     [s0,s0+ns), one class draw then K draws; emit "s hits" per seed.  Tests whether the
 *     first K outputs of that seeded xorshift64 stream behave like K iid draws.
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
static inline uint64_t sm(uint64_t x){ x+=0x9E3779B97F4A7C15ULL;
  x=(x^(x>>30))*0xBF58476D1CE4E5B9ULL; x=(x^(x>>27))*0x94D049BB133111EBULL; return x^(x>>31); }
static int screen(const double*p,int sign){ double W[4];
  for(int k=0;k<NWIT;k++){ freqs(p,WP[k],W); if(sign*(W[2]-W[1])<-1e-9) return 0; } return 1; }
static void initwit(void){ for(int k=0;k<NWIT;k++) for(int j=0;j<NS;j++) WP[k][j]=(double)((WITG[k]>>j)&1); }
static inline void drw(const int*dg,double*p){
  for(int j=0;j<NS;j++){ double u=urnd(); p[j]= dg[j]==0?0.0:(dg[j]==2?1.0:0.02+0.96*u); } }
static inline char cls_of(const double*p){ double W[4]; freqs(p,p,W);
  return W[0]>1-1e-12?'C':(W[1]+W[2]>1-1e-12?'A':(W[3]<1e-12?'Z':'-')); }

int main(int argc,char**argv){
  initwit();
  if(argc>1 && !strcmp(argv[1],"ncensus")){
    int sign=atoi(argv[2]); long lo=atol(argv[3]), hi=atol(argv[4]);
    int K1=atoi(argv[5]), K2=atoi(argv[6]); uint64_t seed=strtoull(argv[7],0,10);
    long nA=0,nemit=0,ndeep=0; double p[NS]; int dg[16];
    for(long g=lo;g<hi;g++){
      long x=g; int nf=0; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; if(dg[j]==1)nf++; }
      if(sign>0){ if(dg[15]!=0) continue; } else { if(dg[0]!=2) continue; }
      rs=sm((uint64_t)g*0x2545F4914F6CDD1DULL ^ sm(seed)); if(!rs) rs=0x9E3779B97F4A7C15ULL;
      drw(dg,p); if(cls_of(p)!='A') continue;
      nA++;
      long a; double m; int h1=0,h2=0;
      for(int r=0;r<K1;r++){ drw(dg,p);
        if(screen(p,sign)){ m=rival_margin(p,sign,&a,-1e-9); if(m>=-1e-9) h1++; } }
      if(h1==0){ ndeep++;
        for(int r=0;r<K2;r++){ drw(dg,p);
          if(screen(p,sign)){ m=rival_margin(p,sign,&a,-1e-9); if(m>=-1e-9){ h2++; if(h2>=4) break; } } } }
      if(h1||h2){ nemit++; printf("%ld %d %d %d %d %d\n",g,nf,h1,h2,K1,K2); }
    }
    fprintf(stderr,"ncensus sign=%d [%ld,%ld) alt-strata %ld deep %ld emitted %ld\n",sign,lo,hi,nA,ndeep,nemit);
  } else if(argc>1 && !strcmp(argv[1],"list")){
    FILE*f=fopen(argv[2],"r"); long shard=atol(argv[3]), nsh=atol(argv[4]);
    uint64_t seed=strtoull(argv[5],0,10);
    char line[128]; long idx=0; double p[NS]; int dg[16];
    while(fgets(line,sizeof line,f)){
      long g; int sign,K;
      if(sscanf(line,"%ld %d %d",&g,&sign,&K)!=3) continue;
      if((idx++)%nsh!=shard) continue;
      long x=g; int nf=0; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; if(dg[j]==1)nf++; }
      rs=sm((uint64_t)g*0x9E3779B97F4A7C15ULL ^ sm(seed*1000003+sign)); if(!rs) rs=12345;
      long nA=0,nC=0,nZ=0,nX=0,scr=0,hit=0;
      for(int r=0;r<K;r++){
        drw(dg,p); char c=cls_of(p);
        if(c=='A')nA++; else if(c=='C')nC++; else if(c=='Z')nZ++; else nX++;
        if(screen(p,sign)){ scr++; long a; double m=rival_margin(p,sign,&a,-1e-9); if(m>=-1e-9) hit++; }
      }
      printf("%ld %d %d %ld %ld %ld %ld %ld %ld\n",g,sign,K,nA,nC,nZ,nX,scr,hit);
      fflush(stdout);
    }
    fclose(f);
  } else if(argc>1 && !strcmp(argv[1],"seedtest")){
    long g=atol(argv[2]); int sign=atoi(argv[3]);
    long s0=atol(argv[4]), ns=atol(argv[5]); int K=atoi(argv[6]);
    int dg[16]; long x=g; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; }
    double p[NS];
    for(long s=s0;s<s0+ns;s++){
      rs=(uint64_t)s*6364136223846793005ULL+1442695040888963407ULL; if(!rs) rs=88172645463325252ULL;
      drw(dg,p);                                   /* census.c's class draw */
      int hit=0;
      for(int r=0;r<K;r++){ drw(dg,p);
        if(screen(p,sign)){ long a; double m=rival_margin(p,sign,&a,-1e-9); if(m>=-1e-9) hit++; } }
      printf("%ld %d\n",s,hit);
    }
  } else { fprintf(stderr,"usage: t2 ncensus|list|seedtest ...\n"); return 1; }
  return 0;
}
