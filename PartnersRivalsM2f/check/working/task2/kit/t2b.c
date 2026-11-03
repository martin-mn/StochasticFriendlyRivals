/* t2b.c -- exhaustive re-census of the two BELOW-the-switch-line wedges.
 *
 *   ccensus <sign> <lo> <hi> <K1> <K2> <seed>
 *     sign=+1 -> wedge W, sign=-1 -> wedge S.  Efficiency below the switch line forces
 *     p0 = 1 (proved, friendly.md 8.1), so candidates have d0 = 2; for sign=+1 the proved
 *     Riv+ subset {p15=0} adds d15 = 0.  Self-play class must be 'C'.
 *     Stage 1: up to K1 draws, early stop once 8 hits (r1 = draws actually used).
 *     Stage 2, only if h1 == 0: up to K2 draws, early stop at 4 hits.
 *     Emits "g nf h1 r1 h2 K2" for every stratum with h1+h2 > 0.
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
static inline void drw(const int*dg,double*p){
  for(int j=0;j<NS;j++){ double u=urnd(); p[j]= dg[j]==0?0.0:(dg[j]==2?1.0:0.02+0.96*u); } }
int main(int argc,char**argv){
  for(int k=0;k<NWIT;k++) for(int j=0;j<NS;j++) WP[k][j]=(double)((WITG[k]>>j)&1);
  if(argc<8||strcmp(argv[1],"ccensus")){ fprintf(stderr,"usage: t2b ccensus sign lo hi K1 K2 seed\n"); return 1; }
  int sign=atoi(argv[2]); long lo=atol(argv[3]), hi=atol(argv[4]);
  int K1=atoi(argv[5]), K2=atoi(argv[6]); uint64_t seed=strtoull(argv[7],0,10);
  long nC=0,ndeep=0,nemit=0; double p[NS],W[4]; int dg[16];
  for(long g=lo;g<hi;g++){
    long x=g; int nf=0; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; if(dg[j]==1)nf++; }
    if(dg[0]!=2) continue;
    if(sign>0 && dg[15]!=0) continue;
    rs=sm((uint64_t)g*0x2545F4914F6CDD1DULL ^ sm(seed+7)); if(!rs) rs=0x9E3779B97F4A7C15ULL;
    drw(dg,p); freqs(p,p,W); if(!(W[0]>1-1e-12)) continue;
    nC++;
    long a; double m; int h1=0,h2=0,r1=0;
    for(int r=0;r<K1;r++){ r1=r+1; drw(dg,p);
      if(screen(p,sign)){ m=rival_margin(p,sign,&a,-1e-9); if(m>=-1e-9){ h1++; if(h1>=8) break; } } }
    if(h1==0){ ndeep++;
      for(int r=0;r<K2;r++){ drw(dg,p);
        if(screen(p,sign)){ m=rival_margin(p,sign,&a,-1e-9); if(m>=-1e-9){ h2++; if(h2>=4) break; } } } }
    if(h1||h2){ nemit++; printf("%ld %d %d %d %d %d\n",g,nf,h1,r1,h2,K2); }
  }
  fprintf(stderr,"ccensus sign=%d [%ld,%ld) classC %ld deep %ld emitted %ld\n",sign,lo,hi,nC,ndeep,nemit);
  return 0;
}
