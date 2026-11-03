/* dimprobe.c -- is the wedge-W friendly-rival set of dimension > 6?
 * Sweeps every support stratum that is (a) efficient below the switch line and (b) has at least
 * DMIN free coordinates, drawing K random interior points and confirming rivalry exactly against
 * all 65536 binary co-players.  Reports the largest number of free coordinates achieved. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "witnesses.h"
#define NS 16
void freqs(const double*p,const double*q,double W[4]);
double rival_margin(const double*p,int sign,long*arg,double earlyquit);
static double WP[NWIT][NS];
static uint64_t rs;
static inline double urnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }
static int screen(const double*p,int sign){ double W[4];
  for(int k=0;k<NWIT;k++){ freqs(p,WP[k],W); if(sign*(W[2]-W[1])<-1e-9) return 0; } return 1; }
int main(int argc,char**argv){
  int DMIN=atoi(argv[1]); int K=atoi(argv[2]); int sign=atoi(argv[3]);
  char want=argv[4][0];                 /* 'C' efficient below, 'A' alternating, 'Z' line only */
  long stripe=atol(argv[5]), nstr=atol(argv[6]);
  for(int k=0;k<NWIT;k++) for(int j=0;j<NS;j++) WP[k][j]=(double)((WITG[k]>>j)&1);
  long TOT=1; for(int i=0;i<16;i++) TOT*=3;
  int dg[16]; double p[NS],W[4]; int bestd=-1; long nhit=0;
  for(long g=stripe;g<TOT;g+=nstr){
    long x=g; int nf=0; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; if(dg[j]==1)nf++; }
    if(nf<DMIN) continue;
    rs=(uint64_t)g*6364136223846793005ULL+1442695040888963407ULL; if(!rs)rs=1;
    for(int j=0;j<16;j++) p[j]= dg[j]==0?0.0:(dg[j]==2?1.0:0.02+0.96*urnd());
    freqs(p,p,W);
    char cls = W[0]>1-1e-12?'C' : (W[1]+W[2]>1-1e-12?'A' : (W[3]<1e-12?'Z':'-'));
    if(cls!=want) continue;
    int hit=0;
    for(int r=0;r<K && !hit;r++){
      for(int j=0;j<16;j++) p[j]= dg[j]==0?0.0:(dg[j]==2?1.0:0.02+0.96*urnd());
      if(!screen(p,sign)) continue;
      long a; if(rival_margin(p,sign,&a,-1e-9)>=-1e-9) hit=1;
    }
    if(hit){ nhit++; if(nf>bestd)bestd=nf;
      printf("%ld %d %c\n",g,nf,cls); }
  }
  fprintf(stderr,"stripe %ld: DMIN=%d want=%c sign=%d -> %ld strata, max free coords %d\n",
          stripe,DMIN,want,sign,nhit,bestd);
  return 0;
}
