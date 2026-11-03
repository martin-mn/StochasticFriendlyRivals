/* run dimprobe's exact per-stratum logic on a single stratum index, verbosely */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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
  long g=atol(argv[1]); int K=atoi(argv[2]); int sign=atoi(argv[3]);
  for(int k=0;k<NWIT;k++) for(int j=0;j<NS;j++) WP[k][j]=(double)((WITG[k]>>j)&1);
  int dg[16]; long x=g; int nf=0;
  for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; if(dg[j]==1)nf++; }
  printf("stratum %ld  free=%d  digits=", g, nf);
  for(int j=0;j<16;j++) putchar(dg[j]==0?'0':(dg[j]==2?'1':'*'));
  putchar('\n');
  rs=(uint64_t)g*6364136223846793005ULL+1442695040888963407ULL; if(!rs) rs=1;
  double p[NS],W[4];
  for(int j=0;j<16;j++) p[j]= dg[j]==0?0.0:(dg[j]==2?1.0:0.02+0.96*urnd());
  freqs(p,p,W);
  char cls = W[0]>1-1e-12?'C':(W[1]+W[2]>1-1e-12?'A':(W[3]<1e-12?'Z':'-'));
  printf("class draw -> W=(%.4f,%.4f,%.4f,%.4f) class=%c\n",W[0],W[1],W[2],W[3],cls);
  int pass=0, hit=0; double best=-1e18;
  for(int r=0;r<K;r++){
    for(int j=0;j<16;j++) p[j]= dg[j]==0?0.0:(dg[j]==2?1.0:0.02+0.96*urnd());
    if(!screen(p,sign)) continue;
    pass++;
    long a; double m=rival_margin(p,sign,&a,-1e-9);
    if(m>best) best=m;
    if(m>=-1e-9){ hit++; if(hit==1) printf("  FIRST HIT at draw %d, margin %.3g\n",r,m); }
  }
  printf("K=%d draws: %d passed the screen, %d confirmed rivals, best margin %.6g\n",K,pass,hit,best);
  return 0;
}
