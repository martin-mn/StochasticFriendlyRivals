/* wit2.c -- harvest arg-min co-players over random CONTINUUM strategies, to grow the witness
 * set from 19 (which covers the 65536 binary non-rivals) to one that also covers the continuum. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define NS 16
void freqs(const double*p,const double*q,double W[4]);
double rival_margin(const double*p,int sign,long*arg,double earlyquit);
static uint64_t rs;
static inline double urnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }
int main(int argc,char**argv){
  long n=atol(argv[1]); long seed=atol(argv[2]); rs=(uint64_t)seed*88172645463325252ULL+12345; if(!rs)rs=7;
  double p[NS];
  for(long i=0;i<n;i++){
    /* random support pattern, then random interior values */
    for(int j=0;j<16;j++){ double u=urnd(); p[j] = u<0.34?0.0 : (u<0.68?1.0 : 0.02+0.96*urnd()); }
    long aP,aM; double mP=rival_margin(p,+1,&aP,-1e18); double mM=rival_margin(p,-1,&aM,-1e18);
    printf("%ld %.17g %ld %.17g\n",aP,mP,aM,mM);
  }
  return 0;
}
