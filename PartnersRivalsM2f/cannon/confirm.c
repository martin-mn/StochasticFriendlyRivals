/* confirm.c -- full 65536-co-player rival confirmation for a list of support strata.
 * usage: confirm <stratumlist> <sign> <K> <stripe> <nstripes>
 * For each stratum in the list, draw K random interior points and report, for each, the exact
 * rival margin = min over ALL 65536 binary memory-two co-players of sign*(W_DC-W_CD), together
 * with the self-play class.  A margin >= -1e-9 is a rival (the box-vertex lemma makes the
 * binary minimisation exact over the whole memory-two continuum of co-players). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define NS 16
void freqs(const double*p,const double*q,double W[4]);
double rival_margin(const double*p,int sign,long*arg,double earlyquit);
static uint64_t rs;
static inline double urnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }
int main(int argc,char**argv){
  FILE*fl=fopen(argv[1],"r"); int sign=atoi(argv[2]); int K=atoi(argv[3]);
  long stripe=atol(argv[4]), nstr=atol(argv[5]);
  long g; long idx=0, nriv=0, nstrat=0; int dg[16]; double p[NS],W[4];
  while(fscanf(fl,"%ld",&g)==1){
    if((idx++)%nstr!=stripe) continue;
    long x=g; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; }
    rs=(uint64_t)g*2862933555777941757ULL+3037000493ULL; if(!rs)rs=88172645463325252ULL;
    int hit=0;
    for(int r=0;r<K;r++){
      for(int j=0;j<16;j++) p[j]= dg[j]==0?0.0 : (dg[j]==2?1.0 : 0.07+0.86*urnd());
      long a; double m=rival_margin(p,sign,&a,-1e-9);
      if(m>=-1e-9){ hit++; nriv++;
        freqs(p,p,W);
        printf("%ld %d %.17g %.17g %.17g %.17g",g,r,W[0],W[1],W[2],W[3]);
        for(int j=0;j<16;j++) printf(" %.17g",p[j]);
        printf("\n"); }
    }
    if(hit) nstrat++;
  }
  fprintf(stderr,"stripe %ld: confirmed-rival strata %ld, points %ld\n",stripe,nstrat,nriv);
  return 0;
}
