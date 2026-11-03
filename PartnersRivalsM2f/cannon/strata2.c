/* strata2.c -- for every one of the 3^16 support strata of the memory-two cube:
 *   (1) classify the eps->0 self-play limit  (C = W_CC 1, A = perfect alternation, Z = W_DD 0)
 *   (2) screen K random interior points of the stratum for rivalry against a 19-co-player
 *       witness set (early exit).  A point failing any witness is definitively NOT a rival,
 *       so the screen has no false negatives; survivors are re-tested against all 65536.
 * Output: one line per stratum that has at least one surviving point, with the self-play
 * class, the number of survivors for T>S and for T<S, and a surviving point index.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define NS 16
void freqs(const double*p,const double*q,double W[4]);

static const long WIT[19]={4480,20817,2827,35980,53456,19532,1024,4512,8205,24576,
                           6208,1088,3072,6144,23387,9766,6464,4256,12288};
static double WITP[19][16];
static uint64_t rs;
static inline double urnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }

/* min over the witness set of sign*(W_DC-W_CD); early exit as soon as it is negative */
static double screen(const double*p,int sign){
  double best=1e18, W[4];
  for(int k=0;k<19;k++){ freqs(p,WITP[k],W); double d=sign*(W[2]-W[1]);
    if(d<best){best=d; if(best<-1e-9) return best;} }
  return best;
}

int main(int argc,char**argv){
  long stripe=atol(argv[1]), nstr=atol(argv[2]); int K=argc>3?atoi(argv[3]):4;
  const char*outpath=argc>4?argv[4]:NULL;
  FILE*f=outpath?fopen(outpath,"w"):stdout;
  for(int k=0;k<19;k++) for(int j=0;j<16;j++) WITP[k][j]=(double)((WIT[k]>>j)&1);
  long TOT=1; for(int i=0;i<16;i++) TOT*=3;
  long nC=0,nA=0,nZ=0, hitP=0,hitM=0, effhitP=0,effhitM=0;
  double p[NS],W[4]; int dg[16];
  for(long g=stripe;g<TOT;g+=nstr){
    long x=g; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; }
    rs=(uint64_t)g*6364136223846793005ULL+1442695040888963407ULL; if(!rs) rs=88172645463325252ULL;
    int sP=0,sM=0; char cls=0;
    for(int r=0;r<K;r++){
      for(int j=0;j<16;j++) p[j]= dg[j]==0?0.0 : (dg[j]==2?1.0 : 0.07+0.86*urnd());
      if(r==0){ freqs(p,p,W);
        if(W[0]>1-1e-12) cls='C'; else if(W[1]+W[2]>1-1e-12) cls='A'; else if(W[3]<1e-12) cls='Z'; else cls='-';
        if(cls=='C')nC++; else if(cls=='A')nA++; else if(cls=='Z')nZ++; }
      if(screen(p,+1)>=-1e-9) sP++;
      if(screen(p,-1)>=-1e-9) sM++;
    }
    if(sP||sM){
      if(sP)hitP++; if(sM)hitM++;
      if(sP&&(cls=='C'||cls=='A'||cls=='Z'))effhitP++;
      if(sM&&(cls=='C'||cls=='A'||cls=='Z'))effhitM++;
      fprintf(f,"%ld %c %d %d\n",g,cls,sP,sM);
    }
  }
  fprintf(stderr,"stripe %ld: C %ld A %ld Z %ld | rivalP-strata %ld rivalM-strata %ld | eff&rivP %ld eff&rivM %ld\n",
          stripe,nC,nA,nZ,hitP,hitM,effhitP,effhitM);
  if(outpath)fclose(f);
  return 0;
}
