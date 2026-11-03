/* strata.c -- enumerate all 3^16 = 43 046 721 support strata of the memory-two cube and
 * classify the eps->0 SELF-PLAY limit of each.  Digit d_j in {0,1,2} means p_j = 0,
 * p_j in (0,1), p_j = 1.  The three efficiency conditions
 *      W_CC = 1            (efficient iff u+v < 1)
 *      W_CD = W_DC = 1/2   (efficient iff u+v > 1)
 *      W_DD = 0            (efficient iff u+v = 1)
 * depend on p only through this pattern -- they say which states carry limit mass, and the
 * support of the eps=0 chain plus the eps-orders of its escapes determine that.  Two
 * independent random interior draws per stratum guard against a non-generic coincidence.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define NS 16
void freqs(const double*p,const double*q,double W[4]);

static uint64_t rs;
static inline double urnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }

int main(int argc,char**argv){
  long stripe=atoi(argv[1]), nstr=atoi(argv[2]);
  const char* outpath = argc>3?argv[3]:NULL;
  FILE*f = outpath?fopen(outpath,"w"):stdout;
  long TOT=1; for(int i=0;i<16;i++) TOT*=3;
  long nCC=0,nALT=0,nDD0=0,ndis=0;
  double p[NS],W1[4],W2[4];
  int dg[16];
  for(long g=stripe;g<TOT;g+=nstr){
    long x=g; for(int j=0;j<16;j++){ dg[j]=x%3; x/=3; }
    rs = (uint64_t)g*6364136223846793005ULL + 1442695040888963407ULL; if(!rs) rs=88172645463325252ULL;
    for(int r=0;r<2;r++){
      for(int j=0;j<16;j++) p[j] = dg[j]==0?0.0 : (dg[j]==2?1.0 : 0.12+0.76*urnd());
      freqs(p,p,r?W2:W1);
    }
    int cc1 = (W1[0]>1-1e-12) , cc2 = (W2[0]>1-1e-12);
    int al1 = (W1[1]+W1[2]>1-1e-12), al2 = (W2[1]+W2[2]>1-1e-12);
    int d01 = (W1[3]<1e-12), d02 = (W2[3]<1e-12);
    if(cc1!=cc2||al1!=al2||d01!=d02){ ndis++; fprintf(stderr,"DISAGREE g=%ld\n",g); }
    if(cc1){nCC++; fprintf(f,"C %ld\n",g);}
    else if(al1){nALT++; fprintf(f,"A %ld\n",g);}
    else if(d01){nDD0++; fprintf(f,"Z %ld\n",g);}
  }
  fprintf(stderr,"stripe %ld/%ld: W_CC=1 %ld, alternating %ld, other W_DD=0 %ld, draw-disagreements %ld\n",
          stripe,nstr,nCC,nALT,nDD0,ndis);
  if(outpath) fclose(f);
  return 0;
}
