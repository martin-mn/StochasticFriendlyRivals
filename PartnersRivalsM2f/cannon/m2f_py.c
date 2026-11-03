/* thin batch driver used by check/m2f.py through ctypes */
#include <stdlib.h>
#define NS 16
void freqs(const double*p,const double*q,double W[4]);
void freqs16(const double*p,const double*q,double w[NS]);
double rival_margin(const double*p,int sign,long*arg,double earlyquit);
static void g2p(unsigned g,double*p){ for(int j=0;j<NS;j++) p[j]=(double)((g>>j)&1); }

void c_freqs(const double*p,const double*q,double*W){ freqs(p,q,W); }
void c_freqs16(const double*p,const double*q,double*w){ freqs16(p,q,w); }
void c_self(const double*p,double*W){ freqs(p,p,W); }
double c_rival(const double*p,int sign,long*arg,double eq){ return rival_margin(p,sign,arg,eq); }
/* batch: n strategies of 16 doubles each (row-major). out: 4 self freqs + margP + margM + argP + argM  */
void c_batch(const double*P,long n,double*out,double eq){
  for(long i=0;i<n;i++){
    const double*p=P+16*i; double W[4]; freqs(p,p,W);
    long aP,aM;
    double mP=rival_margin(p,+1,&aP,eq);
    double mM=rival_margin(p,-1,&aM,eq);
    double*o=out+8*i;
    o[0]=W[0];o[1]=W[1];o[2]=W[2];o[3]=W[3];o[4]=mP;o[5]=mM;o[6]=(double)aP;o[7]=(double)aM;
  }
}
/* rival margin against a supplied witness list only (fast screen) */
void c_batch_wit(const double*P,long n,const long*wit,long nw,double*out){
  for(long i=0;i<n;i++){
    const double*p=P+16*i; double q[16],W[4];
    double bP=1e18,bM=1e18;
    for(long k=0;k<nw;k++){ g2p((unsigned)wit[k],q); freqs(p,q,W);
      double dP=W[2]-W[1]; double dM=-dP;
      if(dP<bP)bP=dP; if(dM<bM)bM=dM; }
    out[2*i]=bP; out[2*i+1]=bM;
  }
}
