#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NS 16
void freqs(const double*p,const double*q,double W[4]);
void freqs16(const double*p,const double*q,double w[NS]);
double rival_margin(const double*p,int sign,long*arg,double earlyquit);
static void g2p(unsigned g,double*p){ for(int j=0;j<NS;j++) p[j]=(double)((g>>j)&1); }
int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: m2f self|rival|pair ...\n");return 1;}
  double p[NS],q[NS],W[4];
  if(!strcmp(argv[1],"self")){            /* self stripe nstripes : binary self-play limits */
    int st=atoi(argv[2]),ns=atoi(argv[3]);
    for(unsigned s=st;s<65536;s+=ns){ g2p(s,p); freqs(p,p,W);
      printf("%u %.17g %.17g %.17g %.17g\n",s,W[0],W[1],W[2],W[3]); }
  } else if(!strcmp(argv[1],"rival")){    /* rival stripe nstripes sign : binary rival test */
    int st=atoi(argv[2]),ns=atoi(argv[3]),sign=atoi(argv[4]);
    for(unsigned s=st;s<65536;s+=ns){ g2p(s,p); long a; double m=rival_margin(p,sign,&a,-1e-9);
      printf("%u %.17g %ld\n",s,m,a); }
  } else if(!strcmp(argv[1],"pair")){     /* pair g1 g2 */
    g2p((unsigned)atol(argv[2]),p); g2p((unsigned)atol(argv[3]),q); freqs(p,q,W);
    printf("%s %s %.17g %.17g %.17g %.17g\n",argv[2],argv[3],W[0],W[1],W[2],W[3]);
  }
  return 0;
}
