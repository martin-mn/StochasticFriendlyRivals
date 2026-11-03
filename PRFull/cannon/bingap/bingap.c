/* bingap.c -- exhaustive census of the binary-pair competitive margins m(p,q) = w_DC - w_CD over all
 * 65536 x 65536 binary memory-two pairs, in the eps -> 0 limit, rigorously:
 *   every edge of a binary pair has leading coefficient 1, so A*_k = #minimum-order spanning trees
 *   rooted at k and D* = sum_k A*_k are integers, and N = D* m = sum_DC A*_k - sum_CD A*_k is an
 *   integer.  The interval engine (certeng.c) encloses D* and m; each integer is taken only when its
 *   enclosure isolates it (otherwise the pair is counted in 'flag' and printed).
 * Symmetry: m(q,p) = -m(p,q) and D*(q,p) = D*(p,q) (the (q,p) chain is the mirror of (p,q)), so
 * only q >= p is computed; --checksym also computes q < p for the task's p and verifies it.
 * Task t of T handles every p with p % T == t (balanced load).
 * Output: one line of totals per task (see README.md for the fields). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
int i_eval(const double*plo,const double*phi,const double*q,double*wlo16,double*whi16,
           double*Wlo4,double*Whi4,int*Ek16,int*E,double*Dlo,double*Dhi,int*supp);
static void g2p(long g,double*q){ for(int j=0;j<16;j++) q[j]=(double)((g>>j)&1); }
/* returns 0 ok (N, D set), 1 not isolated */
static int pair(long p,long q,long long*N,long long*D){
  double pp[16],qq[16],wlo[16],whi[16],Wlo[4],Whi[4],Dlo,Dhi; int Ek[16],E,s;
  g2p(p,pp); g2p(q,qq);
  if(i_eval(pp,pp,qq,wlo,whi,Wlo,Whi,Ek,&E,&Dlo,&Dhi,&s)) return 2;
  long long d=(long long)ceil(Dlo); if(!((double)d<=Dhi && Dhi<(double)(d+1))) return 1;
  double lo=nextafter(Wlo[2]-Whi[1],-INFINITY), hi=nextafter(Whi[2]-Wlo[1],INFINITY);
  double nlo=nextafter(lo*(double)d,-INFINITY), nhi=nextafter(hi*(double)d,INFINITY);
  long long n=(long long)ceil(nlo); if(!((double)n<=nhi && nhi<(double)(n+1))) return 1;
  *N=n; *D=d; return 0;
}
int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage: bingap TASK NTASKS [--checksym] [--pmax K]\n"); return 2; }
  int t=atoi(argv[1]), T=atoi(argv[2]); int checksym=0; long pmax=65536;
  for(int a=3;a<argc;a++){ if(!strcmp(argv[a],"--checksym")) checksym=1; else if(!strcmp(argv[a],"--pmax")&&a+1<argc) pmax=atol(argv[++a]); }
  long long npairs=0,nzero=0,nflag=0,nsymbad=0,nbelow37=0,nbelow1e9=0;
  long long bestN=0,bestD=1; long bp=-1,bq=-1; long long maxD=0; long mp=-1,mq=-1;
  struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
  for(long p=t;p<pmax;p+=T){
    for(long q=p;q<65536;q++){
      long long N,D; int r=pair(p,q,&N,&D); npairs++;
      if(r){ nflag++; if(nflag<=20) fprintf(stderr,"FLAG p=%ld q=%ld code %d\n",p,q,r); continue; }
      if(D>maxD){ maxD=D; mp=p; mq=q; }
      if(N==0){ nzero++; continue; }
      long long a=llabs(N);
      if((__int128)a*37 < (__int128)D) nbelow37++;
      if((double)a/(double)D < 1e-9) nbelow1e9++;
      if(bp<0 || (__int128)a*bestD < (__int128)bestN*D){ bestN=a; bestD=D; bp=p; bq=q; }
    }
    if(checksym){
      for(long q=0;q<p;q+=97){ long long N1,D1,N2,D2;
        if(pair(p,q,&N1,&D1)||pair(q,p,&N2,&D2)){ nflag++; continue; }
        if(!(N1==-N2 && D1==D2)) nsymbad++; }
    }
  }
  clock_gettime(CLOCK_MONOTONIC,&t1); double dt=(t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
  printf("task %d of %d: pairs %lld zero %lld flagged %lld symmetry_violations %lld below_1/37 %lld below_1e-9 %lld "
         "min_nonzero %lld/%lld at p=%ld q=%ld max_Dstar %lld at p=%ld q=%ld seconds %.1f us_per_pair %.3f\n",
         t,T,npairs,nzero,nflag,nsymbad,nbelow37,nbelow1e9,bestN,bestD,bp,bq,maxD,mp,mq,dt,dt/npairs*1e6);
  return (nflag||nsymbad)?1:0;
}
