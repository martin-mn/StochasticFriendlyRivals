/* prodcls.c -- the production classifier of PartnersRivalsM2f/cannon/wfc/wfc.c, re-run on a stored
 * strategy with FULL scans, so that the co-players that decided each bit and their double-precision
 * values can be reported.  It links production's m2f.c unchanged (src/m2f.c, sha256 511bf5b1...).
 *
 * Production (wfc.c classify, rivscreen, rival_margin, nash_scan), with R_ = 1, S_ = cs, T_ = ct of games.dat:
 *   W = freqs(p,p);  below = (UU+VV < 1) with UU = S_, VV = T_-1
 *   eff  = below ? W[0] > 1-1e-12 : W[1]+W[2] > 1-1e-12
 *   sign = (UU-VV < 1) ? +1 : -1
 *   riv  = no binary q with sign*(Wq[2]-Wq[1]) < -1e-9,         Wq = freqs(p,q)
 *          (rivscreen tests 157 of the q with the same expression first; equivalent)
 *   pp   = R_*W[0]+S_*W[1]+T_*W[2]
 *   nash = eff&&riv ? 1 (theorem of Revised2; scanned only for the first 20 per task)
 *                   : no binary q with payq(q,p) > pp+1e-9,   payq = R_*Wq[0]+S_*Wq[1]+T_*Wq[2], Wq = freqs(q,p)
 *   atom = 4 eff + 2 nash + riv
 * The expressions below are written exactly as in wfc.c, so that the compiler (production: GCC 8.5.0,
 * -O3 -march=native on Cannon) makes the same contractions.
 */
#include <math.h>
#define NS 16
void freqs(const double*p,const double*q,double W[4]);

static double R_=1, S_, T_, UU, VV;
static inline double payq(const double*q,const double*p){ double W[4]; freqs(q,p,W);
  return R_*W[0]+S_*W[1]+T_*W[2]; }

/* dout: [0..3] self-play W (double), [4] pp, [5] min margin, [6] max (payq - pp), [7] payq at argmax
 * iout: [0] atom, [1] eff, [2] riv, [3] nash as production sets it, [4] nash by full scan,
 *       [5] argmin margin q, [6] argmax gain q, [7] #q margin < -1e-9, [8] #q with -1e-9 <= margin < 0,
 *       [9] #q payq > pp+1e-9, [10] #q with pp < payq <= pp+1e-9, [11] first q (in 0..65535 order) with
 *       margin < -1e-9 or -1, [12] first q with payq > pp+1e-9 or -1 */
int prod_classify(const double*p,double cs,double ct,double*dout,long*iout){
  R_=1.0; S_=cs; T_=ct; UU=cs; VV=ct-1.0;
  double W[4]; freqs(p,p,W);
  int below=(UU+VV<1.0);
  int eff = below ? (W[0]>1-1e-12) : (W[1]+W[2]>1-1e-12);
  int sign=(UU-VV<1.0)?+1:-1;
  double pp=R_*W[0]+S_*W[1]+T_*W[2];
  double q[NS], Wq[4];
  double mmin=1e18, gmax=-1e18, pqmax=0; long amin=-1, amax=-1, nneg=0, ntol=0, ngain=0, ngtol=0, fneg=-1, fgain=-1;
  for(long g=0; g<65536; g++){
    for(int j=0;j<NS;j++) q[j]=(double)((g>>j)&1);
    freqs(p,q,Wq);
    double d=sign*(Wq[2]-Wq[1]);
    if(d<mmin){ mmin=d; amin=g; }
    if(d < -1e-9){ nneg++; if(fneg<0) fneg=g; } else if(d<0) ntol++;
    double pq=payq(q,p);
    if(pq>pp+1e-9){ ngain++; if(fgain<0) fgain=g; } else if(pq>pp) ngtol++;
    if(pq-pp>gmax){ gmax=pq-pp; amax=g; pqmax=pq; }
  }
  int riv=(nneg==0), nashscan=(ngain==0);
  int nash = (eff&&riv) ? 1 : nashscan;
  for(int k=0;k<4;k++) dout[k]=W[k];
  dout[4]=pp; dout[5]=mmin; dout[6]=gmax; dout[7]=pqmax;
  iout[0]=4*eff+2*nash+riv; iout[1]=eff; iout[2]=riv; iout[3]=nash; iout[4]=nashscan;
  iout[5]=amin; iout[6]=amax; iout[7]=nneg; iout[8]=ntol; iout[9]=ngain; iout[10]=ngtol; iout[11]=fneg; iout[12]=fgain;
  return iout[0];
}

/* the production double values at one co-player g: dout[0] = sign*(W[2]-W[1]) with W = freqs(p,q),
 * dout[1] = payq(q,p), dout[2] = pp, dout[3..6] = freqs(p,q) */
void prod_pair(const double*p,double cs,double ct,long g,double*dout){
  R_=1.0; S_=cs; T_=ct; UU=cs; VV=ct-1.0;
  double W[4], q[NS], Wq[4]; freqs(p,p,W);
  int sign=(UU-VV<1.0)?+1:-1;
  double pp=R_*W[0]+S_*W[1]+T_*W[2];
  for(int j=0;j<NS;j++) q[j]=(double)((g>>j)&1);
  freqs(p,q,Wq);
  dout[0]=sign*(Wq[2]-Wq[1]); dout[1]=payq(q,p); dout[2]=pp;
  for(int k=0;k<4;k++) dout[3+k]=Wq[k];
}
