/* certeng.c -- leading-order eps->0 engine for a pair of memory-two strategies, in two versions:
 *
 *   f_*  : double-precision coefficients (a copy of PartnersRivalsM2f/cannon/m2f.c build + gth),
 *   i_*  : outward-rounded interval coefficients (the scratch ieng.c, made rigorous for interval
 *          inputs and sped up),
 *
 * both extended to return, besides the limit distribution, the INTEGER minimum spanning-tree
 * exponent E and the leading coefficient D* of the tree sum D(eps) = sum_k A_k(eps).
 *
 * Model (identical to m2f.c): state j = 4*x + y, x the last outcome, y the one before; outcome
 * o = 2*(focal defects) + (co-player defects): CC=0 CD=1 DC=2 DD=3, focal action first.  The
 * focal player cooperates at state j with probability eps + (1-2eps) p[j]; the co-player sees the
 * state mirror(j) and cooperates with probability eps + (1-2eps) q[mirror(j)].
 * Leading order of Pr(C) = eps + (1-2eps)x :  (x, eps^0) if x > 0,  (1, eps^1) if x = 0;
 *            of Pr(D) = 1-eps-(1-2eps)x  :  (1-x, eps^0) if x < 1,  (1, eps^1) if x = 1.
 *
 * Exponents: GTH state reduction eliminates states 15..1.  The pivot at step k is
 *   s_k = sum_{j<k} P^{(k)}_{kj}  (P^{(k)} the chain censored on {0..k}),
 * and prod_k s_k = det of (I-P) with row/column 0 deleted = A_0(eps), the sum over spanning trees
 * rooted at 0 (matrix-tree theorem).  So E_0 = sum_k exponent(s_k) and the eps^E_0 coefficient of
 * A_0 is prod_k coefficient(s_k).  GTH back-substitution gives pi_k/pi_0 to leading order
 * (coefficient pc[k], exponent pe[k]), hence E_k = E_0 + pe[k] (the minimum exponent of a tree
 * rooted at k) and E = min_k E_k = E_0 + min_k pe[k].  D* = A_0coef * sum_{pe[k]=min} pc[k].
 * All exponent decisions are integer comparisons that never look at coefficient values, so E,
 * E_k and the support are exact in both versions; only coefficients carry rounding.
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

#define NS 16
#define EINF 1000000000

static inline int swapo(int o){ return ((o&1)<<1)|(o>>1); }
static inline int mirror_state(int j){ return 4*swapo(j>>2)+swapo(j&3); }

/* =============================================================== float engine ============== */
static inline void flead(double x,double*cC,int*eC,double*cD,int*eD){
  if(x>0.0){ *cC=x;   *eC=0; } else { *cC=1.0; *eC=1; }
  if(x<1.0){ *cD=1-x; *eD=0; } else { *cD=1.0; *eD=1; }
}
static void fbuild(const double*p,const double*q,double c[NS][NS],int e[NS][NS]){
  for(int i=0;i<NS;i++) for(int j=0;j<NS;j++){ c[i][j]=0.0; e[i][j]=EINF; }
  for(int j=0;j<NS;j++){
    int jt=mirror_state(j);
    double s[2],t[2]; int es[2],et[2];
    flead(p[j],&s[0],&es[0],&s[1],&es[1]);
    flead(q[jt],&t[0],&et[0],&t[1],&et[1]);
    for(int dF=0;dF<2;dF++) for(int dA=0;dA<2;dA++){
      int nxt=4*(2*dF+dA)+(j>>2);
      double cc=s[dF]*t[dA]; int ee=es[dF]+et[dA];
      if(ee<e[j][nxt]){ e[j][nxt]=ee; c[j][nxt]=cc; }
      else if(ee==e[j][nxt]) c[j][nxt]+=cc;
    }
  }
}
/* w: limit distribution (16); Ek: min tree exponent rooted at k; *E = min_k Ek; *D = D*; *supp: bit k set iff Ek==E */
static void fgth(double c[NS][NS],int e[NS][NS],double*w,int*Ek,int*E,double*D,int*supp){
  int E0=0; double A0=1.0;
  for(int k=NS-1;k>=1;k--){
    int es=EINF; double cs=0;
    for(int j=0;j<k;j++){ int ej=e[k][j];
      if(ej<es){es=ej;cs=c[k][j];} else if(ej==es&&es<EINF){cs+=c[k][j];} }
    E0+=es; A0*=cs;
    for(int i=0;i<k;i++){ if(e[i][k]<EINF){ c[i][k]/=cs; e[i][k]-=es; } }
    for(int i=0;i<k;i++){ if(e[i][k]>=EINF) continue; int eik=e[i][k]; double cik=c[i][k];
      for(int j=0;j<k;j++){ if(e[k][j]>=EINF) continue; int en=eik+e[k][j]; double cn=cik*c[k][j];
        if(en<e[i][j]){e[i][j]=en;c[i][j]=cn;} else if(en==e[i][j]){c[i][j]+=cn;} } }
  }
  double pc[NS]; int pe[NS]; pc[0]=1; pe[0]=0;
  for(int k=1;k<NS;k++){ int ek=EINF; double ck=0;
    for(int i=0;i<k;i++){ if(e[i][k]>=EINF||pe[i]>=EINF) continue; int en=pe[i]+e[i][k]; double cn=pc[i]*c[i][k];
      if(en<ek){ek=en;ck=cn;} else if(en==ek){ck+=cn;} }
    pc[k]=ck; pe[k]=ek; }
  int emin=EINF; double ct=0;
  for(int k=0;k<NS;k++){ if(pe[k]<emin){emin=pe[k];ct=pc[k];} else if(pe[k]==emin){ct+=pc[k];} }
  int s=0;
  for(int k=0;k<NS;k++){
    if(pe[k]==emin){ w[k]=pc[k]/ct; s|=1<<k; } else w[k]=0.0;
    if(Ek) Ek[k]=(pe[k]<EINF)?E0+pe[k]:EINF;
  }
  if(E) *E=E0+emin;
  if(D) *D=A0*ct;
  if(supp) *supp=s;
}
static inline void f_core(const double*p,const double*q,double*w,int*Ek,int*E,double*D,int*supp){
  double c[NS][NS]; int e[NS][NS];
  fbuild(p,q,c,e); fgth(c,e,w,Ek,E,D,supp);
}
static inline void g2p(long g,double*q){ for(int j=0;j<NS;j++) q[j]=(double)((g>>j)&1); }

/* =============================================================== interval engine =========== */
typedef struct { double lo,hi; } I;
/* Outward rounding.  For a correctly rounded op with exact result t>0 and computed r = fl(t) in
 * the normal range, fl(r*(1-2^-51)) <= t <= fl(r*(1+2^-51)) (proof: |r-t| <= u t, u = 2^-53, and
 * (1-4u)(1+u)^2 <= 1, (1+4u)(1-u) >= 1+u).  Outside [2^-1000, 2^1000] fall back to nextafter,
 * which is valid for any single correctly rounded operation. */
#define RLO (1.0-0x1p-51)
#define RHI (1.0+0x1p-51)
static inline double dn(double r){ return (r>0x1p-1000&&r<0x1p1000)? r*RLO : nextafter(r,-INFINITY); }
static inline double up(double r){ return (r>0x1p-1000&&r<0x1p1000)? r*RHI : nextafter(r, INFINITY); }
/* all quantities below are nonnegative intervals */
static inline I iadd(I a,I b){ I r={dn(a.lo+b.lo),up(a.hi+b.hi)}; return r; }
static inline I imul(I a,I b){ I r={dn(a.lo*b.lo),up(a.hi*b.hi)}; return r; }
static inline I idiv(I a,I b){ I r={dn(a.lo/b.hi),up(a.hi/b.lo)}; return r; }
static const I ONE={1.0,1.0};

/* xlo,xhi: an interval inside (0,1) (free), or exactly [0,0] or [1,1] (pinned).  returns -1 if invalid */
static inline int ilead(double xlo,double xhi,I*cC,int*eC,I*cD,int*eD){
  if(xlo==0.0&&xhi==0.0){ *cC=ONE; *eC=1; *cD=ONE; *eD=0; return 0; }
  if(xlo==1.0&&xhi==1.0){ *cC=ONE; *eC=0; *cD=ONE; *eD=1; return 0; }
  if(!(xlo>0.0&&xhi<1.0&&xlo<=xhi)) return -1;
  cC->lo=xlo; cC->hi=xhi; *eC=0;
  /* 1-x: exact if x>=1/2 (Sterbenz); otherwise rounded, enclosed by one ulp each way */
  double a=1.0-xhi, b=1.0-xlo;
  cD->lo=(xhi>=0.5)?a:nextafter(a,-INFINITY); cD->hi=(xlo>=0.5)?b:nextafter(b,INFINITY); *eD=0;
  return 0;
}
/* interval core. wlo/whi: 16-state limit distribution enclosure; W: 4 outcome enclosures (focal
 * view); Ek, E, supp exact integers; D enclosure of D*.  q is a point in [0,1]^16 (usually binary). */
static int i_core(const double*plo,const double*phi,const double*q,I*w,I*W,int*Ek,int*E,I*D,int*supp){
  I c[NS][NS]; int e[NS][NS];
  for(int i=0;i<NS;i++) for(int j=0;j<NS;j++){ c[i][j].lo=c[i][j].hi=0; e[i][j]=EINF; }
  for(int j=0;j<NS;j++){ int jt=mirror_state(j); I s[2],t[2]; int es[2],et[2];
    if(ilead(plo[j],phi[j],&s[0],&es[0],&s[1],&es[1])) return -1;
    if(ilead(q[jt],q[jt],&t[0],&et[0],&t[1],&et[1])) return -2;
    for(int dF=0;dF<2;dF++) for(int dA=0;dA<2;dA++){ int nxt=4*(2*dF+dA)+(j>>2);
      I cc=imul(s[dF],t[dA]); int ee=es[dF]+et[dA];
      if(ee<e[j][nxt]){ e[j][nxt]=ee; c[j][nxt]=cc; } else if(ee==e[j][nxt]) c[j][nxt]=iadd(c[j][nxt],cc); } }
  int E0=0; I A0=ONE;
  for(int k=NS-1;k>=1;k--){ int es=EINF; I cs={0,0};
    for(int j=0;j<k;j++){ int ej=e[k][j]; if(ej<es){es=ej;cs=c[k][j];} else if(ej==es&&es<EINF){cs=iadd(cs,c[k][j]);} }
    E0+=es; A0=imul(A0,cs);
    for(int i=0;i<k;i++) if(e[i][k]<EINF){ c[i][k]=idiv(c[i][k],cs); e[i][k]-=es; }
    for(int i=0;i<k;i++){ if(e[i][k]>=EINF) continue; for(int j=0;j<k;j++){ if(e[k][j]>=EINF) continue;
        int en=e[i][k]+e[k][j]; I cn=imul(c[i][k],c[k][j]);
        if(en<e[i][j]){e[i][j]=en;c[i][j]=cn;} else if(en==e[i][j]) c[i][j]=iadd(c[i][j],cn); } } }
  I pc[NS]; int pe[NS]; pc[0]=ONE; pe[0]=0;
  for(int k=1;k<NS;k++){ int ek=EINF; I ck={0,0};
    for(int i=0;i<k;i++){ if(e[i][k]>=EINF||pe[i]>=EINF) continue; int en=pe[i]+e[i][k]; I cn=imul(pc[i],c[i][k]);
      if(en<ek){ek=en;ck=cn;} else if(en==ek) ck=iadd(ck,cn); }
    pc[k]=ck; pe[k]=ek; }
  int emin=EINF; for(int k=0;k<NS;k++) if(pe[k]<emin) emin=pe[k];
  I ct={0,0}; int first=1;
  for(int k=0;k<NS;k++) if(pe[k]==emin){ if(first){ct=pc[k];first=0;} else ct=iadd(ct,pc[k]); }
  int s=0; int Wset[4]={0,0,0,0};
  for(int o=0;o<4;o++){ W[o].lo=W[o].hi=0; }
  for(int k=0;k<NS;k++){
    if(Ek) Ek[k]=(pe[k]<EINF)?E0+pe[k]:EINF;
    if(pe[k]==emin){ s|=1<<k; I wk=idiv(pc[k],ct); if(wk.hi>1.0) wk.hi=1.0; w[k]=wk; int o=k>>2;
      if(!Wset[o]){ W[o]=wk; Wset[o]=1; } else W[o]=iadd(W[o],wk); }
    else { w[k].lo=w[k].hi=0; }
  }
  for(int o=0;o<4;o++) if(W[o].hi>1.0) W[o].hi=1.0;
  *E=E0+emin; *D=imul(A0,ct); *supp=s;
  return 0;
}
/* sign*(W_DC - W_CD) as an enclosure; *cls: 0 structural zero (no support state with recent outcome
 * CD or DC, so both are exactly 0), 1 certified > 0, -1 certified < 0, 2 enclosure contains 0 */
static inline void imargin(const I*W,int supp,int sign,double*mlo,double*mhi,int*cls){
  int on=0; for(int k=4;k<12;k++) if(supp>>k&1){ on=1; break; }   /* states 4..11: recent CD or DC */
  if(!on){ *mlo=*mhi=0.0; *cls=0; return; }
  double lo,hi;
  if(sign>0){ lo=nextafter(W[2].lo-W[1].hi,-INFINITY); hi=nextafter(W[2].hi-W[1].lo,INFINITY); }
  else      { lo=nextafter(W[1].lo-W[2].hi,-INFINITY); hi=nextafter(W[1].hi-W[2].lo,INFINITY); }
  *mlo=lo; *mhi=hi; *cls=(lo>0)?1:((hi<0)?-1:2);
}

/* =============================================================== exported: single pairs ==== */
void f_eval(const double*p,const double*q,double*w16,int*Ek16,int*E,double*D,int*supp){
  f_core(p,q,w16,Ek16,E,D,supp);
}
/* returns 0 ok, -1 invalid focal interval, -2 invalid co-player */
int i_eval(const double*plo,const double*phi,const double*q,double*wlo16,double*whi16,
           double*Wlo4,double*Whi4,int*Ek16,int*E,double*Dlo,double*Dhi,int*supp){
  I w[NS],W[4],D; int r=i_core(plo,phi,q,w,W,Ek16,E,&D,supp); if(r) return r;
  for(int k=0;k<NS;k++){ wlo16[k]=w[k].lo; whi16[k]=w[k].hi; }
  for(int o=0;o<4;o++){ Wlo4[o]=W[o].lo; Whi4[o]=W[o].hi; }
  *Dlo=D.lo; *Dhi=D.hi; return 0;
}

/* =============================================================== threading helper ========== */
typedef void (*rangefn)(void*ctx,long a,long b);
typedef struct { rangefn fn; void*ctx; long a,b; } TJ;
static void* trun(void*v){ TJ*t=(TJ*)v; t->fn(t->ctx,t->a,t->b); return NULL; }
static void par(rangefn fn,void*ctx,long n,int nth){
  if(nth<=1||n<2*nth){ fn(ctx,0,n); return; }
  if(nth>64) nth=64;
  pthread_t th[64]; TJ tj[64];
  for(int i=0;i<nth;i++){ tj[i].fn=fn; tj[i].ctx=ctx; tj[i].a=n*i/nth; tj[i].b=n*(i+1)/nth; pthread_create(&th[i],NULL,trun,&tj[i]); }
  for(int i=0;i<nth;i++) pthread_join(th[i],NULL);
}
static inline long qget(const long*ql,long i){ return ql?ql[i]:i; }

/* =============================================================== float scan ================ */
typedef struct { const double*p; const long*ql; double*W; int*E; int*supp; double*D; } FS;
static void fs_run(void*v,long a,long b){ FS*c=(FS*)v; double q[NS],w[NS];
  for(long i=a;i<b;i++){ g2p(qget(c->ql,i),q); int E,s; double D; f_core(c->p,q,w,NULL,&E,&D,&s);
    double*o=c->W+4*i; o[0]=o[1]=o[2]=o[3]=0; for(int k=0;k<NS;k++) o[k>>2]+=w[k];
    c->E[i]=E; c->supp[i]=s; c->D[i]=D; } }
/* all binary co-players listed in ql (NULL: all 65536, nq must be 65536) */
void f_scan(const double*p,long nq,const long*ql,double*W4,int*E,int*supp,double*D,int nth){
  FS c={p,ql,W4,E,supp,D}; par(fs_run,&c,nq,nth);
}

/* =============================================================== interval scan ============= */
typedef struct { const double*plo,*phi; int sign; const long*ql; int*cls; double*mlo,*mhi; int*E; int*supp; int err; } IS;
static void is_run(void*v,long a,long b){ IS*c=(IS*)v; double q[NS]; I w[NS],W[4],D; int E,s;
  for(long i=a;i<b;i++){ g2p(qget(c->ql,i),q);
    if(i_core(c->plo,c->phi,q,w,W,NULL,&E,&D,&s)){ c->err=1; c->cls[i]=-9; continue; }
    imargin(W,s,c->sign,&c->mlo[i],&c->mhi[i],&c->cls[i]); c->E[i]=E; if(c->supp) c->supp[i]=s; } }
int i_scan(const double*plo,const double*phi,int sign,long nq,const long*ql,int*cls,double*mlo,double*mhi,int*E,int*supp,int nth){
  IS c={plo,phi,sign,ql,cls,mlo,mhi,E,supp,0}; par(is_run,&c,nq,nth); return c.err?-1:0;
}

/* =============================================================== corner scan (L3) ========== */
/* pat[16]: 0, 1 or 2 (=free).  For each co-player q in the list:
 *   Es[i]   = E of the open stratum against q (computed at the interior point x_free = 1/2; the
 *             exponents are pattern-determined, so any interior point gives the same integer);
 *   cnt[7i+0..6] over the 2^d corners c of the closed stratum cube:
 *       0 npos : E_c == E and sign*(W_DC-W_CD)(c,q) certified > 0
 *       1 nneg : E_c == E and certified < 0
 *       2 nzs  : E_c == E and structural zero (neither CD nor DC recent in the support)
 *       3 nzi  : E_c == E, enclosure contains 0 and D*_hi * max|enclosure| < 1: exactly 0 by
 *                integrality.  For a BINARY pair every edge's leading coefficient is 1, so A*_k is
 *                the number of minimum-order spanning trees rooted at k (an integer), D* their total,
 *                and D* * (W_DC - W_CD) is an integer; |it| < 1 forces it to be 0.
 *       4 nhigh: E_c >  E   (so N_q(c) = 0 = D*(c))
 *       5 nlow : E_c <  E   (impossible by the lemma; counted as a sanity check)
 *       6 nunr : E_c == E, enclosure contains 0 and integrality does not decide (never expected)
 *   mmin[i], mmax[i]: min lower / max upper enclosure end over corners with E_c == E. */
typedef struct { const signed char*pat; int sign; const long*ql; int*Es; int*cnt; double*mmin,*mmax; int d; int fr[NS]; } CS;
static void cs_run(void*v,long a,long b){ CS*c=(CS*)v; double q[NS],x[NS],w[NS]; I iw[NS],W[4],D;
  for(long i=a;i<b;i++){ g2p(qget(c->ql,i),q);
    for(int j=0;j<NS;j++) x[j]=(c->pat[j]==2)?0.5:(double)c->pat[j];
    int E,s; double Df; f_core(x,q,w,NULL,&E,&Df,&s); c->Es[i]=E;
    int*ct=c->cnt+7*i; for(int t=0;t<7;t++) ct[t]=0;
    double mn=INFINITY,mx=-INFINITY;
    for(long m=0;m<(1L<<c->d);m++){
      for(int t=0;t<c->d;t++) x[c->fr[t]]=(double)((m>>t)&1);
      int Ec,sc; i_core(x,x,q,iw,W,NULL,&Ec,&D,&sc);
      if(Ec>E){ ct[4]++; continue; } if(Ec<E){ ct[5]++; continue; }
      double lo,hi; int cl; imargin(W,sc,c->sign,&lo,&hi,&cl);
      if(lo<mn) mn=lo; if(hi>mx) mx=hi;
      if(cl==0) ct[2]++; else if(cl==1) ct[0]++; else if(cl==-1) ct[1]++;
      else if(up(D.hi*fmax(fabs(lo),fabs(hi)))<1.0) ct[3]++; else ct[6]++;
    }
    c->mmin[i]=mn; c->mmax[i]=mx; } }
void corner_scan(const signed char*pat,int sign,long nq,const long*ql,int*Es,int*cnt,double*mmin,double*mmax,int nth){
  CS c; c.pat=pat; c.sign=sign; c.ql=ql; c.Es=Es; c.cnt=cnt; c.mmin=mmin; c.mmax=mmax; c.d=0;
  for(int j=0;j<NS;j++) if(pat[j]==2) c.fr[c.d++]=j;
  par(cs_run,&c,nq,nth);
}
/* per-corner detail for ONE co-player: for corner m (bit t of m = value of the t-th free coordinate):
 * Ec[m], cls[m] (as imargin, 3 = zero by integrality, or 9 when E_c != E), mlo[m], mhi[m], Dhi[m].
 * Returns E of the open stratum. */
int corner_detail(const signed char*pat,int sign,long g,int*Ec,int*cls,double*mlo,double*mhi,double*Dhi){
  double q[NS],x[NS],w[NS]; I iw[NS],W[4],D; int fr[NS],d=0;
  for(int j=0;j<NS;j++){ if(pat[j]==2) fr[d++]=j; x[j]=(pat[j]==2)?0.5:(double)pat[j]; }
  g2p(g,q); int E,s; double Df; f_core(x,q,w,NULL,&E,&Df,&s);
  for(long m=0;m<(1L<<d);m++){
    for(int t=0;t<d;t++) x[fr[t]]=(double)((m>>t)&1);
    int e,sc; i_core(x,x,q,iw,W,NULL,&e,&D,&sc); Ec[m]=e; Dhi[m]=D.hi;
    if(e!=E){ cls[m]=9; mlo[m]=mhi[m]=0; continue; }
    imargin(W,sc,sign,&mlo[m],&mhi[m],&cls[m]);
    if(cls[m]==2 && up(D.hi*fmax(fabs(mlo[m]),fabs(mhi[m])))<1.0) cls[m]=3;
  }
  return E;
}

/* =============================================================== vertex scan (L2) ========== */
/* nv vertices given as enclosures vlo/vhi (nv x 16, each an interval inside the open stratum or
 * a pinned 0/1).  For each co-player: cnt[4i+0..3] = #vertices certified >0, certified <0,
 * structural zero, enclosure containing 0; mmin[i] = min lower end.  Up to cap unresolved
 * (list index, vertex index) pairs are written to un[2*..]; *nun receives the total number. */
typedef struct { const double*vlo,*vhi; long nv; int sign; const long*ql; int*cnt; double*mmin; long*un; long cap; long nun; pthread_mutex_t mu; int err; } VS;
static void vs_run(void*v,long a,long b){ VS*c=(VS*)v; double q[NS]; I w[NS],W[4],D; int E,s;
  for(long i=a;i<b;i++){ g2p(qget(c->ql,i),q); int*ct=c->cnt+4*i; ct[0]=ct[1]=ct[2]=ct[3]=0; double mn=INFINITY;
    for(long k=0;k<c->nv;k++){
      if(i_core(c->vlo+16*k,c->vhi+16*k,q,w,W,NULL,&E,&D,&s)){ c->err=1; continue; }
      double lo,hi; int cl; imargin(W,s,c->sign,&lo,&hi,&cl); if(lo<mn) mn=lo;
      if(cl==1) ct[0]++; else if(cl==-1) ct[1]++; else if(cl==0) ct[2]++;
      else { ct[3]++; pthread_mutex_lock(&c->mu); if(c->nun<c->cap){ c->un[2*c->nun]=i; c->un[2*c->nun+1]=k; } c->nun++; pthread_mutex_unlock(&c->mu); }
    }
    c->mmin[i]=mn; } }
int vertex_scan(const double*vlo,const double*vhi,long nv,int sign,long nq,const long*ql,int*cnt,double*mmin,long*un,long cap,long*nun,int nth){
  VS c; c.vlo=vlo; c.vhi=vhi; c.nv=nv; c.sign=sign; c.ql=ql; c.cnt=cnt; c.mmin=mmin; c.un=un; c.cap=cap; c.nun=0; c.err=0;
  pthread_mutex_init(&c.mu,NULL); par(vs_run,&c,nq,nth); pthread_mutex_destroy(&c.mu);
  *nun=c.nun; return c.err?-1:0;
}

/* =============================================================== self-play ================= */
/* self-play support (exact integer) and float frequencies */
void f_self(const double*p,double*W4,int*supp,int*E){
  double w[NS]; double D; f_core(p,p,w,NULL,E,&D,supp);
  W4[0]=W4[1]=W4[2]=W4[3]=0; for(int k=0;k<NS;k++) W4[k>>2]+=w[k];
}

/* =============================================================== interval scan, full W ======= */
/* enclosures of (w_CC, w_CD, w_DC, w_DD) for every listed co-player, plus E and the support */
typedef struct { const double*plo,*phi; const long*ql; double*Wlo,*Whi; int*E; int*supp; int err; } IW;
static void iw_run(void*v,long a,long b){ IW*c=(IW*)v; double q[NS]; I w[NS],W[4],D; int E,s;
  for(long i=a;i<b;i++){ g2p(qget(c->ql,i),q);
    if(i_core(c->plo,c->phi,q,w,W,NULL,&E,&D,&s)){ c->err=1; continue; }
    for(int o=0;o<4;o++){ c->Wlo[4*i+o]=W[o].lo; c->Whi[4*i+o]=W[o].hi; } c->E[i]=E; c->supp[i]=s; } }
int i_scanW(const double*plo,const double*phi,long nq,const long*ql,double*Wlo,double*Whi,int*E,int*supp,int nth){
  IW c={plo,phi,ql,Wlo,Whi,E,supp,0}; par(iw_run,&c,nq,nth); return c.err?-1:0;
}
