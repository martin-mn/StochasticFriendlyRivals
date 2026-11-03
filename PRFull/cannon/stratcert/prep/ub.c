/* ub.c -- upper bounds on the dimension of the memory-two friendly-rival set, by excluding whole
 * open strata with one binary co-player each (structural lemma, L3 + L4 of check/cert/certify.py).
 *
 * THE CRITERION (exact).  Fix a binary co-player q and a stratum sigma (each coordinate 0, 1 or
 * free).  N_q = D* . margin is multilinear on the closed cube of sigma, so for x in the OPEN cube
 *     N_q(x) = sum_c N_q(c) prod_i (x_i if c_i = 1 else 1 - x_i),
 * a combination of the 2^d corner values with weights that are all strictly positive.  By L3,
 * N_q(c) = D*_c . m(c,q) with D*_c >= 1 when E_c(q) = E_sigma(q), and N_q(c) = 0 when
 * E_c(q) > E_sigma(q); and E_sigma(q) = min_c E_c(q) (a minimum-order tree of the open stratum is a
 * tree of the corner that copies its actions, and pinning never lowers an order).  Hence
 *     q beats every point of the open stratum (margin < 0 throughout)
 *       <=> every corner with E_c = min E has m(c,q) <= 0, and at least one has m(c,q) < 0.
 * (<= is the convex-combination identity; => by continuity and N_q not identically 0.)
 * Here m = w_DC - w_CD for sign +1 (T>S) and w_CD - w_DC for sign -1 (T<S).
 *
 * THE COMPUTATION.  For one q: the row R[c] = (E_c << 2) | pos<<1 | neg over all 65536 binary
 * focal c, from the outward-rounded interval engine i_core (certeng.c, verbatim copy of
 * check/cert/certeng.c): pos/neg = margin (sign +1) certified > 0 / < 0; a structural zero or an
 * integrality zero (D*_hi |enclosure| < 1, D* . m is an integer for a binary pair) sets neither;
 * an undecided enclosure sets BOTH (blocks exclusion in both signs; counted, never expected).
 * Then a ternary DP over all 3^16 strata (index sum_j dig_j 3^j, dig 0 = pinned 0, 1 = pinned 1,
 * 2 = free):  V[sigma] = combine(V[sigma with j->0], V[sigma with j->1]) for the lowest free j,
 * combine(a,b) = the one with smaller E, or a|b on a tie -- so V[sigma] = (min_c E_c)<<2 | flags
 * ORed over the minimum-E corners.  Exclusion for sign +1: flags == neg only (1); for sign -1:
 * flags == pos only (2).
 *
 * Usage:
 *   ub classes CLS.bin [nth]            efficiency class of every open stratum (uint8, 3^16):
 *                                       0 none, 1 C (w_CC=1), 2 A (w_CD+w_DC=1), 3 Z (w_DD=0 only),
 *                                       from the exact self-play support at the centre (L1).
 *   ub cover CLS.bin QLIST OUT [nth]    for every efficient stratum (class 1..3), the least q of QLIST
 *                                       that excludes it for sign +1 and for sign -1.  Writes
 *                                       OUT.idx (uint32 stratum indices, ascending), OUT.qp, OUT.qm
 *                                       (int32 q or -1) and OUT.log (per-q counts).  If OUT.qp/.qm
 *                                       exist with the same OUT.idx they are merged (min).
 *   ub verify TRIPLES.bin OUT [nth]   independent re-check of (stratum, sign, q) triples (int64 x 3 each):
 *                                       direct corner enumeration + interval centre margin; see ver_run.
 *   ub sets QLIST TRIPLES.bin OUT [nth]  bit matrix: which q excludes which (stratum, sign) triple.
 *   ub rowcmp Q FILE                    compare a 65536-byte row file with the interval row of q.
 *   ub row Q [OUT]                      print the row statistics of one q; OUT: write the 65536 row bytes.
 *   ub check PATTERN SIGN Q             direct 2^d corner enumeration for one stratum (debugging:
 *                                       an independent route to the same criterion, no DP).
 */
#include "certeng.c"
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

#define N3 43046721L
static long P3[17];
static int8_t LF[6561];      /* lowest free digit of an 8-digit block, -1 if none */
static uint16_t BN[6561];    /* binary value of an 8-digit block without free digits */

static double now(void){ struct timeval t; gettimeofday(&t,NULL); return t.tv_sec+1e-6*t.tv_usec; }

static void init_tables(void){
  P3[0]=1; for(int j=1;j<=16;j++) P3[j]=3*P3[j-1];
  for(int b=0;b<6561;b++){ int x=b; LF[b]=-1; BN[b]=0;
    for(int j=0;j<8;j++){ int d=x%3; x/=3; if(d==2){ if(LF[b]<0) LF[b]=j; } else BN[b]|=(uint16_t)(d<<j); } }
}

/* ------------------------------------------------------------------ one row */
static void row(long g,uint8_t*R,long*nunres,int*emax){
  double q[NS],x[NS]; I w[NS],W[4],D; g2p(g,q); long nu=0; int em=0;
  for(long c=0;c<65536;c++){
    for(int j=0;j<NS;j++) x[j]=(double)((c>>j)&1);
    int E,s; if(i_core(x,x,q,w,W,NULL,&E,&D,&s)){ fprintf(stderr,"i_core error\n"); exit(3); }
    double lo,hi; int cl; imargin(W,s,+1,&lo,&hi,&cl);
    int f;
    if(cl==1) f=2; else if(cl==-1) f=1; else if(cl==0) f=0;
    else if(up(D.hi*fmax(fabs(lo),fabs(hi)))<1.0) f=0;
    else { f=3; nu++; }
    if(E<0||E>63){ fprintf(stderr,"E out of range %d\n",E); exit(4); }
    if(E>em) em=E;
    R[c]=(uint8_t)((E<<2)|f);
  }
  *nunres=nu; *emax=em;
}

/* ------------------------------------------------------------------ ternary DP */
static void dp(const uint8_t*R,uint8_t*V){
  for(long hi=0;hi<6561;hi++){ long base=hi*6561; int jh=LF[hi];
    for(long lo=0;lo<6561;lo++){ long idx=base+lo; int j=LF[lo];
      if(j<0) j=(jh<0)?-1:8+jh;
      if(j<0){ V[idx]=R[(long)BN[lo]|((long)BN[hi]<<8)]; continue; }
      uint8_t a=V[idx-2*P3[j]], b=V[idx-P3[j]]; int ea=a>>2, eb=b>>2;
      V[idx]=(ea<eb)?a:((eb<ea)?b:(uint8_t)(a|b));
    } }
}

/* ------------------------------------------------------------------ classes */
static uint8_t*CLS;
static void cls_run(void*v,long a,long b){ (void)v; double p[NS],W4[4]; int s,E;
  for(long idx=a;idx<b;idx++){ long x=idx; for(int j=0;j<NS;j++){ int d=x%3; x/=3; p[j]=(d==2)?0.5:(double)d; }
    f_self(p,W4,&s,&E);
    int r[4]={0,0,0,0}; for(int k=0;k<NS;k++) if(s>>k&1) r[k>>2]=1;
    uint8_t c;
    if(!(r[1]||r[2]||r[3])) c=1; else if(!(r[0]||r[3])) c=2; else if(!r[3]) c=3; else c=0;
    CLS[idx]=c; } }


/* ------------------------------------------------------------------ verify (independent route) */
/* For each triple (stratum index, sign, q): enumerate the 2^d corners directly with i_core (no DP),
 * take E_sigma from f_core at the centre (free coordinates 1/2), and evaluate the margin at the
 * centre with the interval engine.  Result codes:
 *    1  every corner has E_c >= E_sigma, min_c E_c = E_sigma, every corner with E_c = E_sigma has
 *       margin <= 0 (certified < 0, structural zero or integrality zero), at least one is < 0, AND
 *       the centre margin enclosure is strictly negative;
 *    2  as 1 but the centre enclosure contains 0 (to be decided exactly by the caller);
 *   -1  the corner criterion fails;  -2  exponent inconsistency (L3 violated);
 *   -3  corner criterion holds but the centre enclosure is >= 0 or a structural zero (L4 violated). */
typedef struct { const int64_t*T; long n; int8_t*res; double*chi; int32_t*emin; } VJ;
static void ver_run(void*vv,long a,long b){ VJ*c=(VJ*)vv; double q[NS],x[NS],w16[NS]; I w[NS],W[4],D;
  for(long i=a;i<b;i++){ long idx=(long)c->T[3*i]; int sign=(int)c->T[3*i+1]; long g=(long)c->T[3*i+2];
    int fr[16],d=0; long y=idx; for(int j=0;j<NS;j++){ int dg=y%3; y/=3; if(dg==2){ fr[d++]=j; x[j]=0.5; } else x[j]=(double)dg; }
    g2p(g,q);
    int Es,ss; double Df; f_core(x,q,w16,NULL,&Es,&Df,&ss);
    I cw[NS],CW[4],CD; int Ec0,sc0; i_core(x,x,q,cw,CW,NULL,&Ec0,&CD,&sc0);
    double clo,chi; int ccl; imargin(CW,sc0,sign,&clo,&chi,&ccl);
    int emn=1000,bad=0,nneg=0,npos=0;
    for(long m=0;m<(1L<<d);m++){ for(int t=0;t<d;t++) x[fr[t]]=(double)((m>>t)&1);
      int E,su; i_core(x,x,q,w,W,NULL,&E,&D,&su); if(E<emn) emn=E;
      if(E<Es){ bad=1; break; } if(E>Es) continue;
      double lo,hi; int cl; imargin(W,su,sign,&lo,&hi,&cl);
      if(cl==2) cl=(up(D.hi*fmax(fabs(lo),fabs(hi)))<1.0)?0:9;
      if(cl==1||cl==9) npos++; if(cl==-1) nneg++; }
    int r;
    if(bad||emn!=Es||Ec0!=Es) r=-2;
    else if(npos||!nneg) r=-1;
    else if(ccl==-1) r=1; else if(ccl==2) r=2; else r=-3;
    c->res[i]=(int8_t)r; c->chi[i]=chi; c->emin[i]=Es; } }


/* ------------------------------------------------------------------ sets (for a minimal cover) */
typedef struct { const long*qs; long a,b; const int64_t*T; long n; uint8_t*bits; long nb; } SJ;
static void* sets_run(void*vv){ SJ*c=(SJ*)vv; uint8_t*R=malloc(65536); uint8_t*V=malloc(N3);
  for(long i=c->a;i<c->b;i++){ long nu; int em; row(c->qs[i],R,&nu,&em); if(nu){ fprintf(stderr,"unresolved in row %ld\n",c->qs[i]); exit(5); }
    dp(R,V); uint8_t*B=c->bits+i*c->nb; memset(B,0,c->nb);
    for(long t=0;t<c->n;t++){ int f=V[c->T[3*t]]&3; int ok=(c->T[3*t+1]>0)?(f==1):(f==2); if(ok) B[t>>3]|=(uint8_t)(1<<(t&7)); } }
  free(R); free(V); return NULL; }

/* ------------------------------------------------------------------ cover */
typedef struct { const long*qs; long nq; long a,b; const uint32_t*tidx; long nt; int32_t*qp,*qm; FILE*log; pthread_mutex_t*mu; double t0; } CJ;
static void amin(int32_t*p,int32_t v){ int32_t o=__atomic_load_n(p,__ATOMIC_RELAXED);
  while(o<0||v<o){ if(__atomic_compare_exchange_n(p,&o,v,0,__ATOMIC_RELAXED,__ATOMIC_RELAXED)) break; } }
static void* cover_run(void*vv){ CJ*c=(CJ*)vv; uint8_t*R=malloc(65536); uint8_t*V=malloc(N3);
  if(!R||!V){ fprintf(stderr,"malloc\n"); exit(2); }
  for(long i=c->a;i<c->b;i++){ long g=c->qs[i]; long nu; int em; double t1=now();
    row(g,R,&nu,&em); double t2=now(); dp(R,V); double t3=now();
    long np=0,nm=0;
    for(long t=0;t<c->nt;t++){ int f=V[c->tidx[t]]&3;
      if(f==1){ np++; amin(&c->qp[t],(int32_t)g); } else if(f==2){ nm++; amin(&c->qm[t],(int32_t)g); } }
    double t4=now();
    pthread_mutex_lock(c->mu);
    fprintf(c->log,"%ld %ld %ld unres %ld emax %d  row %.3f dp %.3f mark %.3f s\n",g,np,nm,nu,em,t2-t1,t3-t2,t4-t3);
    fflush(c->log); pthread_mutex_unlock(c->mu); }
  free(R); free(V); return NULL; }

static void* readall(const char*fn,long*n){ FILE*f=fopen(fn,"rb"); if(!f) return NULL; fseek(f,0,SEEK_END); long s=ftell(f); fseek(f,0,SEEK_SET);
  void*b=malloc(s>0?s:1); if(fread(b,1,s,f)!=(size_t)s){ fprintf(stderr,"read %s\n",fn); exit(2);} fclose(f); *n=s; return b; }
static void writeall(const char*fn,const void*b,long s){ FILE*f=fopen(fn,"wb"); if(!f||fwrite(b,1,s,f)!=(size_t)s){ fprintf(stderr,"write %s\n",fn); exit(2);} fclose(f); }

int main(int argc,char**argv){
  init_tables();
  if(argc>=3&&!strcmp(argv[1],"classes")){
    int nth=argc>3?atoi(argv[3]):3; CLS=malloc(N3); double t0=now();
    par(cls_run,NULL,N3,nth);
    long cnt[4]={0,0,0,0}; for(long i=0;i<N3;i++) cnt[CLS[i]]++;
    writeall(argv[2],CLS,N3);
    printf("classes: none %ld C %ld A %ld Z %ld  (%.1f s, %d threads)\n",cnt[0],cnt[1],cnt[2],cnt[3],now()-t0,nth);
    return 0;
  }
  if(argc>=5&&!strcmp(argv[1],"cover")){
    long n; uint8_t*cls=readall(argv[2],&n); if(!cls||n!=N3){ fprintf(stderr,"bad class file\n"); return 2; }
    int nth=argc>5?atoi(argv[5]):3;
    FILE*fq=fopen(argv[3],"r"); if(!fq){ fprintf(stderr,"no qlist\n"); return 2; }
    long cap=70000,nq=0; long*qs=malloc(cap*sizeof(long)); long g;
    while(fscanf(fq,"%ld",&g)==1){ if(g<0||g>65535){ fprintf(stderr,"bad q %ld\n",g); return 2; } qs[nq++]=g; } fclose(fq);
    long nt=0; for(long i=0;i<N3;i++) if(cls[i]) nt++;
    uint32_t*tidx=malloc(nt*4); long k=0; for(long i=0;i<N3;i++) if(cls[i]) tidx[k++]=(uint32_t)i;
    char fn[4096]; int32_t*qp,*qm; long s1,s2,s0;
    snprintf(fn,sizeof fn,"%s.idx",argv[4]); uint32_t*old=readall(fn,&s0);
    if(old){ if(s0!=nt*4||memcmp(old,tidx,s0)){ fprintf(stderr,"existing %s differs\n",fn); return 2; } free(old);
      snprintf(fn,sizeof fn,"%s.qp",argv[4]); qp=readall(fn,&s1); snprintf(fn,sizeof fn,"%s.qm",argv[4]); qm=readall(fn,&s2);
      if(!qp||!qm||s1!=nt*4||s2!=nt*4){ fprintf(stderr,"bad existing qp/qm\n"); return 2; } }
    else { writeall(fn,tidx,nt*4); qp=malloc(nt*4); qm=malloc(nt*4); for(long t=0;t<nt;t++){ qp[t]=-1; qm[t]=-1; } }
    snprintf(fn,sizeof fn,"%s.log",argv[4]); FILE*log=fopen(fn,"a");
    pthread_mutex_t mu; pthread_mutex_init(&mu,NULL);
    if(nth>64) nth=64; if(nth>nq) nth=(int)(nq>0?nq:1);
    pthread_t th[64]; CJ cj[64]; double t0=now();
    for(int i=0;i<nth;i++){ cj[i]=(CJ){qs,nq,nq*i/nth,nq*(i+1)/nth,tidx,nt,qp,qm,log,&mu,t0}; pthread_create(&th[i],NULL,cover_run,&cj[i]); }
    for(int i=0;i<nth;i++) pthread_join(th[i],NULL);
    fclose(log);
    snprintf(fn,sizeof fn,"%s.qp",argv[4]); writeall(fn,qp,nt*4); snprintf(fn,sizeof fn,"%s.qm",argv[4]); writeall(fn,qm,nt*4);
    long cp=0,cm=0; for(long t=0;t<nt;t++){ cp+=qp[t]>=0; cm+=qm[t]>=0; }
    printf("cover: %ld co-players, %ld efficient strata; excluded sign+ %ld, sign- %ld (cumulative); %.1f s on %d threads\n",nq,nt,cp,cm,now()-t0,nth);
    return 0;
  }
  if(argc>=3&&!strcmp(argv[1],"row")){
    uint8_t*R=malloc(65536); long nu; int em; double t0=now(); row(atol(argv[2]),R,&nu,&em);
    long np=0,nn=0,nz=0; for(long c=0;c<65536;c++){ int f=R[c]&3; np+=f==2; nn+=f==1; nz+=f==0; }
    printf("q %s: pos %ld neg %ld zero %ld unres %ld emax %d (%.3f s)\n",argv[2],np,nn,nz,nu,em,now()-t0);
    if(argc>=4) writeall(argv[3],R,65536);
    return 0;
  }
  if(argc>=5&&!strcmp(argv[1],"check")){
    const char*s=argv[2]; int sign=atoi(argv[3]); long g=atol(argv[4]); if(strlen(s)!=16) return 2;
    int fr[16],d=0; double x[NS],q[NS]; I w[NS],W[4],D; g2p(g,q);
    for(int j=0;j<16;j++){ if(s[j]=='*') fr[d++]=j; x[j]=(s[j]=='1')?1.0:0.0; }
    int emin=1000; long npos=0,nneg=0,nhigh=0;
    int*Ec=malloc(sizeof(int)<<d); int*fl=malloc(sizeof(int)<<d);
    for(long m=0;m<(1L<<d);m++){ for(int t=0;t<d;t++) x[fr[t]]=(double)((m>>t)&1);
      int E,su; i_core(x,x,q,w,W,NULL,&E,&D,&su); double lo,hi; int cl; imargin(W,su,sign,&lo,&hi,&cl);
      if(cl==2) cl=(up(D.hi*fmax(fabs(lo),fabs(hi)))<1.0)?0:9;
      Ec[m]=E; fl[m]=cl; if(E<emin) emin=E; }
    for(long m=0;m<(1L<<d);m++){ if(Ec[m]>emin){ nhigh++; continue; } if(fl[m]==1||fl[m]==9) npos++; if(fl[m]==-1||fl[m]==9) nneg++; }
    printf("%s sign %d q %ld: d %d E %d corners pos %ld neg %ld high %ld -> %s\n",s,sign,g,d,emin,npos,nneg,nhigh,(npos==0&&nneg>0)?"EXCLUDED":"not excluded");
    return 0;
  }
  if(argc>=4&&!strcmp(argv[1],"verify")){
    long n; int64_t*T=readall(argv[2],&n); if(!T||n%24){ fprintf(stderr,"bad triples\n"); return 2; } n/=24;
    int nth=argc>4?atoi(argv[4]):3; VJ c={T,n,malloc(n),malloc(n*8),malloc(n*4)}; double t0=now();
    par(ver_run,&c,n,nth);
    long k[7]={0}; for(long i=0;i<n;i++) k[c.res[i]+3]++;
    char fn[4096]; snprintf(fn,sizeof fn,"%s.res",argv[3]); writeall(fn,c.res,n);
    snprintf(fn,sizeof fn,"%s.chi",argv[3]); writeall(fn,c.chi,n*8);
    snprintf(fn,sizeof fn,"%s.E",argv[3]); writeall(fn,c.emin,n*4);
    printf("verify: %ld triples: ok %ld, centre-unresolved %ld, corner-fail %ld, exponent-fail %ld, centre-contradiction %ld (%.1f s, %d threads)\n",
           n,k[4],k[5],k[2],k[1],k[0],now()-t0,nth);
    return 0;
  }
  if(argc>=5&&!strcmp(argv[1],"sets")){
    FILE*fq=fopen(argv[2],"r"); long nq=0; long*qs=malloc(70000*sizeof(long)); long g;
    while(fq&&fscanf(fq,"%ld",&g)==1) qs[nq++]=g; if(fq) fclose(fq);
    long n; int64_t*T=readall(argv[3],&n); n/=24; long nb=(n+7)/8; uint8_t*bits=calloc(nq*nb,1);
    int nth=argc>5?atoi(argv[5]):3; if(nth>nq) nth=(int)nq; pthread_t th[64]; SJ sj[64]; double t0=now();
    for(int i=0;i<nth;i++){ sj[i]=(SJ){qs,nq*i/nth,nq*(i+1)/nth,T,n,bits,nb}; pthread_create(&th[i],NULL,sets_run,&sj[i]); }
    for(int i=0;i<nth;i++) pthread_join(th[i],NULL);
    writeall(argv[4],bits,nq*nb); printf("sets: %ld co-players x %ld triples (%.1f s)\n",nq,n,now()-t0); return 0;
  }
  if(argc>=4&&!strcmp(argv[1],"rowcmp")){
    /* compare a row file (e.g. from exact arithmetic) with the interval row of q: ub rowcmp Q FILE */
    long n; uint8_t*X=readall(argv[3],&n); if(!X||n!=65536) return 2; uint8_t*R=malloc(65536); long nu; int em;
    row(atol(argv[2]),R,&nu,&em); long diff=0; for(long c=0;c<65536;c++) diff+=(X[c]!=R[c]);
    printf("rowcmp q %s: %ld of 65536 entries differ, %ld unresolved\n",argv[2],diff,nu); return diff?1:0;
  }
  fprintf(stderr,"usage: ub classes CLS.bin [nth] | ub cover CLS.bin QLIST OUT [nth] | ub row Q | ub check PAT SIGN Q\n");
  return 1;
}
