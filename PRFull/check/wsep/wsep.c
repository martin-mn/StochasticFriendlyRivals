/* wsep.c -- exact separation certificate for the friendly rivals of wedge W (referee issue 2).
 *
 * Reuses the verified leading-order engines of ../cert/certeng.c verbatim (#include, no copy), so
 * every number below comes from the same f_core / i_core / imargin code that cert/tests.py checks.
 *
 * For every support stratum of a candidate class it decides, with no tolerance anywhere:
 *   (1) efficiency below the switch line: the self-play limit support lies in the states whose last
 *       outcome is CC (support class 'C').  Exponent-only, exact, constant on the open stratum (L1).
 *   (2) exclusion by ALLD: against ALLD (genome 0) w_DC + w_CC = eps at every eps, so the limit margin
 *       is -w_CD <= 0; it is < 0 on the whole open stratum iff a state with last outcome CD
 *       (states 4..7) lies in the limit support -- again exponent-only and exact (L1).
 *   (3) exclusion by one binary co-player q (L3 + L4): every corner c of the closed stratum cube has
 *       E_c(q) > E(q) (so N_q(c) = 0) or a binary margin certified <= 0 (strictly negative interval,
 *       structural zero, or zero by integrality of D* x margin), at least one corner is certified
 *       strictly negative, and the margin at the centre (all free coordinates 1/2, exact in binary)
 *       is certified strictly negative by the outward-rounded interval engine.  Then N_q < 0 on the
 *       whole open cube, so q strictly out-earns every strategy of the open stratum where T > S.
 *
 * Classes (engine convention: 0, 1, 2 = free; p_0 = 1 throughout, which efficiency forces):
 *   c6  : p6 free                                   (3^14 strata)
 *   c8  : p8 free, p6 in {0,1}                       (2*3^13)
 *   c68 : (p6,p8) in {(0,1),(1,0)}                   (2*3^13)
 *   c9a : p6 = p8 = 0, p9 in {0, free}               (2*3^12)
 *   c9b : p6 = p8 = 1, p9 in {0, 1}                  (2*3^12)
 *   c0  : p0 in {0, free} (control: never efficient) (2*3^15; efficiency only)
 *   file: patterns read from a file (control strata; no pool, full search)
 *   t:TOKENS : a template of 16 tokens, each 0, 1, * (fixed), . (any) or a set like [1*]
 *
 * usage: wsep <class|file:PATH|t:TEMPLATE> <nthreads> <outfile> [stride] [offset]
 * Output: one line per efficient stratum not excluded by ALLD:
 *   cls pattern d E q how nneg nzs nzi nhigh centre_lo centre_hi ntried
 * (how = P pool, F full search, N none found), and a summary block on stdout.
 */
#include "../cert/certeng.c"
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

#define MAXPOOL 512
typedef struct { int nspec; int spec[4]; int nassign; int assign[4][4]; int nrest; int rest[16]; long total;
                 int tmpl; int nv[16]; int vals[16][3]; } Cls;
/* template class "t:TOKENS": 16 tokens, each 0, 1, * (fixed), . (any of 0, 1, free) or a set such as [1*] */
static int mktmpl(const char*t,Cls*C){
  memset(C,0,sizeof(*C)); C->tmpl=1; int j=0; const char*p=t;
  while(*p && j<16){
    if(*p=='['){ p++; while(*p && *p!=']'){ C->vals[j][C->nv[j]++]=(*p=='*')?2:(*p-'0'); p++; } if(*p!=']') return -1; p++; j++; }
    else if(*p=='.'){ C->nv[j]=3; C->vals[j][0]=0; C->vals[j][1]=1; C->vals[j][2]=2; p++; j++; }
    else if(*p=='0'||*p=='1'||*p=='*'){ C->nv[j]=1; C->vals[j][0]=(*p=='*')?2:(*p-'0'); p++; j++; }
    else return -1;
  }
  if(j!=16||*p) return -1;
  C->total=1; for(j=0;j<16;j++){ if(C->nv[j]<1||C->nv[j]>3) return -1; C->total*=C->nv[j]; }
  return 0;
}

static double now(void){ struct timeval t; gettimeofday(&t,NULL); return t.tv_sec+1e-6*t.tv_usec; }
static long ipow3(int n){ long r=1; while(n--) r*=3; return r; }

static int mkcls(const char*name,Cls*C){
  memset(C,0,sizeof(*C));
  int sp[4]={0,0,0,0}; int ns=0;
  #define A(...) do{ int v[]={__VA_ARGS__}; for(int t=0;t<ns;t++) C->assign[C->nassign][t]=v[t]; C->nassign++; }while(0)
  if(!strcmp(name,"c6")){ ns=2; sp[0]=0; sp[1]=6; C->nspec=ns; memcpy(C->spec,sp,sizeof sp); A(1,2); }
  else if(!strcmp(name,"c8")){ ns=3; sp[0]=0; sp[1]=6; sp[2]=8; C->nspec=ns; memcpy(C->spec,sp,sizeof sp); A(1,0,2); A(1,1,2); }
  else if(!strcmp(name,"c68")){ ns=3; sp[0]=0; sp[1]=6; sp[2]=8; C->nspec=ns; memcpy(C->spec,sp,sizeof sp); A(1,0,1); A(1,1,0); }
  else if(!strcmp(name,"c9a")){ ns=4; sp[0]=0; sp[1]=6; sp[2]=8; sp[3]=9; C->nspec=ns; memcpy(C->spec,sp,sizeof sp); A(1,0,0,0); A(1,0,0,2); }
  else if(!strcmp(name,"c9b")){ ns=4; sp[0]=0; sp[1]=6; sp[2]=8; sp[3]=9; C->nspec=ns; memcpy(C->spec,sp,sizeof sp); A(1,1,1,0); A(1,1,1,1); }
  else if(!strcmp(name,"c0")){ ns=1; sp[0]=0; C->nspec=ns; memcpy(C->spec,sp,sizeof sp); A(0); A(2); }
  else return -1;
  #undef A
  int used[16]={0}; for(int t=0;t<C->nspec;t++) used[C->spec[t]]=1;
  for(int j=0;j<16;j++) if(!used[j]) C->rest[C->nrest++]=j;
  C->total=C->nassign*ipow3(C->nrest);
  return 0;
}
static void decode(const Cls*C,long i,signed char*pat){
  if(C->tmpl){ for(int j=0;j<16;j++){ pat[j]=(signed char)C->vals[j][i%C->nv[j]]; i/=C->nv[j]; } return; }
  int a=(int)(i%C->nassign); long r=i/C->nassign;
  for(int t=0;t<C->nspec;t++) pat[C->spec[t]]=(signed char)C->assign[a][t];
  for(int t=0;t<C->nrest;t++){ pat[C->rest[t]]=(signed char)(r%3); r/=3; }
}
static void patstr(const signed char*pat,char*s){ for(int j=0;j<16;j++) s[j]="01*"[(int)pat[j]]; s[16]=0; }

/* SI eq:si-alld closed form (as a stratum condition): competitive against ALLD iff
 * p15 = 0, p13 < 1, and (p5 < 1 or p7 = 0 = p6 p11).  Used only for a cross-check count. */
static int si_alld_ok(const signed char*p){
  if(p[15]!=0) return 0; if(p[13]==1) return 0;
  if(p[5]!=1) return 1;
  return (p[7]==0) && (p[6]==0 || p[11]==0);
}

typedef struct { int E; int ok; int nneg,nzs,nzi,nhigh,nlow; double clo,chi; int why; } Try;
/* why: 0 ok, 1 centre not certified negative, 2 positive corner, 3 unresolved corner, 4 E_c<E (lemma violated),
 *      5 no negative corner although the centre is negative (lemma violated) */
static void try_q(const signed char*pat,const int*fr,int d,long g,Try*T){
  double q[NS],x[NS]; I w[NS],W[4],D; int E,s;
  g2p(g,q); memset(T,0,sizeof(*T));
  for(int j=0;j<NS;j++) x[j]=(pat[j]==2)?0.5:(double)pat[j];
  i_core(x,x,q,w,W,NULL,&E,&D,&s); T->E=E;
  int cl; imargin(W,s,+1,&T->clo,&T->chi,&cl);
  if(cl!=-1){ T->why=1; return; }
  for(long m=0;m<(1L<<d);m++){
    for(int t=0;t<d;t++) x[fr[t]]=(double)((m>>t)&1);
    int Ec,sc; i_core(x,x,q,w,W,NULL,&Ec,&D,&sc);
    if(Ec>E){ T->nhigh++; continue; }
    if(Ec<E){ T->nlow++; T->why=4; return; }
    double lo,hi; int c2; imargin(W,sc,+1,&lo,&hi,&c2);
    if(c2==1){ T->why=2; return; }
    if(c2==-1){ T->nneg++; continue; }
    if(c2==0){ T->nzs++; continue; }
    if(up(D.hi*fmax(fabs(lo),fabs(hi)))<1.0){ T->nzi++; continue; }
    T->why=3; return;
  }
  if(T->nneg==0){ T->why=5; return; }
  T->ok=1;
}

/* ------------------------------------------------------------------ shared run state */
typedef struct {
  const Cls*C; int isfile; signed char (*fpats)[16]; long nfile;
  long total,stride,offset; int nth; FILE*out; pthread_mutex_t mu; long next; int effonly;
  /* totals */
  long n,neff,nalld,nsurv,npool,nfull,nnone,nsi_dis,nlemma,byd_eff[17],byd_surv[17],byd_n[17];
  long poolq[65536];      /* hits per co-player (pool or full) */
  double tfull;
} Run;

static int cmp_d2(const void*a,const void*b){ const double*x=a,*y=b; return (x[0]<y[0])?-1:(x[0]>y[0]); }

static void*worker(void*v){
  Run*R=(Run*)v; long pool[MAXPOOL]; int np=0; char ps[17];
  double*cand=malloc(sizeof(double)*2*65536);
  long n=0,neff=0,nalld=0,nsurv=0,npool=0,nfull=0,nnone=0,nsi=0,nlem=0,be[17]={0},bs[17]={0},bn[17]={0};
  double tfull=0; static const long CH=2048;
  for(;;){
    pthread_mutex_lock(&R->mu); long a=R->next; R->next+=CH; pthread_mutex_unlock(&R->mu);
    long lim=R->isfile?R->nfile:R->total; if(a>=lim) break;
    long b=a+CH; if(b>lim) b=lim;
    for(long i=a;i<b;i++){
      if(!R->isfile && ((i-R->offset)%R->stride)!=0) continue;
      signed char pat[16]; if(R->isfile) memcpy(pat,R->fpats[i],16); else decode(R->C,i,pat);
      int fr[16],d=0; for(int j=0;j<16;j++) if(pat[j]==2) fr[d++]=j;
      double x[NS],w16[NS],q0[NS]; for(int j=0;j<NS;j++){ x[j]=(pat[j]==2)?0.5:(double)pat[j]; q0[j]=0; }
      n++; bn[d]++;
      int E,s; double Df; f_core(x,x,w16,NULL,&E,&Df,&s);
      if(s & 0xFFF0) continue;                      /* not efficient below the switch line */
      neff++; be[d]++;
      if(R->effonly) continue;
      f_core(x,q0,w16,NULL,&E,&Df,&s);
      int alld_excl=(s & 0x00F0)!=0;
      if(alld_excl==si_alld_ok(pat)) nsi++;         /* disagreement with the SI closed form */
      if(alld_excl){ nalld++; continue; }
      nsurv++; bs[d]++;
      Try T; long qfound=-1; char how='N'; int ntried=0;
      if(!R->isfile){
        for(int k=0;k<np;k++){ ntried++; try_q(pat,fr,d,pool[k],&T);
          if(T.why==4||T.why==5) nlem++;
          if(T.ok){ qfound=pool[k]; how='P'; long t=pool[k]; for(int u=k;u>0;u--) pool[u]=pool[u-1]; pool[0]=t; break; } }
      }
      if(qfound<0){
        double t0=now();
        /* full search: every q with a certified-negative centre margin, most negative first */
        double q[NS]; I ww[NS],W[4],D; int nc=0;
        for(long g=0;g<65536;g++){ g2p(g,q); int e2,s2; i_core(x,x,q,ww,W,NULL,&e2,&D,&s2);
          double lo,hi; int cl; imargin(W,s2,+1,&lo,&hi,&cl); if(cl==-1){ cand[2*nc]=hi; cand[2*nc+1]=(double)g; nc++; } }
        qsort(cand,nc,2*sizeof(double),cmp_d2);
        for(int k=0;k<nc;k++){ long g=(long)cand[2*k+1]; ntried++; try_q(pat,fr,d,g,&T);
          if(T.why==4||T.why==5) nlem++;
          if(T.ok){ qfound=g; how='F'; break; } }
        if(qfound>=0 && !R->isfile){ if(np<MAXPOOL) np++; for(int u=np-1;u>0;u--) pool[u]=pool[u-1]; pool[0]=qfound; }
        tfull+=now()-t0;
      }
      if(how=='P') npool++; else if(how=='F') nfull++; else nnone++;
      patstr(pat,ps);
      pthread_mutex_lock(&R->mu);
      if(qfound>=0) R->poolq[qfound]++;
      fprintf(R->out,"%s\t%s\t%d\t%d\t%ld\t%c\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%d\n",
              R->isfile?"file":"c", ps, d, T.E, qfound, how, T.nneg, T.nzs, T.nzi, T.nhigh,
              qfound>=0?T.clo:0.0, qfound>=0?T.chi:0.0, ntried);
      pthread_mutex_unlock(&R->mu);
    }
  }
  pthread_mutex_lock(&R->mu);
  R->n+=n; R->neff+=neff; R->nalld+=nalld; R->nsurv+=nsurv; R->npool+=npool; R->nfull+=nfull; R->nnone+=nnone;
  R->nsi_dis+=nsi; R->nlemma+=nlem; R->tfull+=tfull;
  for(int k=0;k<17;k++){ R->byd_eff[k]+=be[k]; R->byd_surv[k]+=bs[k]; R->byd_n[k]+=bn[k]; }
  pthread_mutex_unlock(&R->mu);
  free(cand); return NULL;
}

#ifndef WSEP_NO_MAIN
int main(int argc,char**argv){
  if(argc<4){ fprintf(stderr,"usage: wsep <class|file:PATH|t:TEMPLATE> <nthreads> <outfile> [stride] [offset] [effonly]\n"); return 2; }
  static Run R; memset(&R,0,sizeof R); Cls C;
  R.nth=atoi(argv[2]); if(R.nth<1) R.nth=1; if(R.nth>3 && !getenv("WSEP_ALLOW_MORE")) R.nth=3;
  R.stride=(argc>4)?atol(argv[4]):1; R.offset=(argc>5)?atol(argv[5]):0; R.effonly=(argc>6)?atoi(argv[6]):0;
  if(!strncmp(argv[1],"file:",5)){
    FILE*f=fopen(argv[1]+5,"r"); if(!f){ perror("pattern file"); return 2; }
    char buf[256]; long cap=1024; R.fpats=malloc(cap*16); R.nfile=0;
    while(fgets(buf,sizeof buf,f)){ if(strlen(buf)<16) continue; if(R.nfile==cap){ cap*=2; R.fpats=realloc(R.fpats,cap*16); }
      for(int j=0;j<16;j++) R.fpats[R.nfile][j]=(buf[j]=='*')?2:(buf[j]-'0'); R.nfile++; }
    fclose(f); R.isfile=1; R.total=R.nfile;
  } else if(!strncmp(argv[1],"t:",2)){
    if(mktmpl(argv[1]+2,&C)){ fprintf(stderr,"bad template %s\n",argv[1]); return 2; }
    R.C=&C; R.total=C.total;
  } else {
    if(mkcls(argv[1],&C)){ fprintf(stderr,"unknown class %s\n",argv[1]); return 2; }
    R.C=&C; R.total=C.total;
  }
  R.out=fopen(argv[3],"w"); if(!R.out){ perror("out"); return 2; }
  pthread_mutex_init(&R.mu,NULL);
  double t0=now();
  pthread_t th[3]; for(int i=0;i<R.nth;i++) pthread_create(&th[i],NULL,worker,&R);
  for(int i=0;i<R.nth;i++) pthread_join(th[i],NULL);
  fclose(R.out);
  double el=now()-t0;
  printf("class %s total %ld stride %ld processed %ld\n",argv[1],R.total,R.stride,R.n);
  printf("efficient %ld\nexcluded_by_ALLD %ld\nSI_alld_disagreements %ld\nsurvivors %ld\nexcluded_pool %ld\nexcluded_full %ld\nNOT_EXCLUDED %ld\nlemma_violations %ld\n",
         R.neff,R.nalld,R.nsi_dis,R.nsurv,R.npool,R.nfull,R.nnone,R.nlemma);
  printf("seconds %.2f (full-search %.2f thread-s) threads %d\n",el,R.tfull,R.nth);
  printf("by_d n/eff/surv:"); for(int k=0;k<17;k++) if(R.byd_n[k]) printf(" %d:%ld/%ld/%ld",k,R.byd_n[k],R.byd_eff[k],R.byd_surv[k]); printf("\n");
  printf("coplayers:"); int nq=0; for(long g=0;g<65536;g++) if(R.poolq[g]){ nq++; } printf(" %d distinct\n",nq);
  for(long g=0;g<65536;g++) if(R.poolq[g]) printf("q %ld %ld\n",g,R.poolq[g]);
  return 0;
}
#endif
