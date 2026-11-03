/* wfc.c -- the Wright-Fisher process with pairwise comparison over the CONTINUOUS memory-two cube
 *          [0,1]^16, with every sampled resident classified EXACTLY in the limit eps -> 0.
 *
 * THE PROCESS is that of DiskM2WF/Cannon/dw1/sf.f, the code behind Figures 4d-f and 5d-f of
 * PartnersRivals1/Revised2 (and that of Nowak, arXiv:2608.06147): N individuals, one random
 * co-player each (m = 1), payoffs from the eps-perturbed sixteen-state pair chain solved by GTH,
 * synchronous update, mutation with probability mu, otherwise Fermi imitation of a random role
 * model (itself included) with intensity beta/(2m) on the accumulated payoff, the second half of
 * each replicate sampled.  The generation loop follows sf.f statement for statement.  Exactly three
 * things change, because they must:
 *
 *   (1) MUTATION.  A mutant's sixteen coordinates are drawn independently from the three-outcome
 *       law of arXiv:2608.06147 SI 3.1: 0 with probability alpha, 1 with probability alpha, uniform
 *       on (0,1) with probability 1 - 2 alpha; the pair chain then plays the clamped strategy
 *       eps + (1 - 2 eps) p, so the atoms sit at p_min = eps and p_max = 1 - eps as in that paper.
 *       alpha = 1/2 is uniform over the 65536 binary genomes, i.e. EXACTLY sf.f's kernel (imut = 1),
 *       which is how this code is validated against the binary run; alpha = 1/3 is the continuum run.
 *       The initial population is drawn from the same law (at alpha = 1/2, sf.f's 16 fair coins).
 *
 *   (2) STRATEGIES are vectors, not 16-bit codes: a table of live strategies with reference counts,
 *       each carrying a never-reused 32-bit serial.  The pair cache (2^20 slots, sf.f's hash) is
 *       keyed by the ordered serial pair, so a cached value can never belong to a dead strategy.
 *
 *   (3) ABUNDANCE is sampled, not histogrammed.  sf.f adds every individual of every sampled
 *       generation into a 65536-bin histogram; over the continuum there is no finite histogram.
 *       Instead NT sample times are placed evenly over the sampled half, with a random offset, and
 *       at each NI individuals are drawn uniformly.  The fraction of samples in an atom is then an
 *       unbiased estimator of that atom's share of the individual-generations of the half, the
 *       quantity sf.f's histogram measures, with binomial error at most 1/(2 sqrt(NT NI)) per
 *       replicate.  The efficiency is NOT sampled: it is accumulated exactly, as in sf.f.
 *
 * and one thing is added: after each replicate the DISTINCT sampled strategies are classified in the
 * limit eps -> 0 by this project's engine m2f.c (the classification of Revised2, carried to the
 * continuum):
 *     efficient   self-play limit W_CC = 1 below the switch line, W_CD = W_DC = 1/2 above it
 *     competitive min over the 65536 binary co-players of sign(T-S)(W_DC - W_CD) >= 0   (box-vertex
 *                 lemma: the worst co-player over the whole cube is binary)
 *     Nash        max over the 65536 binary co-players q of pi(q,p) <= pi(p,p)          (the same
 *                 lemma for the co-player's best reply: an MDP whose optimum over the clamped box
 *                 is at a clamped vertex; the eps -> 0 limit commutes with the finite max)
 *     comparisons with tolerance 1e-9, as Revised2.  An efficient competitive strategy is Nash by
 *     the theorem of Revised2, so its Nash scan is skipped (atom 111) -- except for the first
 *     CHECK111 of them in each task, which are scanned anyway and must come out Nash.
 *
 * Modes
 *   ./wfc run  <task> <alpha_permille> <itend> <NT> <NI> [nrep] [seedbase] [rep0]
 *        one task = one game of games.dat (line <task>, 1-based) and replicates rep0..rep0+nrep-1
 *        (default 1..10), replicate r seeded seedbase + 10(task-1) + r.  A replicate's result does
 *        not depend on how many share its task: the cache decides only whether a pair value is
 *        recomputed, never what it is, and serials keep their order under any offset.
 *   ./wfc binall <task> 500
 *        every binary genome at that game, one byte each (B<task>.bin): the check against
 *        PartnersRivals/Cannon/atoms512.bin, the exact atom table of Revised2
 *   ./wfc neutral <task> <alpha_permille> <M> [seed]
 *        the alpha-measure of the seven atoms at that game: M independent draws from the law,
 *        each classified -- the continuum analogue of |a|/65536, the enrichment's denominator
 *
 * Output (run):  stdout, one line per replicate
 *   REP gid ipt u v Emax pay ef alpha irep nsamp ndist a000 a001 a010 a011 a100 a101 a110 a111
 *       nhit nmis nnashscan nrivscan t_sim t_class  inW ex111 nfr111 ncand111 total
 *   inW = 1 at games of wedge W, where ex111 is the EXACT number of individual-generations of the sampled half
 *   held by friendly rivals (nfr111 distinct of them, ncand111 candidates tested) out of total = (itend/2) N
 *   where aXYZ counts samples with efficient=X, Nash=Y, competitive=Z (a101 must be 0);
 *   and a binary file S<task>.bin holding every distinct sampled strategy of every replicate:
 *   int32 irep, int32 count, int32 atom, double p[16]   -- for re-analysis without re-simulating.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "witnesses.h"

#define NS 16
void freqs(const double*p,const double*q,double W[4]);
double rival_margin(const double*p,int sign,long*arg,double earlyquit);

/* ------------------------------------------------------------------ parameters of sf.f
 * The defaults are the large-population run (N = 1000, beta = 100, mu = 1e-2: Revised2's F1, Figure 4).  The small-population
 * run F2 (N = 100, beta = 3, mu = 1e-4, 1e8 generations: Revised2's SI Figure 4) is the same source compiled with
 *     -DNPOP=100 -DBETAV=3.0 -DMUV=1e-4
 * and nothing else changes; every output file begins with a '#' line naming the values it was built with. */
#ifndef NPOP
#define NPOP 1000
#endif
#ifndef BETAV
#define BETAV 100.0
#endif
#ifndef MUV
#define MUV 1e-2
#endif
#define N     NPOP
#define MPL   1
static const double BETA = BETAV, MU = MUV, EPS = 1e-4, THR = 25.0;
#define CHECK111 20

/* ------------------------------------------------------------------ RNG: xoshiro256** */
static uint64_t X0,X1,X2,X3;
static inline uint64_t rotl(uint64_t x,int k){ return (x<<k)|(x>>(64-k)); }
static inline uint64_t nxt(void){
  uint64_t r=rotl(X1*5,7)*9, t=X1<<17;
  X2^=X0; X3^=X1; X1^=X2; X0^=X3; X2^=t; X3=rotl(X3,45); return r; }
static void seedrng(uint64_t s){
  uint64_t *d[4]={&X0,&X1,&X2,&X3};
  for(int k=0;k<4;k++){ s+=0x9E3779B97F4A7C15ull; uint64_t z=s;
    z=(z^(z>>30))*0xBF58476D1CE4E5B9ull; z=(z^(z>>27))*0x94D049BB133111EBull; *d[k]=z^(z>>31); } }
/* strictly inside (0,1): k < 2^52 so k + 0.5 is exact and the largest value is 1 - 2^-53.  (The 53-bit form
   (k + 0.5) 2^-53 rounds to exactly 1.0 at k = 2^53 - 1, which would let ridx(n) return n.) */
static inline double u01(void){ return ((double)(nxt()>>12)+0.5)*(1.0/4503599627370496.0); }
static inline int ridx(int n){ return (int)(u01()*n); }

/* ------------------------------------------------------------------ the game */
static double R_=1, S_, T_, EMAX, UU, VV;

/* ------------------------------------------------------------------ strategy table */
static inline int swapo(int o){ return ((o&1)<<1)|(o>>1); }
static inline int mirror_state(int j){ return 4*swapo(j>>2)+swapo(j&3); }
#define MAXS (2*N+8)
static double RAW[MAXS][NS];      /* p in [0,1]^16, pinned coordinates exactly 0.0 or 1.0 */
static double CLP[MAXS][NS];      /* eps + (1-2eps) p, focal view                          */
static double CLM[MAXS][NS];      /* the same, read at the mirrored state: the co-player view */
static uint32_t SER[MAXS]; static int REF[MAXS]; static int FREEL[MAXS]; static int nfree;
static uint32_t nextser;
static double ALPHA;
/* EXACT abundance of atom 111 in wedge W.  There 111 is rare both in the alpha-measure (3.9e-6) and in the
 * population, so 4000 samples per replicate cannot resolve it and the most-enriched map would be decided by 0-3
 * sampled friendly rivals.  A friendly rival in W needs p[0] = 1 (W_CC = 1 in self-play) and p[15] = 0 (every
 * T>S rival has it), both exact necessary conditions; a mutant meeting both gets the exact self-play test at birth
 * (C111 = 1 if efficient), its individual-generations in the sampled half are counted exactly (AB), and at its
 * death the exact rivalry test decides whether AB goes to ex111.  Every friendly rival is Nash (Revised2's
 * theorem), so ex111 is the exact 111 count, the quantity sf.f's histogram gives. */
static int GAMEW;
static unsigned char C111[MAXS]; static long long AB[MAXS];
static long long ex111, ncand111, nfr111;

static void table_reset(void){ nfree=0; for(int k=MAXS-1;k>=0;k--){ FREEL[nfree++]=k; REF[k]=0; } }
static void drawalpha(double*p){
  for(int j=0;j<NS;j++){ double r=u01();
    p[j] = r<ALPHA ? 0.0 : (r<2.0*ALPHA ? 1.0 : u01()); } }
static int newmutant(void){
  if(nfree<=0){ fprintf(stderr,"strategy table full\n"); exit(3); }
  int id=FREEL[--nfree]; SER[id]=nextser++; REF[id]=0;
  drawalpha(RAW[id]);
  for(int j=0;j<NS;j++) CLP[id][j]=EPS+(1.0-2.0*EPS)*RAW[id][j];
  for(int j=0;j<NS;j++) CLM[id][j]=CLP[id][mirror_state(j)];
  C111[id]=0; AB[id]=0;
  if(GAMEW && RAW[id][0]==1.0 && RAW[id][15]==0.0){ double W[4]; freqs(RAW[id],RAW[id],W); C111[id]=(W[0]>1-1e-12); }
  return id; }
static void settle111(int id);
static inline void release(int id){ if(--REF[id]==0){ if(C111[id]) settle111(id); FREEL[nfree++]=id; } }

/* ------------------------------------------------------------------ eps-perturbed pair chain, GTH
 * payf3 of sf.f for arbitrary (not binary) strategies: state s = 4*recent + before; from s the chain
 * moves to 4*onew + (s>>2); outcome onew = 2*(focal D) + (co-player D). */
static void pair_eps(int a,int b,double w[4]){
  double g[NS][NS]; memset(g,0,sizeof g);
  const double *x=CLP[a], *y=CLM[b];
  for(int s=0;s<NS;s++){ int r=s>>2; double xs=x[s], ys=y[s];
    g[s][r]=xs*ys; g[s][4+r]=xs*(1-ys); g[s][8+r]=(1-xs)*ys; g[s][12+r]=(1-xs)*(1-ys); }
  for(int n=NS-1;n>=1;n--){
    double s=0; for(int j=0;j<n;j++) s+=g[n][j];
    double inv=1.0/s;
    for(int i=0;i<n;i++){ double pin=g[i][n];
      if(pin!=0.0){ pin*=inv; g[i][n]=pin; for(int j=0;j<n;j++) g[i][j]+=pin*g[n][j]; } } }
  double v[NS], t=1.0; v[0]=1.0;
  for(int n=1;n<NS;n++){ double s=0; for(int i=0;i<n;i++) s+=v[i]*g[i][n]; v[n]=s; t+=s; }
  t=1.0/t;
  w[0]=(v[0]+v[1]+v[2]+v[3])*t;   w[1]=(v[4]+v[5]+v[6]+v[7])*t;
  w[2]=(v[8]+v[9]+v[10]+v[11])*t; w[3]=(v[12]+v[13]+v[14]+v[15])*t;
}

/* ------------------------------------------------------------------ pair cache (sf.f's hash) */
#define LCB 20
#define NCB (1u<<LCB)
static uint64_t CK[NCB]; static double CVL[NCB], CVH[NCB];
static long long nhit,nmis;
static void cache_clear(void){ for(uint32_t k=0;k<NCB;k++) CK[k]=UINT64_MAX; }

/* ------------------------------------------------------------------ the eps -> 0 classification */
static double WP[NWIT][NS];
static void initwit(void){ for(int k=0;k<NWIT;k++) for(int j=0;j<NS;j++) WP[k][j]=(double)((WITG[k]>>j)&1); }
static int rivscreen(const double*p,int sign){
  double W[4];
  for(int k=0;k<NWIT;k++){ freqs(p,WP[k],W); if(sign*(W[2]-W[1]) < -1e-9) return 0; }
  return 1; }
/* recent Nash violators, tried first: consecutive samples at one game are refuted by the same few */
#define NVW 64
static unsigned VW[NVW]; static int nvw;
static void vw_push(unsigned g){
  for(int k=0;k<nvw;k++) if(VW[k]==g){ for(int m=k;m>0;m--) VW[m]=VW[m-1]; VW[0]=g; return; }
  if(nvw<NVW) nvw++;
  for(int m=nvw-1;m>0;m--) VW[m]=VW[m-1]; VW[0]=g; }
static inline double payq(const double*q,const double*p){ double W[4]; freqs(q,p,W);
  return R_*W[0]+S_*W[1]+T_*W[2]; }
static int nash_scan(const double*p,double pp){
  double q[NS];
  for(int k=0;k<nvw;k++){ unsigned g=VW[k]; for(int j=0;j<NS;j++) q[j]=(double)((g>>j)&1);
    if(payq(q,p)>pp+1e-9){ vw_push(g); return 0; } }
  for(unsigned g=0;g<65536u;g++){ for(int j=0;j<NS;j++) q[j]=(double)((g>>j)&1);
    if(payq(q,p)>pp+1e-9){ vw_push(g); return 0; } }
  return 1; }
static long nnashscan, nrivscan, n111checked;
static void settle111(int id){
  if(AB[id]>0){ ncand111++; long a;
    if(rivscreen(RAW[id],+1) && rival_margin(RAW[id],+1,&a,-1e-9)>=-1e-9){ ex111+=AB[id]; nfr111++; } }
  C111[id]=0; AB[id]=0; }
/* atom = 4*efficient + 2*Nash + competitive */
static int classify(const double*p){
  double W[4]; freqs(p,p,W);
  int below=(UU+VV<1.0);
  int eff = below ? (W[0]>1-1e-12) : (W[1]+W[2]>1-1e-12);
  int sign=(UU-VV<1.0)?+1:-1;
  int riv=0; long a;
  if(rivscreen(p,sign)){ nrivscan++; riv = rival_margin(p,sign,&a,-1e-9) >= -1e-9; }
  double pp=R_*W[0]+S_*W[1]+T_*W[2];
  int nash;
  if(eff&&riv){
    nash=1;
    if(n111checked<CHECK111){ n111checked++; nnashscan++;
      if(!nash_scan(p,pp)){ fprintf(stderr,"THEOREM VIOLATED: efficient competitive, not Nash\n"); exit(4); } }
  } else { nnashscan++; nash=nash_scan(p,pp); }
  return 4*eff+2*nash+riv; }

/* ------------------------------------------------------------------ sample table (dedupe by serial) */
#define LSH 16
#define NSH (1u<<LSH)
static uint32_t SHK[NSH]; static int SHC[NSH]; static double SHP[NSH][NS]; static int nsh;
static void samp_clear(void){ for(uint32_t k=0;k<NSH;k++){ SHK[k]=UINT32_MAX; SHC[k]=0; } nsh=0; }
static void samp_add(int id){
  uint32_t k=(SER[id]*2654435761u)>>(32-LSH);
  while(SHK[k]!=UINT32_MAX && SHK[k]!=SER[id]) k=(k+1)&(NSH-1);
  if(SHK[k]==UINT32_MAX){ SHK[k]=SER[id]; memcpy(SHP[k],RAW[id],sizeof(SHP[k])); nsh++;
    if(nsh>NSH/2){ fprintf(stderr,"sample table half full: lower NT*NI\n"); exit(5); } }
  SHC[k]++; }

/* ------------------------------------------------------------------ games.dat */
static int read_game(int task,int*gid,int*ipt){
  FILE*f=fopen("games.dat","r"); if(!f){ perror("games.dat"); exit(2); }
  char line[512]; int k=0;
  while(fgets(line,sizeof line,f)){
    if(line[0]=='#') continue;
    int g,ic; double cr,cs,ct,u,v; char tag[64];
    if(sscanf(line,"%d %d %lf %lf %lf %lf %lf %63s",&g,&ic,&cr,&cs,&ct,&u,&v,tag)!=8) continue;
    if(++k==task){ fclose(f);
      if(cr!=1.0){ fprintf(stderr,"expected cR = 1\n"); exit(2); }
      *gid=g; R_=cr; S_=cs; T_=ct; UU=cs; VV=ct-1.0;
      if(sscanf(tag,"V_p%d",ipt)!=1) *ipt=-1;
      EMAX = fmax(R_,0.5*(S_+T_)); return 0; } }
  fclose(f); fprintf(stderr,"task %d beyond games.dat\n",task); exit(2); }

static double now(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+1e-9*ts.tv_nsec; }

int main(int argc,char**argv){
  if(argc<4){ fprintf(stderr,"usage: wfc run <task> <alpha_pm> <itend> <NT> <NI> [nrep] [seedbase]\n"
                              "       wfc neutral <task> <alpha_pm> <M> [seed]\n"); return 1; }
  initwit();
  int task=atoi(argv[2]);
  { const char*sl=strchr(argv[3],'/');                /* "1/3" exactly, or permille as before ("500") */
    ALPHA = sl ? atof(argv[3])/atof(sl+1) : atoi(argv[3])/1000.0; }
  if(ALPHA<=0.0||ALPHA>0.5){ fprintf(stderr,"alpha must be in (0, 0.5]\n"); return 1; }
  /* binall ignores alpha; pass 500 */
  int gid,ipt; read_game(task,&gid,&ipt);
  GAMEW = (UU+VV<1.0) && (UU-VV<1.0);

  if(!strcmp(argv[1],"binall")){
    /* every binary genome at this game, one byte each: the validation against atoms512.bin */
    static unsigned char out[65536]; double p[NS]; double t0=now();
    for(unsigned g=0;g<65536u;g++){ for(int j=0;j<NS;j++) p[j]=(double)((g>>j)&1); out[g]=(unsigned char)classify(p); }
    char fn[64]; snprintf(fn,sizeof fn,"B%d.bin",task); FILE*fo=fopen(fn,"wb"); fwrite(out,1,65536,fo); fclose(fo);
    long cnt[8]={0}; for(unsigned g=0;g<65536u;g++) cnt[out[g]]++;
    printf("BIN %d %d",gid,ipt); for(int a=0;a<8;a++) printf(" %ld",cnt[a]); printf(" %.1f\n",now()-t0);
    return 0; }
  if(!strcmp(argv[1],"neutral")){
    long M=atol(argv[4]); uint64_t sd = argc>5 ? strtoull(argv[5],0,10) : 777000000ull+task;
    seedrng(sd); long cnt[8]={0}; double p[NS]; double t0=now();
    for(long k=0;k<M;k++){ drawalpha(p); cnt[classify(p)]++; }
    printf("NEU %d %d %.12g %.12g %.12g %.4f %ld",gid,ipt,UU,VV,EMAX,ALPHA,M);
    for(int a=0;a<8;a++) printf(" %ld",cnt[a]);
    printf(" %ld %ld %.1f\n",nnashscan,nrivscan,now()-t0);
    return 0; }

  if(strcmp(argv[1],"run")||argc<7){ fprintf(stderr,"bad mode/arguments\n"); return 1; }
  long itend=atol(argv[4]); int NT=atoi(argv[5]), NI=atoi(argv[6]);
  int nrep = argc>7 ? atoi(argv[7]) : 10;
  uint64_t sbase = argc>8 ? strtoull(argv[8],0,10) : 9000000ull;
  int rep0 = argc>9 ? atoi(argv[9]) : 1;
  const char *spath = argc>10 ? argv[10] : NULL;       /* where the strategy dump goes: named by the caller,
                                                          so tasks at different alpha can never share a file */      /* first replicate index: a 1-replicate task can run
                                                 replicate r with exactly the seed it has in a 10-replicate one */
  long init=itend/2; const double be=BETA*0.5/MPL;

  printf("# wfc N=%d beta=%g mu=%g eps=%g alpha=%.10f itend=%ld NT=%d NI=%d task=%d reps=%d..%d\n",
         N,BETA,MU,EPS,ALPHA,itend,NT,NI,task,rep0,rep0+nrep-1);
  char bn[256]; if(spath) snprintf(bn,sizeof bn,"%s",spath); else snprintf(bn,sizeof bn,"S%d_%d.bin",task,rep0);
  FILE*fb=fopen(bn,"wb"); if(!fb){ perror(bn); return 2; }
  static int s[N], sn[N]; static double f[N];
  cache_clear(); nhit=nmis=0; nextser=0; n111checked=0;

  for(int irep=rep0;irep<rep0+nrep;irep++){
    seedrng(sbase+(uint64_t)(task-1)*10+irep);   /* stride 10 whatever nrep is */
    table_reset(); samp_clear(); nvw=0;
    for(int k=0;k<MAXS;k++){ C111[k]=0; AB[k]=0; }
    ex111=ncand111=nfr111=0;
    long ns0=nnashscan, nr0=nrivscan;
    for(int i=0;i<N;i++){ s[i]=newmutant(); REF[s[i]]++; }
    double z=0, ef=0; double t0=now();
    /* sample times: NT evenly spaced over the second half, random offset */
    long H=itend-init; double off=u01(); long nexts=0; long tnext=init+1+(long)floor(off*(double)H/NT);

    for(long it=1;it<=itend;it++){
      for(int i=0;i<N;i++) f[i]=0.0;
      for(int i=0;i<N;i++){
        int si=s[i]; int j;
        do { j=ridx(N); } while(j==i);
        int sj=s[j];
        uint32_t a=SER[si], b=SER[sj]; int isw = (a>b);
        int ia = isw? sj:si, ib = isw? si:sj;
        uint64_t key = ((uint64_t)SER[ia]<<32)|SER[ib];
        uint64_t ih = key^(key>>17); ih*=1000003ull; ih^=ih>>11;
        uint32_t slot=(uint32_t)(ih&(NCB-1));
        double vl,vh;
        if(CK[slot]==key){ vl=CVL[slot]; vh=CVH[slot]; nhit++; }
        else { double w[4]; pair_eps(ia,ib,w);
          vl=R_*w[0]+S_*w[1]+T_*w[2]; vh=R_*w[0]+T_*w[1]+S_*w[2];
          CK[slot]=key; CVL[slot]=vl; CVH[slot]=vh; nmis++; }
        if(!isw){ f[i]+=vl; f[j]+=vh; } else { f[i]+=vh; f[j]+=vl; }
      }
      if(it>init){
        z+=1.0; double efg=0; for(int i=0;i<N;i++){ efg+=f[i]; if(C111[s[i]]) AB[s[i]]++; } ef+=efg;
        while(nexts<NT && it==tnext){
          for(int k=0;k<NI;k++) samp_add(s[ridx(N)]);
          nexts++; tnext=init+1+(long)floor(((double)nexts+off)*(double)H/NT);
        }
      }
      for(int i=0;i<N;i++){
        if(u01()<MU){ sn[i]=newmutant(); }
        else {
          int j=ridx(N);
          double xf=be*(f[j]-f[i]), prob;
          if(xf>=THR) prob=1.0; else if(xf<=-THR) prob=0.0;
          else if(xf>=0.0) prob=1.0/(1.0+exp(-xf)); else { double ex=exp(xf); prob=ex/(1.0+ex); }
          sn[i] = (u01()<prob) ? s[j] : s[i];
        }
        REF[sn[i]]++;
      }
      for(int i=0;i<N;i++){ release(s[i]); s[i]=sn[i]; }
    }
    double tsim=now()-t0;
    for(int i=0;i<N;i++) release(s[i]);

    double pay=ef/(double)N*0.5/MPL/z, effi=pay/EMAX;
    long cnt[8]={0}; long nsamp=0; double t1=now();
    for(uint32_t k=0;k<NSH;k++) if(SHK[k]!=UINT32_MAX){
      int at=classify(SHP[k]); cnt[at]+=SHC[k]; nsamp+=SHC[k];
      int32_t hdr[3]={irep,SHC[k],at}; fwrite(hdr,sizeof hdr,1,fb); fwrite(SHP[k],sizeof(double),NS,fb); }
    double tcl=now()-t1;
    printf("REP %d %d %.12g %.12g %.12g %.12g %.10f %.4f %d %ld %d",
           gid,ipt,UU,VV,EMAX,pay,effi,ALPHA,irep,nsamp,nsh);
    int order[8]={0,1,2,3,4,5,6,7};
    for(int k=0;k<8;k++) printf(" %ld",cnt[order[k]]);
    printf(" %lld %lld %ld %ld %.1f %.1f %d %lld %lld %lld %.0f\n",nhit,nmis,nnashscan-ns0,nrivscan-nr0,tsim,tcl,
           GAMEW,ex111,nfr111,ncand111,z*(double)N);
    fflush(stdout);
    if(cnt[5]){ fprintf(stderr,"atom 101 nonempty\n"); return 4; }
  }
  fclose(fb);
  return 0;
}
