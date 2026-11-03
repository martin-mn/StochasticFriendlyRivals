/* cover.c -- which of a given list of co-players excludes (L3 + L4, exactly as wsep's try_q) each
 * stratum of a pattern file.  Output per stratum: pattern, then the success bitmask over the list
 * as 4 hex words (list of at most 256 co-players).  Used for a greedy set cover of the certificate.
 * usage: cover <qfile> <patternfile> <nthreads<=3> <outfile> */
#define WSEP_NO_MAIN
#include "wsep.c"
typedef struct { signed char (*pat)[16]; long n; long*q; int nq; uint64_t (*mask)[4]; long next; pthread_mutex_t mu; } CV;
static void*cvw(void*v){ CV*c=(CV*)v;
  for(;;){ pthread_mutex_lock(&c->mu); long a=c->next; c->next+=256; pthread_mutex_unlock(&c->mu);
    if(a>=c->n) break; long b=a+256; if(b>c->n) b=c->n;
    for(long i=a;i<b;i++){ int fr[16],d=0; for(int j=0;j<16;j++) if(c->pat[i][j]==2) fr[d++]=j;
      for(int k=0;k<c->nq;k++){ Try T; try_q(c->pat[i],fr,d,c->q[k],&T); if(T.ok) c->mask[i][k>>6]|=1ULL<<(k&63); } } }
  return NULL; }
int main(int argc,char**argv){
  if(argc<5){ fprintf(stderr,"usage: cover <qfile> <patternfile> <nthreads> <outfile>\n"); return 2; }
  static CV c; memset(&c,0,sizeof c); c.q=malloc(sizeof(long)*256);
  FILE*f=fopen(argv[1],"r"); long g; while(c.nq<256 && fscanf(f,"%ld",&g)==1) c.q[c.nq++]=g; fclose(f);
  f=fopen(argv[2],"r"); char buf[256]; long cap=1<<16; c.pat=malloc(cap*16);
  while(fgets(buf,sizeof buf,f)){ if(strlen(buf)<16) continue; if(c.n==cap){cap*=2; c.pat=realloc(c.pat,cap*16);}
    for(int j=0;j<16;j++) c.pat[c.n][j]=(buf[j]=='*')?2:(buf[j]-'0'); c.n++; }
  fclose(f); c.mask=calloc(c.n,sizeof *c.mask); pthread_mutex_init(&c.mu,NULL);
  int nth=atoi(argv[3]); if(nth<1) nth=1; if(nth>3) nth=3; double t0=now();
  pthread_t th[3]; for(int i=0;i<nth;i++) pthread_create(&th[i],NULL,cvw,&c); for(int i=0;i<nth;i++) pthread_join(th[i],NULL);
  FILE*o=fopen(argv[4],"w"); char ps[17];
  for(long i=0;i<c.n;i++){ patstr(c.pat[i],ps); fprintf(o,"%s\t%016llx\t%016llx\t%016llx\t%016llx\n",ps,
      (unsigned long long)c.mask[i][0],(unsigned long long)c.mask[i][1],(unsigned long long)c.mask[i][2],(unsigned long long)c.mask[i][3]); }
  fclose(o); printf("strata %ld coplayers %d seconds %.1f\n",c.n,c.nq,now()-t0); return 0; }
