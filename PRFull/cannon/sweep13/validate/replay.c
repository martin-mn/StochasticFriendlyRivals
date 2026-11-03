/* replay.c (scratchpad validation tool): re-draw the kernel stream of a psweep task with the PATCHED
 * census.c's own kdraw() and locate each archived W point (pt file, %.17g, exact doubles) in it.
 * usage: replay <n> <seed> <alpha> <ptfile> [maxpts]     prints: seed npts nfound idx_1 idx_2 ... */
#define main census_main
#include "../census.c"
#undef main
int main(int argc,char**argv){
  long n=atol(argv[1]); long seed=atol(argv[2]); long am; char tag[48];
  if(parse_alpha(argv[3],&am,tag,sizeof tag)) return 2;
  static double W[20000][NS]; int nw=0; char buf[4096];
  FILE*f=fopen(argv[4],"r"); if(!f) return 3;
  while(fgets(buf,sizeof buf,f)&&nw<20000){ if(buf[0]!='W') continue; char*s=buf+1;
    for(int j=0;j<NS;j++) W[nw][j]=strtod(s,&s); nw++; }
  fclose(f);
  rs=(uint64_t)seed*2862933555777941757ULL+3037000493ULL; if(!rs)rs=7;
  double p[NS]; long npin=0; int next=0; long idx[20000];
  for(long i=0;i<n && next<nw;i++){ kdraw(p,&npin);
    if(!memcmp(p,W[next],sizeof p)){ idx[next++]=i; } }
  printf("%ld %d %d",seed,nw,next); for(int k=0;k<next;k++) printf(" %ld",idx[k]); printf("\n");
  return next==nw?0:1; }
