/* kthr.c (scratchpad): exact test of the patched kdraw() pin thresholds.  The xorshift state is set so that
 * the next k53() returns a chosen integer k (by inverting the xorshift step), then kdraw() is called and the
 * value of coordinate 0 read.  usage: kthr <alpha> */
#define main census_main
#include "../census.c"
#undef main
static inline uint64_t unl(uint64_t y,int a){ uint64_t x=y; for(int k=a;k<64;k+=a) x^=y<<k; return x; }
static inline uint64_t unr(uint64_t y,int a){ uint64_t x=y; for(int k=a;k<64;k+=a) x^=y>>k; return x; }
static inline uint64_t bwd(uint64_t s){ s=unl(s,17); s=unr(s,7); s=unl(s,13); return s; }
int main(int argc,char**argv){
  long am; char tag[48]; if(parse_alpha(argv[1],&am,tag,sizeof tag)) return 2;
  uint64_t Q=AQ, T[]={0,Q-1,Q,2*AN*Q-1,2*AN*Q,AD*Q-1,AD*Q,(1ull<<53)-1};
  const char*nm[]={"0","AQ-1","AQ","2AQ-1","2AQ","3AQ-1","3AQ (redraw)","2^53-1 (redraw)"};
  for(int t=0;t<8;t++){ rs=bwd((T[t]<<11)|0x5a5); long np=0,r0=NREDRAW; double p[NS]; kdraw(p,&np);
    printf("alpha %s  k = %-16s %16llu -> p0 = %-22.17g redraws %ld\n",tag,nm[t],(unsigned long long)T[t],p[0],NREDRAW-r0); }
  return 0; }
