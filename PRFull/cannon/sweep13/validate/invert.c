/* invert.c (scratchpad): recover the psweep seed of an archived draw from one interior coordinate.
 * The interior value v = k/2^53 fixes the top 53 bits of the xorshift64 state; the 2^11 candidates are
 * filtered by the draw's 16 pin decisions, then the state is walked back to a state of the form
 * seed*M + C with a small seed.   usage: invert <alpha_permille> <maxback> v0 ... v15 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static inline uint64_t fwd(uint64_t s){ s^=s<<13; s^=s>>7; s^=s<<17; return s; }
static inline uint64_t unl(uint64_t y,int a){ uint64_t x=y; for(int k=a;k<64;k+=a) x^=y<<k; return x; }
static inline uint64_t unr(uint64_t y,int a){ uint64_t x=y; for(int k=a;k<64;k+=a) x^=y>>k; return x; }
static inline uint64_t bwd(uint64_t s){ s=unl(s,17); s=unr(s,7); s=unl(s,13); return s; }
static inline double u(uint64_t s){ return (double)((s>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }
int main(int argc,char**argv){
  double al=atoi(argv[1])/1000.0; long maxback=atol(argv[2]); double v[16]; int jint=-1;
  for(int j=0;j<16;j++){ v[j]=strtod(argv[3+j],0); if(v[j]!=0.0&&v[j]!=1.0&&jint<0) jint=j; }
  const uint64_t M=2862933555777941757ULL, C=3037000493ULL;
  uint64_t Minv=M; for(int i=0;i<6;i++) Minv*=2-M*Minv;   /* Newton: inverse mod 2^64 */
  uint64_t k=(uint64_t)(v[jint]*9007199254740992.0);
  for(uint64_t lo=0;lo<2048;lo++){
    uint64_t sb=(k<<11)|lo;                       /* state whose urnd() is the interior value of coord jint */
    /* replay the whole draw forward from the state before coord 0: find it by walking back */
    uint64_t s=bwd(sb); if(!(u(s)>=2*al)) continue; /* the pin decision of coord jint: interior */
    /* walk back through coords jint-1..0 (each: pin decision, plus a value if interior -- none here
       unless another interior coord precedes; general case handled by the forward check below) */
    uint64_t st=s; int ok=1;
    for(int j=jint-1;j>=0;j--){ st=bwd(st); double uu=u(st);
      if(v[j]==0.0){ if(!(uu<al)) {ok=0;break;} } else if(v[j]==1.0){ if(!(uu>=al&&uu<2*al)) {ok=0;break;} } else {ok=0;break;} }
    if(!ok) continue;
    uint64_t s0=bwd(st);                           /* state before the draw */
    /* forward check of the whole draw */
    uint64_t f=s0; for(int j=0;j<16&&ok;j++){ f=fwd(f); double uu=u(f);
      if(v[j]==0.0) ok=uu<al; else if(v[j]==1.0) ok=(uu>=al&&uu<2*al); else { ok=uu>=2*al; f=fwd(f); ok=ok&&(u(f)==v[j]); } }
    if(!ok) continue;
    printf("candidate low bits %llu: draw starts after state %016llx\n",(unsigned long long)lo,(unsigned long long)s0);
    uint64_t t=s0; for(long b=0;b<=maxback;b++){ uint64_t seed=(t-C)*Minv;
      if(seed<100000000ULL){ printf("  seed %llu, %ld RNG calls before the draw\n",(unsigned long long)seed,b); }
      t=bwd(t); }
  }
  return 0; }
