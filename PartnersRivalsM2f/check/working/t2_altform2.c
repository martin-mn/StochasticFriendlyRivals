/* altform2.c -- (1) validate the Freidlin-Wentzell reduction SUPPORT-WISE against freqs16,
 *               (2) test the reachability closed form for Eff_alt, Eff_line, Eff_CC.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define NS 16
void freqs16(const double*p,const double*q,double w[NS]);
static inline int swapo(int o){ return ((o&1)<<1)|(o>>1); }
static inline int mirror_state(int j){ return 4*swapo(j>>2) + swapo(j&3); }
static uint64_t rs;
static inline double urnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0; }
#define COREMASK 0x0FF0u
static int ord[NS][4],nxt[NS][4]; static unsigned reach[NS],adj0[NS];

static int analyse(const int*d,unsigned cls[NS],int kesc[NS],unsigned land[NS]){
  for(int s=0;s<NS;s++){ int a=d[s],b=d[mirror_state(s)];
    int fC=(a==0),fD=(a==2),cC=(b==0),cD=(b==2);
    int co[4]={fC+cC,fC+cD,fD+cC,fD+cD};
    for(int o=0;o<4;o++){ord[s][o]=co[o];nxt[s][o]=4*o+(s>>2);} }
  for(int s=0;s<NS;s++){adj0[s]=0;for(int o=0;o<4;o++)if(ord[s][o]==0)adj0[s]|=1u<<nxt[s][o];}
  for(int s=0;s<NS;s++) reach[s]=adj0[s]|(1u<<s);
  for(int it=0;it<NS;it++){int ch=0;for(int s=0;s<NS;s++){unsigned r=reach[s];
      for(int t=0;t<NS;t++)if(r>>t&1)r|=reach[t]; if(r!=reach[s]){reach[s]=r;ch=1;}} if(!ch)break;}
  int n=0;unsigned done=0;
  for(int s=0;s<NS;s++){ if(done>>s&1)continue; unsigned scc=0;
    for(int t=0;t<NS;t++) if((reach[s]>>t&1)&&(reach[t]>>s&1)) scc|=1u<<t;
    done|=scc; unsigned out=0; for(int t=0;t<NS;t++) if(scc>>t&1) out|=adj0[t];
    if(out&~scc) continue;
    int k=99;unsigned L=0;
    for(int t=0;t<NS;t++) if(scc>>t&1) for(int o=0;o<4;o++){ int u=nxt[t][o]; if(scc>>u&1)continue;
      if(ord[t][o]<k){k=ord[t][o];L=1u<<u;} else if(ord[t][o]==k) L|=1u<<u; }
    cls[n]=scc;kesc[n]=k;land[n]=L;n++; }
  return n;
}
int main(int argc,char**argv){
  long stripe=atol(argv[1]),nstr=atol(argv[2]);
  long TOT=1;for(int i=0;i<16;i++)TOT*=3;
  long suppMis=0, altMis=0, lineMis=0, ccMis=0, nA=0,nL=0,nC=0;
  long isolated=0, isolAlt=0;
  double p[NS],w[NS];int d[16];unsigned cls[NS],land[NS];int kesc[NS];
  for(long g=stripe;g<TOT;g+=nstr){
    long x=g;for(int j=0;j<16;j++){d[j]=x%3;x/=3;}
    int n=analyse(d,cls,kesc,land);
    unsigned cadj[NS];
    for(int i=0;i<n;i++){cadj[i]=0;for(int t=0;t<NS;t++)if(land[i]>>t&1)for(int j2=0;j2<n;j2++)
        if(reach[t]&cls[j2]) cadj[i]|=1u<<j2;}
    unsigned cr[NS];for(int i=0;i<n;i++)cr[i]=cadj[i]|(1u<<i);
    for(int it=0;it<n;it++){int ch=0;for(int i=0;i<n;i++){unsigned r=cr[i];
        for(int j2=0;j2<n;j2++)if(r>>j2&1)r|=cr[j2]; if(r!=cr[i]){cr[i]=r;ch=1;}} if(!ch)break;}
    unsigned recmask=0,cdone=0;int nrec=0;
    for(int i=0;i<n;i++){ if(cdone>>i&1)continue; unsigned scc=0;
      for(int j2=0;j2<n;j2++) if((cr[i]>>j2&1)&&(cr[j2]>>i&1)) scc|=1u<<j2;
      cdone|=scc; unsigned out=0;for(int j2=0;j2<n;j2++)if(scc>>j2&1)out|=cadj[j2];
      if(out&~scc)continue; recmask|=scc;nrec++; }
    unsigned supp=0;for(int i=0;i<n;i++)if(recmask>>i&1)supp|=cls[i];
    /* truth support */
    rs=(uint64_t)g*6364136223846793005ULL+1442695040888963407ULL;if(!rs)rs=88172645463325252ULL;
    for(int j=0;j<16;j++)p[j]=d[j]==0?0.0:(d[j]==2?1.0:0.12+0.76*urnd());
    freqs16(p,p,w);
    unsigned tsupp=0;for(int j=0;j<NS;j++) if(w[j]>1e-13) tsupp|=1u<<j;
    if(tsupp!=supp){ if(suppMis<10) fprintf(stderr,"SUPP g=%ld pred=%04x true=%04x\n",g,supp,tsupp); suppMis++; }
    /* ---- closed forms ---- */
    /* alternating closed classes present? */
    unsigned altcls=0; int haveAlt=0;
    for(int i=0;i<n;i++) if((cls[i]&~COREMASK)==0){haveAlt=1;altcls|=cls[i];}
    /* NONALT = union of non-alternating closed classes */
    unsigned nonalt=0; for(int i=0;i<n;i++) if(cls[i]&~COREMASK) nonalt|=cls[i];
    int safe1 = !(reach[1]&nonalt), safe13 = !(reach[13]&nonalt);
    int predA = haveAlt && safe1 && safe13;
    int truthA = ((tsupp&~COREMASK)==0);
    if(truthA)nA++;
    if(predA!=truthA){ if(altMis<10)fprintf(stderr,"ALT g=%ld pred=%d truth=%d\n",g,predA,truthA); altMis++; }
    /* isolated non-core recurrent SCC not reachable from core?  count for diagnostics */
    if(haveAlt && safe1 && safe13 && (supp&~COREMASK)) isolated++;
    if(haveAlt) isolAlt++;
    /* Eff_CC : support in states 0..3 ; Eff_line : support avoids 12..15 */
    int truthC = ((tsupp&~0x000Fu)==0); if(truthC)nC++;
    int truthL = !(tsupp&0xF000u); if(truthL)nL++;
    /* predicted: p0=1 and from states 4,8 (escapes of {0}) only class {0} reachable */
    int predC = (d[0]==2) && !(reach[4]&~1u & ~0u ? 0:0);
    { unsigned non0=0; for(int i=0;i<n;i++) if(cls[i]!=1u) non0|=cls[i];
      predC = (d[0]==2) && !(reach[4]&non0) && !(reach[8]&non0); }
    if(predC!=truthC){ if(ccMis<10)fprintf(stderr,"CC g=%ld pred=%d truth=%d\n",g,predC,truthC); ccMis++; }
  }
  fprintf(stderr,"stripe %ld: suppMis %ld | altMis %ld (truthA %ld) | ccMis %ld (truthC %ld) truthL %ld | haveAlt %ld isolatedNonCore %ld\n",
          stripe,suppMis,altMis,nA,ccMis,nC,nL,isolAlt,isolated);
  return 0;
}
