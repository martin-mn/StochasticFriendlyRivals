/* fw.c -- EXACT Freidlin-Wentzell account of the eps->0 self-play limit support.
 *   V(s,t) = min total eps-order of a path s->t in the 16-state self-play graph
 *            (edge s -> 4o+(s>>2) costs the number of players that must err);
 *   closed classes K_1..K_n of the eps=0 chain;  c(i,j) = min_{s in K_i, t in K_j} V(s,t);
 *   W(j)   = min weight of an in-tree to K_j in the complete digraph on classes (Edmonds);
 *   SUPPORT of the limit = union of the K_j attaining min_j W(j).
 * Tested state-by-state against freqs16 on all 3^16 support strata.
 * Then: Eff_alt = support inside states 4..11;  Eff_line = support avoids 12..15;
 *       Eff_CC  = support inside states 0..3.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define NS 16
#define INF 1000000
void freqs16(const double*p,const double*q,double w[NS]);
static inline int swapo(int o){return((o&1)<<1)|(o>>1);}
static inline int mirror_state(int j){return 4*swapo(j>>2)+swapo(j&3);}
static uint64_t rs;
static inline double urnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return (double)((rs>>11)&0x1FFFFFFFFFFFFFull)/9007199254740992.0;}

/* ---- Edmonds minimum spanning arborescence rooted at 'root' (edges u->v) ---- */
#define MAXV 40
typedef struct{int u,v,w;} E;
static int edmonds(int n,E*ed,int m,int root){
  int res=0; int inw[MAXV],inu[MAXV],id[MAXV],vis[MAXV];
  E buf[MAXV*MAXV]; memcpy(buf,ed,sizeof(E)*m);
  while(1){
    for(int i=0;i<n;i++) inw[i]=INF,inu[i]=-1;
    for(int e=0;e<m;e++){int u=buf[e].u,v=buf[e].v,w=buf[e].w;
      if(u!=v&&w<inw[v]){inw[v]=w;inu[v]=u;}}
    for(int i=0;i<n;i++) if(i!=root&&inu[i]<0) return INF;
    int cnt=0; for(int i=0;i<n;i++) id[i]=-1,vis[i]=-1;
    inw[root]=0;
    for(int i=0;i<n;i++){
      res+=inw[i]; int v=i;
      while(vis[v]!=i&&id[v]<0&&v!=root){vis[v]=i;v=inu[v];}
      if(v!=root&&id[v]<0){
        for(int u=inu[v];u!=v;u=inu[u]) id[u]=cnt;
        id[v]=cnt++; }
    }
    if(cnt==0) return res;
    for(int i=0;i<n;i++) if(id[i]<0) id[i]=cnt++;
    int m2=0;
    for(int e=0;e<m;e++){int u=buf[e].u,v=buf[e].v,w=buf[e].w;
      if(id[u]!=id[v]){buf[m2].u=id[u];buf[m2].v=id[v];buf[m2].w=w-inw[v];m2++;} }
    m=m2; n=cnt; root=id[root];
  }
}

static int ord[NS][4],nxt[NS][4]; static unsigned reach[NS],adj0[NS]; static int V[NS][NS];
static int build(const int*d,unsigned cls[NS]){
  for(int s=0;s<NS;s++){int a=d[s],b=d[mirror_state(s)];
    int fC=(a==0),fD=(a==2),cC=(b==0),cD=(b==2); int co[4]={fC+cC,fC+cD,fD+cC,fD+cD};
    for(int o=0;o<4;o++){ord[s][o]=co[o];nxt[s][o]=4*o+(s>>2);} }
  for(int s=0;s<NS;s++){adj0[s]=0;for(int o=0;o<4;o++)if(ord[s][o]==0)adj0[s]|=1u<<nxt[s][o];}
  for(int s=0;s<NS;s++)reach[s]=adj0[s]|(1u<<s);
  for(int it=0;it<NS;it++){int ch=0;for(int s=0;s<NS;s++){unsigned r=reach[s];
    for(int t=0;t<NS;t++)if(r>>t&1)r|=reach[t]; if(r!=reach[s]){reach[s]=r;ch=1;}} if(!ch)break;}
  for(int i=0;i<NS;i++)for(int j=0;j<NS;j++)V[i][j]=(i==j)?0:INF;
  for(int s=0;s<NS;s++)for(int o=0;o<4;o++){int t=nxt[s][o];if(ord[s][o]<V[s][t])V[s][t]=ord[s][o];}
  for(int k=0;k<NS;k++)for(int i=0;i<NS;i++)for(int j=0;j<NS;j++)
    if(V[i][k]+V[k][j]<V[i][j])V[i][j]=V[i][k]+V[k][j];
  int n=0;unsigned done=0;
  for(int s=0;s<NS;s++){if(done>>s&1)continue;unsigned scc=0;
    for(int t=0;t<NS;t++)if((reach[s]>>t&1)&&(reach[t]>>s&1))scc|=1u<<t;
    done|=scc;unsigned out=0;for(int t=0;t<NS;t++)if(scc>>t&1)out|=adj0[t];
    if(out&~scc)continue; cls[n++]=scc;}
  return n;
}
int main(int argc,char**argv){
  long stripe=atol(argv[1]),nstr=atol(argv[2]);int verbose=argc>3?atoi(argv[3]):0;
  long TOT=1;for(int i=0;i<16;i++)TOT*=3;
  long suppMis=0,tA=0,tL=0,tC=0,pA=0,pL=0,pC=0,mA=0,mL=0,mC=0,maxn=0;
  double p[NS],w[NS];int d[16];unsigned cls[NS];
  for(long g=stripe;g<TOT;g+=nstr){
    long x=g;for(int j=0;j<16;j++){d[j]=x%3;x/=3;}
    int n=build(d,cls); if(n>maxn)maxn=n;
    int c[NS][NS];
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){int b=INF;
      if(i==j) b=0;
      else { for(int s=0;s<NS;s++) if(cls[i]>>s&1) for(int t=0;t<NS;t++) if(cls[j]>>t&1)
               if(V[s][t]<b) b=V[s][t]; }
      c[i][j]=b;}
    /* in-tree to root r  ==  arborescence rooted at r in the reversed graph */
    E ed[NS*NS];int m=0;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++) if(i!=j&&c[i][j]<INF){ed[m].u=j;ed[m].v=i;ed[m].w=c[i][j];m++;}
    int Wb=INF,Wv[NS];
    for(int r=0;r<n;r++){Wv[r]=edmonds(n,ed,m,r); if(Wv[r]<Wb)Wb=Wv[r];}
    unsigned supp=0;for(int r=0;r<n;r++)if(Wv[r]==Wb)supp|=cls[r];
    rs=(uint64_t)g*6364136223846793005ULL+1442695040888963407ULL;if(!rs)rs=88172645463325252ULL;
    for(int j=0;j<16;j++)p[j]=d[j]==0?0.0:(d[j]==2?1.0:0.12+0.76*urnd());
    freqs16(p,p,w);
    unsigned ts=0;for(int j=0;j<NS;j++)if(w[j]>1e-13)ts|=1u<<j;
    if(ts!=supp){suppMis++;if(verbose&&suppMis<25)fprintf(stderr,"SUPP g=%ld pred=%04x true=%04x n=%d\n",g,supp,ts,n);}
    unsigned TG[3]={0x0FF0u,0x0FFFu,0x000Fu};
    long*pc[3]={&pA,&pL,&pC},*tc[3]={&tA,&tL,&tC},*mm[3]={&mA,&mL,&mC};
    for(int k=0;k<3;k++){int pr=supp&&!(supp&~TG[k]),tr=ts&&!(ts&~TG[k]);
      if(pr)(*pc[k])++; if(tr)(*tc[k])++; if(pr!=tr)(*mm[k])++;}
  }
  fprintf(stderr,"stripe %ld: suppMis %ld maxNclasses %ld | ALT t %ld p %ld m %ld | LINE t %ld p %ld m %ld | CC t %ld p %ld m %ld\n",
    stripe,suppMis,maxn,tA,pA,mA,tL,pL,mL,tC,pC,mC);
  return 0;
}
