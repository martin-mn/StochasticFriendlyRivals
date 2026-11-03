/* effcc.c -- the closed form  Eff_CC = { p0 = 1 and K >= E + 2 }  checked, as a SET EQUALITY,
 * against the regenerated list of W_CC = 1 support strata. */
#include <stdio.h>
#include <stdlib.h>
static int MIR[16]; static inline int swapo(int o){return ((o&1)<<1)|(o>>1);}
#define INF 1000
int main(int argc,char**argv){
  for(int j=0;j<16;j++) MIR[j]=4*swapo(j>>2)+swapo(j&3);
  long TOT=1; for(int i=0;i<16;i++) TOT*=3;
  unsigned char*truth=calloc((TOT+7)/8,1);
  FILE*f=fopen(argv[1],"r"); long g,nt=0; while(fscanf(f,"%ld",&g)==1){truth[g>>3]|=1u<<(g&7);nt++;} fclose(f);
  long stripe=atol(argv[2]),nstr=atol(argv[3]);
  long npred=0,fp=0,fn=0; int dg[16],C[16][16],D[16][16]; unsigned A0[16];
  for(g=stripe;g<TOT;g+=nstr){
    long x=g; for(int j=0;j<16;j++){dg[j]=x%3;x/=3;}
    int istruth=(truth[g>>3]>>(g&7))&1, pred=0;
    if(dg[0]==2){
      for(int i=0;i<16;i++)for(int j=0;j<16;j++)C[i][j]=INF;
      for(int j=0;j<16;j++){unsigned m0=0;int dj=dg[j],dm=dg[MIR[j]];
        for(int aF=0;aF<2;aF++){int cF=(aF==0)?(dj==0):(dj==2);
          for(int aA=0;aA<2;aA++){int cA=(aA==0)?(dm==0):(dm==2);
            int nx=4*(2*aF+aA)+(j>>2),c=cF+cA; if(c<C[j][nx])C[j][nx]=c; if(c==0)m0|=1u<<nx;}}
        A0[j]=m0;}
      for(int i=0;i<16;i++)for(int j=0;j<16;j++)D[i][j]=(i==j)?0:C[i][j];
      for(int k=0;k<16;k++)for(int i=0;i<16;i++){int d=D[i][k];if(d>=INF)continue;
        for(int j=0;j<16;j++){int v=d+D[k][j];if(v<D[i][j])D[i][j]=v;}}
      unsigned F[16];
      for(int s=0;s<16;s++){unsigned S=1u<<s,q=0;while(S!=q){q=S;for(int j=0;j<16;j++)if(S>>j&1)S|=A0[j];}F[s]=S;}
      unsigned Z[16]; int m=0,i0=-1,done=0;
      for(int s=0;s<16;s++){ if(done>>s&1)continue;
        int rec=1; for(int t=0;t<16;t++) if((F[s]>>t&1)&&!((F[t]>>s)&1)){rec=0;break;}
        if(!rec)continue; Z[m]=F[s];done|=F[s];if(s==0)i0=m;m++; }
      if(m==1) pred=1;                                   /* K = infinity, E = 0 */
      else{
        int V[16][16];
        for(int a=0;a<m;a++)for(int b=0;b<m;b++){ if(a==b){V[a][b]=0;continue;} int best=INF;
          for(int u=0;u<16;u++)if(Z[a]>>u&1)for(int v=0;v<16;v++)if(Z[b]>>v&1) if(D[u][v]<best)best=D[u][v];
          V[a][b]=best;}
        int K=INF; for(int a=0;a<m;a++) if(a!=i0 && V[i0][a]<K) K=V[i0][a];
        int idx[16],n=0; for(int a=0;a<m;a++) if(a!=i0) idx[n++]=a;
        static int fd[1<<15]; int NS2=1<<n;
        for(int S=0;S<NS2;S++) fd[S]=INF*32; fd[0]=0;
        for(int S=0;S<NS2;S++){ if(fd[S]>=INF*32)continue;
          for(int t=0;t<n;t++){ if(S>>t&1)continue; int best=V[idx[t]][i0];
            for(int u=0;u<n;u++) if(S>>u&1){int v=V[idx[t]][idx[u]]; if(v<best)best=v;}
            int nv=fd[S]+best; if(nv<fd[S|(1<<t)])fd[S|(1<<t)]=nv; } }
        int E=fd[NS2-1]-(m-1);
        pred = (K>=E+2);
      }
    }
    if(pred)npred++;
    if(pred&&!istruth)fp++; if(!pred&&istruth)fn++;
  }
  fprintf(stderr,"stripe %ld: closed-form %ld  false-pos %ld  false-neg %ld  (truth total %ld)\n",stripe,npred,fp,fn,nt);
  return 0;
}
