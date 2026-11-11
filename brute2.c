/* brute2.c: fast total k-clique counter over graph6 stream (n<=31), by
   recursive clique extension (each clique visited exactly once). */
#include <stdio.h>
#include <stdint.h>
static uint32_t adj[32]; static unsigned long long tot[33];
static void ext(uint32_t cand,int k){ tot[k]++; while(cand){int v=__builtin_ctz(cand); cand&=cand-1; ext(cand&adj[v],k+1);} }
int main(void){
  char line[256]; int n=0;
  while(fgets(line,sizeof line,stdin)){
    n=line[0]-63; for(int i=0;i<n;i++) adj[i]=0; int bit=0;
    for(int j=1;j<n;j++) for(int i=0;i<j;i++){ int c=line[1+bit/6]-63; if(c>>(5-bit%6)&1){adj[i]|=1u<<j;adj[j]|=1u<<i;} bit++; }
    /* cand restricted to higher-index vertices so each clique counted once */
    for(int i=0;i<n;i++) adj[i]&= ~((2u<<i)-1);
    ext((n==32)?~0u:((1u<<n)-1),0);
  }
  for(int k=0;k<=n;k++) printf("%llu%c",tot[k],k==n?'\n':' '); return 0; }
