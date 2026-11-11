/* brute.c: read graph6 lines (n<=16) from stdin, output total k-cliques per k.
   Independent check: geng -q N | ./brute */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
int main(void){
  char line[256]; unsigned long long tot[17]={0}; int n=-1;
  while(fgets(line,sizeof line,stdin)){
    n=line[0]-63; uint32_t adj[16]={0}; int bit=0;
    for(int j=1;j<n;j++) for(int i=0;i<j;i++){ int c=line[1+bit/6]-63; if(c>>(5-bit%6)&1){adj[i]|=1u<<j;adj[j]|=1u<<i;} bit++; }
    for(uint32_t S=0;S<(1u<<n);S++){ int ok=1; for(int v=0;v<n&&ok;v++) if(S>>v&1) if((S&~adj[v]&~(1u<<v))) ok=0; if(ok) tot[__builtin_popcount(S)]++; }
  }
  unsigned long long s=0; for(int k=0;k<=n;k++) s+=tot[k];
  printf("%d %llu",n,s); for(int k=0;k<=n;k++) printf(" %llu",tot[k]); printf("\n"); return 0; }
