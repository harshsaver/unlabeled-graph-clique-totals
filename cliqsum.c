/* cliqsum.c -- weighted Burnside enumeration of total k-cliques over all
 * unlabeled simple graphs on n vertices.
 *
 * For a permutation g of cycle type lambda, the graphs fixed by g are unions of
 * edge-orbits of g (E(lambda) of them).  A vertex set S is a clique in such a
 * graph iff every edge-orbit meeting a pair inside S is present, so
 *     #{G fixed by g : S clique in G} = 2^(E(lambda) - h(lambda,S)),
 * where h = number of edge orbits meeting pairs of S.  By weighted Burnside
 * (sum over orbits of an invariant f = (1/|Gamma|) sum_g sum_{x in Fix g} f(x)):
 *     T(n,k) = sum_lambda (1/z_lambda) sum_{|S|=k} 2^(E - h(S)).
 * This program, for each partition lambda of n, outputs the exact table
 * cnt[k][h] = #{S : |S|=k, h(S)=h}.  Big-integer aggregation is in aggregate.py.
 *
 * Edge orbits: inside a c-cycle, pairs at cyclic distance delta (1<=delta<=c/2)
 * form one orbit each.  Between an a-cycle and a b-cycle, pair (i,j) lies in
 * orbit (j-i) mod gcd(a,b).  So h(S) = sum_cycles #{distinct distances in S_C}
 *   + sum_{A<B} |R_A - R_B|  (difference set in Z_gcd, R = residues of S mod gcd).
 * h is invariant under rotating S inside each cycle independently, so each
 * cycle's part of S is enumerated up to rotation (necklaces weighted by orbit
 * size); equal-length cycles are enumerated as multisets; fixed points are
 * handled in closed form; the (unique) cycle longer than n/2, if any, is
 * streamed rather than stored.
 *
 * Usage: cliqsum n [kmax]   Output: "P c1 c2 ..." then "k h cnt" lines then "END".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 62
typedef unsigned long long u64;

static int n, kmax;
static int gcdt[MAXN+1][MAXN+1];
static u64 binom[MAXN+1][MAXN+1];

typedef struct { uint64_t mask; int orb, pop, internal; } Neck;
static Neck *necks[MAXN+1]; static long nneck[MAXN+1];

static inline uint64_t fullmask(int c){ return c>=64?~0ULL:((1ULL<<c)-1); }
static inline uint64_t rotl(uint64_t m, int t, int c){ t%=c; if(!t) return m; return ((m<<t)|(m>>(c-t)))&fullmask(c); }
static inline uint64_t fold(uint64_t m, int c, int d){ uint64_t r=0, f=fullmask(d); for(int i=0;i<c;i+=d) r|=(m>>i)&f; return r; }
static inline uint64_t negate(uint64_t r, int d){ uint64_t nr=r&1; for(int i=1;i<d;i++) if(r>>i&1) nr|=1ULL<<(d-i); return nr; }
static inline int diffcount(uint64_t A, uint64_t negB, int d){
    uint64_t acc=0, f=fullmask(d);
    while(A){ int x=__builtin_ctzll(A); A&=A-1; acc|=rotl(negB,x,d); if(acc==f) break; }
    return __builtin_popcountll(acc);
}

/* iterate canonical necklaces of length c with popcount <= kmax */
typedef void (*neck_cb)(const Neck*, void*);
static void iterate_necks(int c, neck_cb cb, void *ctx){
    int pmax = kmax<c?kmax:c;
    uint64_t lim = fullmask(c);
    for(int pop=0; pop<=pmax; pop++){
        uint64_t m = pop==0?0:fullmask(pop);
        for(;;){
            int ismin=1, orb=c;
            for(int t=1;t<c;t++){ uint64_t r=rotl(m,t,c); if(r<m){ismin=0;break;} if(r==m && orb==c) orb=t; }
            if(ismin){
                Neck N; N.mask=m; N.orb=orb; N.pop=pop;
                uint64_t seen=0;
                for(int i=0;i<c;i++) if(m>>i&1) for(int j=i+1;j<c;j++) if(m>>j&1){ int d=j-i; if(c-d<d) d=c-d; seen|=1ULL<<d; }
                N.internal=__builtin_popcountll(seen);
                cb(&N,ctx);
            }
            if(pop==0 || m==lim) break;
            uint64_t u=m&-m, v=m+u; if(v==0) break; m=v+(((v^m)/u)>>2); if(m>lim) break;
        }
    }
}
typedef struct { Neck *a; long cnt, cap; } Store;
static void store_cb(const Neck *N, void *ctx){ Store *s=ctx; if(s->cnt==s->cap){ s->cap=s->cap?2*s->cap:64; s->a=realloc(s->a,sizeof(Neck)*s->cap);} s->a[s->cnt++]=*N; }

static int cyc[MAXN], ncyc, nfix;
static const Neck *chosen[MAXN];
static long choice[MAXN];
static u64 *table; static int HMAX;

static void leaf(u64 w, int k0, int h0, int q){
    for(int m=0;m<=nfix && k0+m<=kmax;m++){
        int h=h0+m*(m-1)/2+m*q;
        table[(size_t)(k0+m)*(HMAX+1)+h] += w*binom[nfix][m];
    }
}
static int cross(const Neck *N, int c, int i){
    int h=0;
    if(!N->pop) return 0;
    for(int j=0;j<i;j++){ const Neck *M=chosen[j]; if(M->pop){ int d=gcdt[c][cyc[j]];
        h+=diffcount(fold(N->mask,c,d), negate(fold(M->mask,cyc[j],d),d), d);} }
    return h;
}
static u64 run_multinomial(int runstart, int i){
    u64 mn=1; int pos=0, a=runstart;
    while(a<=i){ int b=a; while(b+1<=i && choice[b+1]==choice[a]) b++; int len=b-a+1; mn*=binom[pos+len][len]; pos+=len; a=b+1; }
    return mn;
}
static void dfs(int i, u64 w, int k0, int h0, int q, int runstart);
typedef struct { u64 w; int k0,h0,q; } StreamCtx;
static void stream_cb(const Neck *N, void *vctx){
    StreamCtx *s=vctx; int i=ncyc-1, c=cyc[i];
    if(s->k0+N->pop>kmax) return;
    chosen[i]=N;
    int h=s->h0+N->internal+cross(N,c,i);
    leaf(s->w*(u64)N->orb, s->k0+N->pop, h, s->q+(N->pop>0));
}
static void dfs(int i, u64 w, int k0, int h0, int q, int runstart){
    if(i==ncyc){ leaf(w,k0,h0,q); return; }
    int c=cyc[i];
    if(2*c>n){ /* unique, last, streamed */
        StreamCtx s={w,k0,h0,q}; iterate_necks(c,stream_cb,&s); return;
    }
    int newrun = (i==0 || cyc[i-1]!=c);
    if(newrun) runstart=i;
    long start = newrun?0:choice[i-1];
    int runend = (i==ncyc-1 || cyc[i+1]!=c);
    for(long t=start;t<nneck[c];t++){
        const Neck *N=&necks[c][t];
        if(k0+N->pop>kmax) continue;
        chosen[i]=N; choice[i]=t;
        int h=h0+N->internal+cross(N,c,i);
        u64 w2=w*(u64)N->orb;
        if(runend) w2*=run_multinomial(runstart,i);
        dfs(i+1,w2,k0+N->pop,h,q+(N->pop>0),runstart);
    }
}

static int part[MAXN+1];
static void do_partition(int len){
    nfix=0; ncyc=0;
    for(int i=len-1;i>=0;i--){ if(part[i]==1) nfix++; else cyc[ncyc++]=part[i]; }
    memset(table,0,sizeof(u64)*(size_t)(kmax+1)*(HMAX+1));
    dfs(0,1,0,0,0,0);
    printf("P"); for(int i=0;i<len;i++) printf(" %d",part[i]); printf("\n");
    for(int k=0;k<=kmax;k++) for(int h=0;h<=HMAX;h++){ u64 v=table[(size_t)k*(HMAX+1)+h]; if(v) printf("%d %d %llu\n",k,h,v); }
    printf("END\n");
}
static void gen(int rem, int maxp, int len){
    if(rem==0){ do_partition(len); return; }
    for(int p=(rem<maxp?rem:maxp);p>=1;p--){ part[len]=p; gen(rem-p,p,len+1); }
}

int main(int argc,char**argv){
    n=atoi(argv[1]); kmax=argc>2?atoi(argv[2]):n; if(kmax>n) kmax=n;
    if(n>MAXN){fprintf(stderr,"n too big\n");return 1;}
    if(kmax==n && n>60){fprintf(stderr,"counts could overflow\n");return 1;}
    for(int a=0;a<=MAXN;a++) for(int b=0;b<=MAXN;b++){ int x=a,y=b; while(y){int t=x%y;x=y;y=t;} gcdt[a][b]=x; }
    for(int a=0;a<=MAXN;a++){ binom[a][0]=1; for(int b=1;b<=a;b++) binom[a][b]=binom[a-1][b-1]+(b<=a-1?binom[a-1][b]:0); }
    for(int c=2;2*c<=n;c++){ Store s={0,0,0}; iterate_necks(c,store_cb,&s); necks[c]=s.a; nneck[c]=s.cnt; }
    HMAX=n*(n-1)/2;
    table=malloc(sizeof(u64)*(size_t)(kmax+1)*(HMAX+1));
    if(n==0){ printf("P\n0 0 1\nEND\n"); return 0; }
    gen(n,n,0);
    return 0;
}
