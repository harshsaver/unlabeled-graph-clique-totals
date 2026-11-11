# OEIS submission text (paste into the edit forms on oeis.org)

Each entry gets: a b-file (upload the file in this repo), the formula, the program, and a link.
Data lines may hold at most ~260 characters, so only the first terms are shown there; the b-file has the rest.

---

## A401214 (triangle T(n,k))

**b-file:** `b401214.txt` (rows n = 0..34)

**Formula:**
T(n,k) = Sum_{partitions p of n} (1/z_p) * Sum_{k-subsets S of the vertices} 2^(E(p) - h(p,S)), where for any permutation g of cycle type p, E(p) = Sum_i floor(p_i/2) + Sum_{i<j} gcd(p_i,p_j) is the number of edge orbits of g, h(p,S) is the number of those orbits containing a pair of vertices of S, and z_p = Prod_i i^m_i * m_i! (weighted Burnside lemma).

**Link:** _[Your Name]_, <a href="https://github.com/harshsaver/unlabeled-graph-clique-totals">Weighted Burnside computation of total cliques</a>, GitHub (2026).

**Program:** (Python) — see below.

**Extensions:** Rows n = 10..34 from _[your OEIS name]_, Nov 2026

---

## A401215 (total cliques)

**b-file:** `b401215.txt` (n = 0..34)

**Data (extend to):**
```
1,2,7,23,97,442,2833,25614,397194,11379248,632101893,67443396800,13594172910933,5123168574712765,3597930283198608997,4710547737087550722120,11521447475590437425190429,52793025789398255403767485402,454564072638403065266397219769159
```

**Formula:** a(n) = Sum_{k=0..n} T(n,k), with T(n,k) as in A401214.

**Extensions:** a(11)-a(34) from _[your OEIS name]_, Nov 2026

---

## A401216 (total triangles)

**b-file:** `b401216.txt` (n = 1..60)

**Data (extend to):**
```
0,0,1,8,57,468,5126,92020,2976692,182783912,21163085164,4554841114104,1808273858763608,1322971478637440480,1788051469200476074624,4480902768602048580542784,20906477656440674659968086912,182341312978022657239502820397184,2984323129772262227223174393095831712
```

**Formula:** a(n) = T(n,3) with T as in A401214; asymptotically a(n) ~ A000088(n)*binomial(n,3)/8.

**Extensions:** a(11)-a(60) from _[your OEIS name]_, Nov 2026

---

## Program (all three entries)

```python
# Python reference implementation for A401214 (T(n,k)), A401215, A401216.
# Weighted Burnside: T(n,k) = Sum_{partitions p of n} (1/z_p) Sum_{|S|=k} 2^(E(p)-h(p,S)).
from math import gcd, factorial, prod
from itertools import combinations
from collections import Counter
from sympy.utilities.iterables import partitions

def T(n, k):
    tot = 0
    for q in partitions(n):
        c = [l for l, m in q.items() for _ in range(m)]
        cyc = [(i, j) for i, l in enumerate(c) for j in range(l)]  # vertex = (cycle, position)
        E = sum(l//2 for l in c) + sum(gcd(c[i], c[j]) for i in range(len(c)) for j in range(i))
        def orbit(u, v):
            (a, x), (b, y) = sorted([u, v])
            if a == b:
                d = (y - x) % c[a]; return (a, a, min(d, c[a] - d))
            return (a, b, (y - x) % gcd(c[a], c[b]))
        s = sum(2**(E - len({orbit(u, v) for u, v in combinations(S, 2)}))
                for S in combinations(cyc, k))
        tot += s * factorial(n) // prod(l**m * factorial(m) for l, m in q.items())
    return tot // factorial(n)

if __name__ == '__main__':
    print([[T(n, k) for k in range(n+1)] for n in range(8)])       # A401214
    print([sum(T(n, k) for k in range(n+1)) for n in range(10)])    # A401215
    print([T(n, 3) for n in range(1, 13)])                          # A401216
```
