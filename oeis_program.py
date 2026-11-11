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
