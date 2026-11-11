# Total cliques over all unlabeled graphs — new terms for OEIS A401214 / A401215 / A401216

**Result.** New exact values of three OEIS sequences (created Oct 2026, tagged `more`;
as of 2026-11-11 the OEIS lists only up to n = 9 / 10):

| Sequence | Meaning | OEIS had | Now computed |
|---|---|---|---|
| A401214 | T(n,k) = total number of k-cliques over all unlabeled simple graphs on n vertices (triangle) | rows n ≤ 9 | rows n ≤ 31 |
| A401215 | total number of cliques (any size, incl. empty) over all unlabeled graphs on n vertices | n ≤ 10 | n ≤ 31 |
| A401216 | total number of triangles over all unlabeled graphs on n vertices | n ≤ 10 | n ≤ 60 |

(Equivalently: total independent sets, by complementation.)  Full data: `b401214.txt`,
`b401215.txt`, `b401216.txt` (OEIS b-file format), raw rows in `rows/`, `tri/`.

First new terms:

```
A401215: a(11) = 67443396800, a(12) = 13594172910933, a(13) = 5123168574712765,
         a(14) = 3597930283198608997, a(15) = 4710547737087550722120, ...
A401216: a(11) = 21163085164, a(12) = 4554841114104, a(13) = 1808273858763608,
         a(14) = 1322971478637440480, a(15) = 1788051469200476074624, ...
A401214 row 10: 12005168, 120051680, 270116280, 182783912, 42903490, 4039723, 194133, 7160, 326, 20, 1
A401214 row 11: 1018997864, 11208976504, 28022441260, 21163085164, 5471208332, 535219963,
                22894354, 560511, 12409, 416, 22, 1
```

## Method (new formula)

Weighted Burnside: for a group Γ acting on X and a Γ-invariant f,
Σ_{orbits O} f(O) = (1/|Γ|) Σ_g Σ_{x ∈ Fix(g)} f(x).
Take X = labeled graphs on [n], Γ = S_n, f(G) = #k-cliques of G. A graph fixed by g is a
union of edge-orbits of g; it contains S as a clique iff it contains all h(g,S) edge-orbits
that meet pairs inside S. Hence

  **T(n,k) = Σ_{λ ⊢ n} (1/z_λ) Σ_{|S|=k} 2^{E(λ) − h(λ,S)}**,

where E(λ) = Σ⌊c_i/2⌋ + Σ_{i<j} gcd(c_i,c_j) is the number of edge orbits, and
h(λ,S) = Σ_cycles #{distinct cyclic distances inside S∩C} + Σ_{A<B} |R_A − R_B|, where R_A is
the set of residues of S∩A mod gcd(|A|,|B|), with the difference set taken in Z_gcd.

`cliqsum.c` evaluates the inner sum exactly for each cycle type. It enumerates S∩C up to
rotation (h does not change when each cycle is rotated separately), groups equal-length cycles
as multisets, and handles fixed points in closed form. `aggregate.py` does the big-integer
sum. Row 30 takes a few minutes on a laptop; A401216 alone (k ≤ 3) reaches n = 60 in about a minute.

## Verification

1. **Reproduces every known OEIS term** (A401214 rows 0–9, A401215 n ≤ 10, A401216 n ≤ 10),
   and column k=0 reproduces A000088 (number of graphs), e.g. g(16) = 64001015704527557894928.
2. **Independent brute force** (`brute2.c` + nauty `geng`): all 1,018,997,864 graphs on 11
   vertices were generated and every clique counted. The resulting row 11 is identical to the
   Burnside output (`brute11/total.txt`). The same check passes for n ≤ 10.
3. **Integrality**: each Σ_λ (n!/z_λ)·(...) must be divisible by n!. This is asserted for
   every value and never fails, which a bug would almost surely break.
4. **Internal consistency**: two independent runs (full enumeration vs. k ≤ 3 pruned
   enumeration) agree on columns 0–3 for every overlapping n. Also T(n,1) = n·g(n) and
   T(n,2) = C(n,2)·g(n)/2 (complementation) hold exactly.
5. **Asymptotics**: T(n,3) / (g(n)·C(n,3)/8) → 1 (1.00696 at n=11, 1.0000115 at n=20,
   1 + 1e-11 at n=40), as expected since almost all graphs are asymmetric.

## Reproduce

```
clang -O3 -o cliqsum cliqsum.c
python3 aggregate.py 12          # prints: n, row sum (A401215), T(n,0..n)
python3 aggregate.py 40 3        # columns k<=3 only (k=3 is A401216)
geng -q 11 | ./brute2            # independent brute force (needs nauty), ~10 min
```
