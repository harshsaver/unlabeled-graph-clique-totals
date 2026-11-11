#!/usr/bin/env python3
"""Build OEIS-style b-files from rows/ and tri/ outputs."""
import os

def rd(p):
    return list(map(int, open(p).read().split())) if os.path.exists(p) and os.path.getsize(p) > 0 else None

rows = {}
for n in range(0, 63):
    r = rd(f'rows/row{n}.txt')
    if r: rows[n] = r[2:]          # T(n,0..n)
N = max(n for n in rows if all(m in rows for m in range(n + 1)))

with open('b401215.txt', 'w') as f:   # row sums, offset 0
    for n in range(N + 1):
        f.write(f'{n} {sum(rows[n])}\n')

with open('b401214.txt', 'w') as f:   # triangle flattened, offset 0
    i = 0
    for n in range(N + 1):
        for v in rows[n]:
            f.write(f'{i} {v}\n'); i += 1

with open('b401216.txt', 'w') as f:   # column k=3, offset 1
    for n in range(1, 61):
        if n in rows:
            v = rows[n][3] if n >= 3 else 0
        else:
            v = rd(f'tri/tri{n}.txt')[2 + 3]
        f.write(f'{n} {v}\n')
print('rows complete to n =', N)
