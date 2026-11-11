#!/usr/bin/env python3
"""Run ./cliqsum n kmax and aggregate exactly:
T(n,k) = (1/n!) * sum_lambda (n!/z_lambda) * sum_h cnt[k][h] * 2^(E(lambda)-h).
Prints row n of the triangle T(n,0..kmax)."""
import subprocess, sys
from math import gcd, factorial
from collections import Counter

def E(parts):
    e = sum(c // 2 for c in parts)
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            e += gcd(parts[i], parts[j])
    return e

def class_size(n, parts):
    z = 1
    for c, m in Counter(parts).items():
        z *= c ** m * factorial(m)
    return factorial(n) // z

def row(n, kmax=None, exe='./cliqsum'):
    kmax = n if kmax is None else min(kmax, n)
    p = subprocess.Popen([exe, str(n), str(kmax)], stdout=subprocess.PIPE, text=True)
    tot = [0] * (kmax + 1)
    parts = None; e = 0; cs = 0; acc = None
    for line in p.stdout:
        if line.startswith('P'):
            parts = list(map(int, line.split()[1:])); e = E(parts); cs = class_size(n, parts)
            acc = [0] * (kmax + 1)
        elif line.startswith('END'):
            for k in range(kmax + 1):
                tot[k] += cs * acc[k]
        else:
            k, h, c = map(int, line.split())
            acc[k] += c << (e - h)
    p.wait()
    assert p.returncode == 0
    f = factorial(n)
    out = []
    for t in tot:
        assert t % f == 0, "non-integer result: bug"
        out.append(t // f)
    return out

if __name__ == '__main__':
    n = int(sys.argv[1]); kmax = int(sys.argv[2]) if len(sys.argv) > 2 else None
    r = row(n, kmax)
    print(n, sum(r), *r)
