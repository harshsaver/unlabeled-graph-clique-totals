#!/bin/bash
# Independent brute force for n=11: all 1,018,997,864 graphs via nauty geng, split 10 ways
cd "$(dirname "$0")"
seq 0 9 | xargs -P 5 -I{} sh -c 'geng -q 11 {}/10 | ./brute2 > brute11/part{}.txt'
python3 -c "
import glob
t=None
for f in glob.glob('brute11/part*.txt'):
    v=list(map(int,open(f).read().split()))
    t=v if t is None else [a+b for a,b in zip(t,v)]
print(11,sum(t),*t)" > brute11/total.txt
