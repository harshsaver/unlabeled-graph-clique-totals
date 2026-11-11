#!/bin/bash
# column k=3 (A401216) and full k<=3 columns, n=11..60
cd "$(dirname "$0")"
seq 60 -1 11 | xargs -P 2 -I{} sh -c 'python3 aggregate.py {} 3 > tri/tri{}.txt'
