#!/bin/bash
# full rows of the triangle T(n,k), n=0..34
cd "$(dirname "$0")"
seq 34 -1 0 | xargs -P 8 -I{} sh -c 'python3 aggregate.py {} > rows/row{}.txt'
