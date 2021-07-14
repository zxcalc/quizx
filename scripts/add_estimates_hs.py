#!/bin/env python3

import csv
import sys

def terms(t):
    count = 7.0 ** (t // 6)
    t = t % 6
    count *= 2.0 ** (t // 2)
    t = t % 2
    if t == 1: count *= 2.0
    return count

with open(sys.argv[1]) as f:
    r = csv.reader(f)
    data = list(r)
    # data = data[1:]
    print('"q","ccz","seed","terms","time","tcounts","unsimp_terms","naive_terms"')
    for d in data:
        d.append(str(int(sum(terms(int(t)) for t in d[5].split(',')))))
        d.append(str(terms(int(d[1]) * 4 * 7)))
        print('"' + '","'.join(d) + '"')
