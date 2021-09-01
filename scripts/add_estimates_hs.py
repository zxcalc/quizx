#!/bin/env python3

import csv
import sys
from bigfloat import BigFloat, log10

def terms(t):
    count = BigFloat('7.0') ** (t // 6)
    t = t % 6
    count *= BigFloat('2.0') ** (t // 2)
    t = t % 2
    if t == 1: count *= 2.0
    return count

with open(sys.argv[1]) as f:
    r = csv.reader(f)
    data = list(r)
    # data = data[1:]
    print('"q","ccz","seed","terms","time","tcounts","unsimp_terms","naive_terms","log10_unsimp_terms","log10_naive_terms"')
    for d in data:
        unsimp_terms = int(sum(terms(int(t)) for t in d[5].split(',')))
        naive_terms = terms(int(d[1]) * 4 * 7)
        d.append(str(unsimp_terms))
        d.append(str(naive_terms))
        d.append(str(log10(unsimp_terms)))
        d.append(str(log10(naive_terms)))
        print('"' + '","'.join(d) + '"')
