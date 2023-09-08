#!/usr/bin/env python

# calculate and print expfilter strength

import numpy as np

# exp filter
def expfilter_fac(inds, alp, s):
    N = len(inds)-1
    f = np.exp( -alp * np.power(inds/N, s) )
    return f


def print_arrays(inds, *argv):
    n = len(inds)
    print('# inds array1 array2 ...')
    for i in range(n):
        print(inds[i], end=' ')
        for arg in argv:
            print(arg[i], end=' ')
        print('')

# use 10 points
n = 10
inds = np.array(range(n))

alp = 36
s = 32
f32 = expfilter_fac(inds, alp, s)
#data = np.column_stack((inds, f))
#print(data)
#print_arrays(inds, f)

alp = 36
s = 40
f40 = expfilter_fac(inds, alp, s)

alp = 36
s = 64
f64 = expfilter_fac(inds, alp, s)

# print all 3 results
print_arrays(inds, f32, f40, f64)
print('# for plotting use:')
print('# ./filters.py > f.txt ; tgraph.py -m -c 1:2 f.txt -c 1:3 f.txt -c 1:4 f.txt')
