#!/usr/bin/env python
from __future__ import division
import numpy as np

path = '/Users/abieler/pyComa/results'

f = open(path + '/testUserFile_1d.txt', 'w')
for r in np.arange(np.log10(2000), np.log10(1e8), np.log10(1e8) / 1000):
    f.write("%e, %e\n" % (10**r, 1 / 10**r * 1000))
f.close()



f = open(path + '/testUserFile_2d.txt', 'w')
for r in np.arange(np.log10(2000), np.log10(1e8), np.log10(1e8) / 201):
    rr = 10**r
    print '%.2e' % rr
    for phi in np.arange(0, 181, 1):
        x = rr * np.cos(phi/180*np.pi)
        y = rr * np.sin(phi/180*np.pi)
        n = 1e10 * (1/rr)**3 * np.exp(-phi/10)
        f.write("%e,%e,%e\n" % (x, y, n))
f.close()