#!/usr/bin/env python

import numpy as np

f = open('/Users/abieler/pyLOS/testComa.txt', 'w')
for r in np.arange(np.log10(2000), np.log10(1e8), np.log10(1e8) / 1000):
    f.write("%e, %e\n" % (10**r, 1 / 10**r * 1000))
f.close()