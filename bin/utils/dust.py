#!/usr/bin/env python
from __future__ import division
import numpy as np


def dustModel(Q, v, aMin, aMax)
    '''
    Qt   : total dust production rate
    v    : dust outflow velocity
    aMin : minimum dust bin size
    aMax : maximum dust bin size
    '''
    na = 100   # number of dust bins
    nr = 2000  # number of spatial bins
    rMin = 2000.0 # comet radius
    rMax = 1e8    # max radius to calculate    

    #Get the dust distribution first    
    aDist, qDust = dustDistribution(Qt, v, aMin, aMax, na)
    
    #Caluclate the dust density
    lrmin = np.log10(rMin)
    lrmax = np.log10(rMax)
    dr = (rMax - rMin) / nr
    r = 10**(np.arange(rmin, rmax, dr))

    n = qDust / (4 * np.pi * r**2 * v)
    return r, aDist, n
    

def dustDistribution(Qt, v, aMin, aMax, na)
    
    p  = 4.0  #size distribution slope (index)
    pm = p-1
    dustNorm = pm*(aMin*aMax)^pm)/(aMax^pm - aMin^pm)

    laMin = np.log10(aMin)
    laMax = np.log10(aMax)
    da = (laMax - laMin) / na
    aDust = 10**(np.arange(laMin, laMax, da))
    
    qDust = np.arange(na)
    qDust = Qt*dustNorm*aDust^(-p)
    return aDust, Qdust
    

