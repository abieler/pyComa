#!/usr/bin/env python
from __future__ import division
import numpy as np


def dustModel(Qt, v)
    '''
    Qt   : total dust production rate
    v    : dust outflow velocity
    aMin : minimum dust bin size
    aMax : maximum dust bin size
    '''
    na = 100   # number of dust bins
    nr = 2000  # number of spatial bins
    aMin = 1e-9 # minimum dust bin size - 1 nm
    aMax = 1.0  # maximum dust bin size - 1 m
    rMin = 2000.0 # comet radius
    rMax = 1e8    # max radius to calculate    

    #Get the initial dust distribution at the surface    
    aDust, qDust = dustDistribution(Qt, v, aMin, aMax, na)
    
    #Calculate the dust density
    lrmin = np.log10(rMin)
    lrmax = np.log10(rMax)
    dr = (rMax - rMin) / nr
    r = 10**(np.arange(rmin, rmax, dr))

    n = qDust / (4 * np.pi * r**2 * v)  # uniform radial expansion model

    return r, n, aDist
    

def dustDistribution(Qt, aMin, aMax, na)
    
    laMin = np.log10(aMin)
    laMax = np.log10(aMax)
    da = (laMax - laMin) / na
    aDust = 10**(np.arange(laMin, laMax, da))

    dustDensity = 1000.0 (kg/m^3)   # for now use water density

    dustNorm = 3.0/4.0/np.pi/dustDensity/(log(aMax)-log(aMin))
    
    qDust = np.arange(na)
    qDust = Qt*dustNorm*aDust^(-4)
    return aDust, Qdust
    


