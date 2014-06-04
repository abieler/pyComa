#!/usr/bin/env python
from __future__ import division
import numpy as np

dustDensity = 1000.0 # (kg/m^3)   # for now use water density

def dustModel(Qt, v):
    '''
    Qt   : total dust production rate
    v    : dust outflow velocity
    aMin : minimum dust bin size
    aMax : maximum dust bin size
    '''
    na = 10   # number of dust bins
    nr = 20   # number of spatial bins
    aMin = 1e-9 # minimum dust bin size - 1 nm
    aMax = 1.0  # maximum dust bin size - 1 m
    rMin = 2000.0 # comet radius
    rMax = 1e9    # max radius to calculate    

    #Get the initial dust distribution at the surface    
    aDust, qDust = dustDistribution(Qt, aMin, aMax, na)
    na = aDust.shape[0]
    
    # Create the spatial array
    lrMin = np.log10(rMin)
    lrMax = np.log10(rMax)
    r = 10**(np.linspace(lrMin, lrMax, num=nr))

    #Calculate the dust density
    nDensity = np.zeros((nr,na))
    mDensity = np.zeros((nr,na))

    for i in range(nr):
        for j in range(na):
            # simple radial expansion model
            nDensity[i,j] = qDust[j] / (4.0 * np.pi * r[i]**2 * v)  
            mDensity[i,j] = 4.0/3.0*np.pi*aDust[j]**3*nDensity[i,j]
            # Constant model just for testing
            # nDensity[i,j] = 1.0
            # mDensity[i,j] = 1.0
            # 1/r^2 model just for testing
            nDensity[i,j] = 1.0/r[i]**2
            mDensity[i,j] = 1.0/r[i]**2

    return r, nDensity, mDensity, aDust
    

def dustDistribution(Qt, aMin, aMax, na):
    
    laMin = np.log10(aMin)
    laMax = np.log10(aMax)
    aDust = 10**(np.linspace(laMin, laMax, num=na))

    dustNorm = 3.0/4.0/np.pi/dustDensity/(np.log(aMax)-np.log(aMin))
    # qDust = Qt*dustNorm*aDust**(-4)  # Standard model
    qDust = np.zeros(na)+1.0          # Constant value for testing

    return aDust, qDust
    
