#!/usr/bin/env python
from __future__ import division
import numpy as np
import scipy.constants as sc

def fluxDensity(columnDensity, radius, iFOV, args):
    # returns: 
    # Stotal which is aveBright integrated over dust radius (a)
    # average brightness as a function of pixel x, pixel y, a, freq
    # and the array wavelenght array

    rAU = 1.3   # AU
    rH  = rAU * sc.au # m
    
    # F = frequency in Hz
    nF = 1
    minF = 600e9
    maxF = 600e9 
    F = np.linspace(minF, maxF, num=nF)

    # get dimensions and create array    
    (nX, nY, nR) = columnDensity.shape
    result = np.zeros((nX, nY, nF)) 

    aveBright = aveBrightness(columnDensity, radius, F, args)
    
    sA = 2.0*np.pi*(1.0-np.cos(iFOV))

    sTotal = sA * np.trapz(aveBright, x=radius, axis=2) 

    # returns: 
    # Stotal which is aveBright integrated over dust radius (a)
    # average brightness as a function of pixel x, pixel y, a, freq
    # and the array wavelenght array
    return sTotal, aveBright, F


def aveBrightness(columnDensity, radius, F, args):
    # returns the average brightness as a function of pixel x, pixel y, a, freq

    (nX, nY, nR) = columnDensity.shape
    nF = F.shape[0]
    result = np.zeros((nX, nY, nR, nF)) 

    # Calculate relavant parameters
    B = plank(radius, F)
    Q = qabs(radius, F)
    
    for i in range(nX):
        for j in range(nY):
            for k in range(nR):
                for l in range(nF):
                    result[i,j,k,l] = np.pi * radius[k]**2 * B[k,l] * Q[k,l] * columnDensity[i,j,k]

    # returns the average brightness as a function of pixel x, pixel y, a, freq
    return result

def qabs(radius, F):
    nF = F.shape[0]
    nR = radius.shape[0]
    result = np.zeros((nR, nF))

    # For testing use olivine = 0.8687
    result = result + 0.8687

    return result

def plank(radius, F):
    nF = F.shape[0]
    nR = radius.shape[0]
    result = np.zeros((nR, nF))
    
    T = dustT(radius)              

    #black body radiation: 2 h freq^3 / c^2 / (exp(h freq / kT)-1)
    for i in range(nR):
        for j in range(nF):
            result[i,j] = 2.0 * sc.h * F[j]**3 / sc.c**2 / (np.exp(sc.h * F[j] / sc.k / T[i])-1)

    return result

def dustT(radius):
    result = np.zeros_like(radius)
    # For testing use an equilibirum temp 244K
    result = result + 244.0

    return result
   



