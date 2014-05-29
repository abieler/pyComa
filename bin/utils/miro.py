#!/usr/bin/env python
from __future__ import division
import numpy as np
import spice
import sqlite3

def calculateBrightness(N_oversampleX, N_oversampleY, ccd, args):
    
    ccdFinal = np.zeros(19)
    result = []

    ll = 0
    for k in range(7):              # loop over first 7 pixels
        ccdFinal[k] = np.mean(ccd[k * N_oversampleX:(k + 1) * N_oversampleX, :])

    for k in range(14, 19):          # loop over last 5 pixels
        ccdFinal[k] = np.mean(ccd[k * N_oversampleX:(k + 1) * N_oversampleX, :])

    ll = int(N_oversampleY * 0.3 // 2)
    for k in range(7, 8):
        ccdFinal[k] = np.mean(ccd[k * N_oversampleX:(k + 1) * N_oversampleX, ll:-ll])

    ll = int(N_oversampleY * 0.2 // 2)
    for k in range(13, 14):
        ccdFinal[k] = np.mean(ccd[k * N_oversampleX:(k + 1) * N_oversampleX, ll:-ll])

    ll = int(N_oversampleY * 0.5 // 2)
    for k in range(8, 13):
        ccdFinal[k] = np.mean(ccd[k * N_oversampleX:(k + 1) * N_oversampleX, ll:-ll])

    # calculate brightness if alice_spectra was selected
    if args.iInstrumentSelector == 7:
        gFactors, wavelengths = get_gfactor_from_db(args)
        for gFactor in gFactors:
            ccdF = ccdFinal * gFactor / (4 * np.pi) * pixelFOV
            result.append(ccdF)

        result.append(ccdFinal)
    else:
        wavelengths = []
        result.append(ccdFinal)

    return np.array(result), wavelengths


def calculate_column(nRay, dTravel, pixelSize=10**-6,
                     iFOV=9.39*10**-6, gFactor=2*10**-7):

    #N = np.trapz(iFOV * np.array(dTravel)**2 * nRay, dTravel)
    N = np.trapz(nRay, dTravel)
    return N


def qabs(a,freq,nc):



