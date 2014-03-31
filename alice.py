#!/usr/bin/env python
from __future__ import division
import numpy as np
import spice
import sqlite3

pixelFOV = np.array([9.39,
                    9.39,
                    9.39,
                    9.39,
                    9.39,
                    9.39,
                    9.39,
                    6.10,
                    4.69,
                    4.69,
                    4.69,
                    4.69,
                    4.69,
                    7.51,
                    9.39,
                    9.39,
                    9.39,
                    9.39,
                    9.39
                    ]) * 1e-6


def calculateBrightness(N_oversampleX, N_oversampleY, ccd, gFactor):


    print 'ccd mean:', np.mean(ccd)
    ccdFinal = np.zeros(19)
    
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

    ccdFinal = ccdFinal * gFactor / 4 * np.pi * pixelFOV

    return ccdFinal


def calculate_column(nRay, dTravel, pixelSize=10**-6,
                    iFOV=9.39*10**-6, gFactor=2*10**-7):

    #N = np.trapz(iFOV * np.array(dTravel)**2 * nRay, dTravel)
    N = np.trapz(nRay, dTravel)
    return N


def get_v_sun(kernelMetaFile, utcStartTime):
    '''
    compute radial component of CG towards/away from sun for given utcTime.
    v_sun: velocity component on the comet-sun-line.
    '''

    spice.furnsh(kernelMetaFile)
    et = spice.str2et(utcStartTime)
    state, lightTime = spice.spkezr("CHURYUMOV-GERASIMENKO",
                                    et, "J2000", "NONE", "SUN")

    x = state[0]
    y = state[1]
    z = state[2]
    vx = state[3]
    vy = state[4]
    vz = state[5]

    r_length = np.sqrt(x**2 + y**2 + z**2)

    r_hat = np.zeros(3)
    v = np.zeros(3)

    r_hat[0] = x / r_length
    r_hat[1] = y / r_length
    r_hat[2] = z / r_length

    v[0] = r_hat[0] * vx
    v[1] = r_hat[1] * vy
    v[2] = r_hat[2] * vz

    v_sun = v[0] + v[1] + v[2]

    return v_sun

def get_gfactor_from_db(v_sun=8.8, gasTemp=100, species='CO'):

    db = sqlite3.connect('test.sqlite')
    cur = db.cursor()

    DBqueryLowerV = 'SELECT v_sun, gFactor from gFactors WHERE name\
                    = "%s" AND gasTemp = %i AND v_sun < %f \
                    ORDER BY v_sun DESC' % (species, gasTemp, v_sun)
    DBqueryHigherV = 'SELECT v_sun, gFactor from gFactors WHERE name\
                     = "%s" AND gasTemp = %i AND v_sun > %f \
                     ORDER BY v_sun ASC' % (species, gasTemp, v_sun)

    cur.execute(DBqueryHigherV)
    v_high, gFactorHigh = cur.fetchall()[0]

    cur.execute(DBqueryLowerV)
    v_low, gFactorLow = cur.fetchall()[0]
    db.close()

    gFactorInterpolated = np.interp([v_sun], [v_low, v_high],
                                    [gFactorLow, gFactorHigh*2])
    gFactor = gFactorInterpolated
    print 'high :', v_high, gFactorHigh*2
    print 'low  :', v_low, gFactorLow
    print 'sun:', v_sun, gFactorInterpolated

    return gFactor