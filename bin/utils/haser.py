#!/usr/bin/env python
from __future__ import division
import numpy as np

# test comment
def haserModel(Q, v, tp, td=0, tg=0):
    '''
    Q   : gas production rate of parent
    v   : radial outflow velocity parent
    tp  : lifetime parent
    td  : lifetime daughter (if any)
    tg  : lifetime granddaughtter (if any)
    '''

    rmin = np.log10(2000)
    rmax = np.log10(1e9)
    dr = (rmax - rmin) / 2000
    r = 10**(np.arange(rmin, rmax, dr))

    gamma_p = v * tp

    if tp > 0 and td == 0 and tg == 0:
        n = Q / (4 * np.pi * r**2 * v) * np.exp(-r/gamma_p)
    elif tp > 0 and td > 0 and tg == 0:
        gamma_d = v * td
        if (gamma_d == gamma_p):
            gamma_d = gamma_d - (gamma_d / 1000)
        n = Q / (4 * np.pi * r**2 * v) * (gamma_d / (gamma_p - gamma_d)) * (np.exp(-r/gamma_p) - np.exp(-r/gamma_d))
    elif tp > 0 and td > 0 and tg > 0: 
        gamma_d = v * td
        gamma_g = v * tg
        if gamma_p == gamma_d:
            gamma_d = gamma_d - (gamma_d / 5000)
        if gamma_p == gamma_g:
            gamma_g = gamma_g - (gamma_g / 5000)
        if gamma_g == gamma_d:
            gamma_g = gamma_g - (gamma_g / 5000)
        A = (gamma_g/(gamma_p - gamma_g)) * (1 + (gamma_d/(gamma_p-gamma_d)))
        B = -A + (gamma_g * gamma_d / ((gamma_p - gamma_d) * (gamma_d - gamma_g)))
        C = -A - B

        n = Q / (4 * np.pi * r**2 * v) * (A * np.exp(-r/gamma_p) + B * np.exp(-r/gamma_g) + C * np.exp(-r/gamma_d))
    return r, n


