#!/usr/bin/env python
from __future__ import division
import numpy as np

def haserModel(Q, v, tp, td=False):
    '''
    Q   : gas production rate of parent
    v   : radial outflow velocity parent
    tp  : lifetime parent
    td  : lifetime daughter (if any)
    '''
    
    rmin = np.log10(2000)
    rmax = np.log10(1e8)
    dr = (rmax - rmin)/ 2000
    r = 10**(np.arange(rmin, rmax, dr))
    
    gamma_p = v * tp
    
    if td == None:
        n = Q / (4 * np.pi * r**2 * v) * np.exp(-r/gamma_p)
    else:
        gamma_d = v * td
        n = Q / (4 * np.pi * r**2 * v) * (gamma_d / (gamma_p - gamma_d)) * (np.exp(-r/gamma_p) - np.exp(-r/gamma_d))
    
    return r, n


