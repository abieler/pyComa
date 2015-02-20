#!/usr/bin/env python
from __future__ import division

def get_specs():

    PhiX = 2.5 
    PhiY = 2.5 
    iFOV = 0.1 

    instrumentSpecs = {'name' : 'sun',
                       'nPixelsX' : 400,
                       'nPixelsY' : 400,
                       'nOversampleX' : 1,
                       'nOversampleY' : 1,
                       'PhiX' : PhiX,
                       'PhiY' : PhiY,
                       'iFOV' : iFOV, 
                       'PixelSize' : 1,
                       'computedQuantity' : 'column density [#/m2]',
                       'frame':'ROS_ALICE'
                      }

    instrumentSpecs['PixelFOV'] = 0 

    return instrumentSpecs
