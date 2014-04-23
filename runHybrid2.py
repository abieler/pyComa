'''
python wrapper for the hybrid2 ICES tool. combines
the line of sight capability with the hybrid2 model.
calculates column density of electrons between rosetta
space craft and the earth.
'''
import numpy as np
import subprocess
import argparse
import os
import matplotlib.pyplot as plt

import spice
from cmdline_args import cmdline_args
#from create_ray_noCython import createRay
from createRay import createRay
from data_loaders import load_hybrid2_data

parser = argparse.ArgumentParser()
args = cmdline_args(parser)

pathToExecutable = '/Users/ices/www-v4.1/htdocs/ICES/Models/Hybrid2'

spice.furnsh(args.StringKernelMetaFile)
Et = spice.str2et(args.StringUtcStartTime)
rEarth, lightTime = spice.spkpos("EARTH", Et, "67P/C-G_CSO", "NONE", "ROSETTA")        # s/c coordinates in CSO frame of reference
rRosetta, lightTime = spice.spkpos("ROSETTA", Et, "67P/C-G_CSO", "NONE", "CHURYUMOV-GERASIMENKO")
rEarth = np.array(rEarth)                           # transform km to m
rRosetta = np.array(rRosetta) * 1000
print 'Distance from comet: %.2e' % (np.sqrt(np.sum(rRosetta ** 2)))

p = rEarth / np.sqrt(np.sum(rEarth**2))             # normalized vector pointing from s/c to Earth in cso coordinates
rRay = np.array([value for value in rRosetta])      # position of rosetta s/c in cso coordinates

xTravel = np.array(createRay(rRay, p))
print xTravel
x = xTravel[:, 0]
y = xTravel[:, 1]
z = xTravel[:, 2]
with open(args.StringRuntimeDir + '/' + 'traj.dat', 'w') as file:
    file.write('#START\n')
    for xx, yy, zz in zip(x, y, z):
        file.write('%e %e %e\n' % (xx, yy, zz))
#os.chdir(args.StringRuntimeDir)
os.system(pathToExecutable + '/aikef.py %s run . traj.dat' % (args.StringHybridCase))

##############################################################
# load data and perform LOS
##############################################################

xTravelRay, density = load_hybrid2_data('orbit-output.txt')
xTravelRay = np.array([np.array([xx,yy,zz]) for xx,yy,zz in zip(xTravelRay[0,:], xTravelRay[1,:], xTravelRay[2,:])])
dTravelRay = np.sqrt(np.sum((xTravelRay[0] - xTravelRay)**2,axis=1))

columnDensity = np.trapz(density, dTravelRay)
with open('result.txt','w') as f:
    f.write('Electron column density between s/c and Earth in [#/m3] as calculated from the Hybrid2 model.\n')
    f.write('Case : %s\n' % (args.StringHybridCase))
    f.write('Spice: %s\n' % (args.StringKernelMetaFile))
    f.write('Date : %s\n' % (args.StringUtcStartTime))
    f.write('ColumnDensity: %e\n' % (columnDensity))
print "columnDensity: %.2e" %(columnDensity)
