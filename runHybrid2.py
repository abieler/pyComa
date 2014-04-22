#!/opt/local/bin/python2.7
import numpy as np
import subprocess
import argparse
import os

import spice
from cmdline_args import cmdline_args
from create_ray_noCython import createRay
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
os.chdir(args.StringRuntimeDir)
os.system(pathToExecutable + '/aikef.py %s run . traj.dat' % (args.StringHybridCase))

#runTimePath = '/Users/ices/www-v4.1/htdocs/ICES/Runtime/hybrid2TestAndre'
#os.system('/Users/ices/www-v4.1/htdocs/ICES/Models/Hybrid2/aikef.py CG_2.5_au_02 run . traj.dat')

##############################################################
# load data and perform LOS
##############################################################

xTravelRay, density = load_hybrid2_data('output-orbit.txt')

xTravelRay = np.array([np.array([xx, yy, zz]) for xx, yy, zz in zip(xTravelRay[:, 0], xTravelRay[:, 1], xTravelRay[:, 2])])
dTravelRay = np.sqrt((xTravelRay[0] - xTravelRay)**2)

columnDensity = np.trapz(dTravel, density)
