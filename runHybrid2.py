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
import time

import spice
from utils.cmdline_args import cmdline_args
from utils.createRay import createRay
from utils.data_loaders import load_hybrid2_data
import utils.spice_functions as spice_functions
from utils.data_plotting import plot_result_insitu

parser = argparse.ArgumentParser()
args = utils.cmdline_args(parser)

pathToExecutable = '/Users/ices/www-dev/htdocs/ICES/Models/LoS/pyComa'

if args.StringMeasurement == 'LOS':
    #####################################################################
    # get position of s/c and pointing towards Earth from SPICE
    #####################################################################
    spice.furnsh(args.StringKernelMetaFile)
    Et = spice.str2et(args.StringUtcStartTime)
    rEarth, lightTime = spice.spkpos("EARTH", Et, "67P/C-G_CSO", "NONE", "ROSETTA")
    rRosetta, lightTime = spice.spkpos("ROSETTA", Et, "67P/C-G_CSO", "NONE", "CHURYUMOV-GERASIMENKO")
    rEarth = np.array(rEarth)
    rRosetta = np.array(rRosetta) * 1000  # transform km to m
    print 'Distance from comet: %.2e' % (np.sqrt(np.sum(rRosetta ** 2)))

    # p = normalized vector pointing from s/c to Earth in cso coordinates
    p = rEarth / np.sqrt(np.sum(rEarth**2))

    # rRay = position of rosetta s/c in cso coordinates
    rRay = np.array([value for value in rRosetta])

    ######################################################################
    # build line of sight ray xTravel, and extract x, y, z coordinates
    #######################################################################
    xTravel = np.array(createRay(rRay, p))
    x = xTravel[:, 0]
    y = xTravel[:, 1]
    z = xTravel[:, 2]

elif args.StringMeasurement == 'insitu':

    x, y, z, r, dates = spice_functions.get_coordinates(args.StringUtcStartTime, args.StringKernelMetaFile,
                                                             'ROSETTA', '67P/C-G_CSO', "None", "CHURYUMOV-GERASIMENKO",
                                                             args.StringUtcStopTime, args.nDeltaT)

#############################################################################
# write coordinates to traj.dat file and execute aikef.py
##############################################################################
with open('traj.dat', 'w') as file:
    file.write('#START\n')
    for xx, yy, zz in zip(x, y, z):
        file.write('%e %e %e\n' % (xx, yy, zz))
print 'running aikef script now'
os.system(pathToExecutable + '/aikef_andre.py %s run . traj.dat' % (args.StringHybridCase))
print 'done running script'

##################################################################
# load output from aikef.py (hybrid2) and perform LOS calculation
##################################################################
time.sleep(3)
xTravelRay, density = load_hybrid2_data('orbit-output.txt')

if args.StringMeasurement == 'LOS':

    xTravelRay = np.array([np.array([xx, yy, zz]) for xx, yy, zz in zip(xTravelRay[0, :],
                                                                        xTravelRay[1, :],
                                                                        xTravelRay[2, :])])
    dTravelRay = np.sqrt(np.sum((xTravelRay[0] - xTravelRay)**2, axis=1))
    columnDensity = np.trapz(density, dTravelRay)

    with open(args.StringOutputDir + '/result.txt', 'w') as f:
        f.write('Electron column density between s/c and Earth in [#/m3] as calculated from the Hybrid2 model.\n')
        f.write('Case : %s\n' % (args.StringHybridCase))
        f.write('Spice: %s\n' % (args.StringKernelMetaFile))
        f.write('Date : %s\n' % (args.StringUtcStartTime))
        f.write('ColumnDensity: %.3e\n' % (columnDensity))
    print "electron column density: %.3e" % (columnDensity)

elif args.StringMeasurement == 'insitu':
    with open(args.StringOutputDir + '/' + 'electrons' + '.out', 'w') as f:
        f.write('Local number densities for the rosetta spacecraft at selected dates. Comet is at (0,0,0) with the sun on the positive x axis.(inf,0,0)\n')
        f.write('DSMC case: %s\n' % (os.path.split(args.StringDataFileDSMC)[0].split('/')[-1]))
        f.write('spice kernel: %s\n' % (args.StringKernelMetaFile.split('/')[-1]))
        f.write('date,x[m],y[m],z[m],distance_from_center[m],numberDensity [1/m3]\n')
        for dd, xx, yy, zz, rr, nn in zip(dates, x, y, z, r, density):
            f.write("%s,%e,%e,%e,%e,%e\n" % (dd, xx, yy, zz, rr, nn))

    plot_result_insitu(args)
