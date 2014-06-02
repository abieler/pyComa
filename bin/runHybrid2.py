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
args = cmdline_args(parser)

pathToExecutable = '/Users/ices/www-dev/htdocs/ICES/Models/LoS/pyComa/bin'
#pathToExecutable = '/Users/abieler/pyComa/bin'

if args.StringMeasurement == 'LOS':
    print 'LOS hybrid case...'
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
    print 'in-situ case hybrid...'
    if args.iPointingCase == 0:
        x, y, z, r, dates = spice_functions.get_coordinates(args.StringUtcStartTime,
                                                            args.StringKernelMetaFile,
                                                            'ROSETTA',
                                                            '67P/C-G_CSO',
                                                            "None",
                                                            "CHURYUMOV-GERASIMENKO",
                                                            args.StringUtcStopTime,
                                                            args.nDeltaT)
    elif args.iPointingCase == 2:
        x, y, z = load_user_trajectory(args)
        r = np.array([np.sqrt(xx**2 + yy**2 + zz**2) for xx, yy, zz in zip(x,y,z)])

#############################################################################
# write coordinates to traj.dat file and execute aikef.py
##############################################################################
with open('traj.dat', 'w') as file:
    file.write('#START\n')
    for xx, yy, zz in zip(x, y, z):
        file.write('%e %e %e\n' % (xx, yy, zz))
print 'running aikef script now'
os.system(pathToExecutable + '/aikef_for_pyComa.py %s run . traj.dat' % (args.StringHybridCase))
print 'done running script'

##################################################################
# load output from aikef.py (hybrid2) and perform LOS calculation
##################################################################
xTravelRay, B_total, B, U_e, nElectrons, nIons, U_ions = load_hybrid2_data('orbit-output.txt')

B_hat = [B/np.sqrt((B**2).sum())]  # B normalized
v_perp_B = np.array([np.cross(b_hat, np.cross(b_hat, v)) for b_hat, v in zip(B_hat, U_e)])  # v component perpendicular to B
v_perp_B_ion = np.array([np.cross(b_hat, np.cross(b_hat, v)) for b_hat, v in zip(B_hat, U_ions)])  # v component perpendicular to B
v_perp_B_abs = np.sqrt((v_perp_B**2).sum())  # length of v component perpendicular to B
v_perp_B_ion_abs = np.sqrt((v_perp_B_ion**2).sum()) 

Q_ELECTRON = 1.602176 * 10**-19
M_ELECTRON = 9.109383 * 10**-31
AMU = 1.6605 * 10**-27            # atomic mass unit [kg]
Na = 6.02214129 * 10**23          # avogadro's number 1/mol
mH2O = 18.01528                   # water molecule mass [g/mol]
mIon_kg = mH2O / Na / 1000        # mass of ion molecule in kg
M_ION = mIon_kg
print 'M_ION: ', M_ION
EPSILON_0 = 8.8541878 * 10**-12

rg_electron = M_ELECTRON * v_perp_B_abs / Q_ELECTRON / B_total
rg_ion = M_ION * v_perp_B_ion_abs / Q_ELECTRON / B_total

wg_electron = Q_ELECTRON * B_total / M_ELECTRON
wg_ion = Q_ELECTRON * B_total / M_ION

wp_electron = np.sqrt(nElectrons * Q_ELECTRON**2 / (M_ELECTRON * EPSILON_0))
wp_ion = np.sqrt(nIons * Q_ELECTRON**2 / (M_ION * EPSILON_0))



if args.StringMeasurement == 'LOS':

    xTravelRay = np.array([np.array([xx, yy, zz]) for xx, yy, zz in zip(xTravelRay[0, :],
                                                                        xTravelRay[1, :],
                                                                        xTravelRay[2, :])])
    dTravelRay = np.sqrt(np.sum((xTravelRay[0] - xTravelRay)**2, axis=1))
    columnDensity = np.trapz(density, dTravelRay)

    with open(args.StringOutputDir + '/result.txt', 'w') as f:
        f.write(('Electron column density between s/c and Earth in [#/m3] as' +
                 ' calculated from the Hybrid2 model.\n'))
        f.write('Case : %s\n' % (args.StringHybridCase))
        f.write('Spice: %s\n' % (args.StringKernelMetaFile))
        f.write('Date : %s\n' % (args.StringUtcStartTime))
        f.write('ColumnDensity: %.3e\n' % (columnDensity))
    print "electron column density: %.3e" % (columnDensity)

elif args.StringMeasurement == 'insitu':
    with open(args.StringOutputDir + '/' + 'electrons' + '.out', 'w') as f:
        f.write(('Local plasma parameters at rosetta spacecraft for selected dates.'
                ' Comet center is at (0,0,0) with the sun on the positive x axis.\n'))
        f.write('AIKEF Hybrid case: %s\n' % args.StringHybridCase)
        if args.iPointingCase == 0:
            f.write('spice kernel: %s\n' % (args.StringKernelMetaFile.split('/')[-1]))
        elif args.iPointingCase == 2:
            f.write('User defined trajectory: %s' % args.StringUserTrajectoryFile)
        f.write('date x[m] y[m] z[m] distance_from_center[m] numberDensity_e[1/m3] numberDensity_ions[1/m3] B_total[nT] gyroRadius_e[m] gyroFreq_e[rad/s] plasmaFreq_e[rad/s] gyroRadius_ion gyroFreq_ion plasmaFreq_ion\n')
        for dd, xx, yy, zz, rr, ne, bb, rge, wge, wpe, rgi, wgi, wpi  in zip(dates, x, y, z, r, nElectrons, nIons, B_total,
                                                                             rg_electron, wg_electron, wp_electron,
                                                                             rg_ion, wg_ion, wp_ion ):
            f.write("%s %e %e %e %e %e %e %e %e %e %e %e %e %e\n" % (dd, xx, yy, zz, rr, ne, bb, rge, wge, wpe, rgi, wgi, wpi))

    plot_result_insitu(args)
