from __future__ import division
import os
import sys
import datetime
import argparse
import numpy as np
import matplotlib.pyplot as plt
try:
    import matplotlib.tri as mtri
    triImport = True
except:
    print "Error importing tirangulation module matplotlib.tri"
    triImport = False
import matplotlib

import spice
from data_loaders import *
import spice_functions
from data_plotting import plot_in_situ
from haser import haserModel


parser = argparse.ArgumentParser()
parser.add_argument("--iModelCase", type=int, choices=[0, 1, 2], help='0: dsmc model, 1: haser model, 2: user model')
parser.add_argument("--iPointingCase", type=int, choices=[0, 1], help='0: spice pointing, 1: user pointing')
parser.add_argument("--iInstrumentSelector", type=int, choices=[1, 2, 3, 4, 5, 6])
parser.add_argument("--StringOutputDir", type=str)

parser.add_argument("--StringDSMCdir", type=str)
parser.add_argument("--StringDataFileDSMC", type=str)
parser.add_argument("--IsDust", type=int, choices=[0, 1], help='1 for dust case, 0 for gas case')
parser.add_argument("--DustSizeMin", type=float)
parser.add_argument('--DustSizeMax', type=float)


parser.add_argument("--QHaser", type=float)
parser.add_argument("--vHaser", type=float)
parser.add_argument("--tpHaser", type=float)
parser.add_argument("--tdHaser", type=float)

parser.add_argument("--StringUserDataFile", type=str)                          # file to upload from user which contains user coma model
parser.add_argument("--UserDelimiter", type=str)                         # delimiter used in datafile
parser.add_argument("--iUserNrOfHeaderRows", type=int)                   # number of header lines in datafile
parser.add_argument("--iUserDim", type=int)                               # number of dimenisons of user coma model

parser.add_argument("--StringKernelMetaFile", type=str)
parser.add_argument("--StringUtcStartTime", type=str)
parser.add_argument("--StringUtcStopTime", type=str)
parser.add_argument("--nDeltaT", type=int)

parser.add_argument("--UserR", type=float)                               # Distance in km from nucleus center
parser.add_argument("--UserPhaseAngle", type=float)
parser.add_argument("--UserLatitude", type=float)
parser.add_argument("--UserAlpha", type=float)
parser.add_argument("--UserBeta", type=float)
parser.add_argument("--UserGamma", type=float)

args = parser.parse_args()

print '--' * 20
print "mkFile:", args.StringKernelMetaFile
print "start:", args.StringUtcStartTime
print "stop:", args.StringUtcStopTime
print 'deltaT:', args.nDeltaT
print 'outDir:', args.StringOutputDir
print '--' * 20

x_SC, y_SC, z_SC, r_SC, dates_SC = spice_functions.get_coordinates(args.StringUtcStartTime, args.StringKernelMetaFile,
                                                         'ROSETTA', 'J2000', "None", "CHURYUMOV-GERASIMENKO",
                                                         args.StringUtcStopTime, args.nDeltaT)

os.system('rm ' + args.StringOutputDir + '/*.out' )

iDim = get_iDim(args)

if args.iModelCase == 0:
    path = os.path.split(args.StringDataFileDSMC)[0]
    filenames = [path + '/' + filename for filename in os.listdir(path) if (filename.split('.')[-1].lower() == 'dat') and 'Dust' not in filename]
    print filenames
    print '----------------------------------------------------'

elif args.iModelCase == 1:
    filenames = ['Haser']
elif args.iModelCase == 2:
    filenames = [args.StringUserDataFile]

for filename in filenames:

    ########################################################
    # load data
    ########################################################
    if args.iModelCase == 0:
        print 'dsmc case'
        if args.IsDust:
            print 'dust'
            NumberDensityIndices, allSizeIntervals = getAllDustIntervalIndices(filename, iDim)
            x, y, n = loadDustData(allSizeIntervals, NumberDensityIndices, iDim, filename)
        else:
            print 'gas'
            x, y, n = loadGasData(filename, iDim)
    elif args.iModelCase == 1:
        print 'haser case'
        x, n = haserModel(args.QHaser, args.vHaser, args.tpHaser, args.tdHaser)
        y = None
    elif args.iModelCase == 2:
        print 'user case'
        x, y, n = loadGasData(args.StringUserDataFile, iDim, True, args.UserDelimiter, args.iUserNrOfHeaderRows)

    ##############################################################
    # triangulation and interpolation for 2d case
    if iDim == 1:
        pass
    elif iDim == 2:
        Triangles = mtri.Triangulation(x, y)
        Interpolator = mtri.LinearTriInterpolator(Triangles, n)

    print 'interpolation done'
    #############################################################

    if iDim == 1:
        n_SC = np.interp(np.sqrt(x_SC**2 + y_SC**2 + z_SC**2), x, n)
    elif iDim == 2:
        n_SC = Interpolator.__call__(x_SC, np.sqrt(y_SC**2 + z_SC**2))
    elif iDim == 3:
        n_SC = 0

    ############################################
    # write results to file
    ############################################
    if args.iModelCase == 0:
        species = filename.split('.')[-2]
    elif args.iModelCase == 1:
        species = 'Haser'
    elif args.iModelCase == 2:
        species = 'User'

    file = open(args.StringOutputDir + '/' + species + '.out', 'w')
    file.write('Local number densities for the rosetta spacecraft at selected dates. Comet is at (0,0,0) with the sun on the positive x axis.(inf,0,0)\n')
    file.write('DSMC case: %s\n' % (os.path.split(args.StringDataFileDSMC)[0].split('/')[-1]))
    file.write('spice kernel: %s\n' % (args.StringKernelMetaFile.split('/')[-1]))
    file.write('date,x[m],y[m],z[m],distance_from_center[m],numberDensity [1/m3]\n')
    for dd, xx, yy, zz, rr, nn in zip(dates_SC, x_SC, y_SC, z_SC, r_SC, n_SC):
        file.write("%s,%e,%e,%e,%e,%e\n" % (dd, xx, yy, zz, rr, nn))
    file.close()
    print 'done'

#######################################################
# plot results
#######################################################
plot_in_situ(args)