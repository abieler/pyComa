from __future__ import division
import os
import sys
import datetime
import argparse
import time
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
from utils.data_loaders import *
import utils.spice_functions as spice_functions
from utils.data_plotting import plot_result_insitu
from utils.haser import haserModel
from utils.cmdline_args import cmdline_args

t0 = time.time()
parser = utils.argparse.ArgumentParser()
args = cmdline_args(parser)

print '--' * 20
print "mkFile:", args.StringKernelMetaFile
print "start:", args.StringUtcStartTime
print "stop:", args.StringUtcStopTime
print 'deltaT:', args.nDeltaT
print 'outDir:', args.StringOutputDir
print '--' * 20

if args.iPointingCase == 0:
    x_SC, y_SC, z_SC, r_SC, dates_SC = spice_functions.get_coordinates(args.StringUtcStartTime, args.StringKernelMetaFile,
                                                             'ROSETTA', '67P/C-G_CSO', "None", "CHURYUMOV-GERASIMENKO",
                                                             args.StringUtcStopTime, args.nDeltaT)
elif args.iPointingCase == 1:
    print 'user pointing is not a valid choice for in situ tool!'
    print 'for user defined trajectory select iPointingCase = 2'
    print 'exiting now'
    sys.exit()
elif args.iPointingCase == 2:
    print 'upload user data file...'
    x_SC, y_SC, z_SC = load_user_trajectory(args)
    dates_SC = np.arange(len(x_SC))
    r_SC = np.sqrt(x_SC**2 + y_SC**2 + z_SC**2)

os.system('rm ' + args.StringOutputDir + '/*.out')

iDim = get_iDim(args)

if args.iModelCase == 0:
    path = os.path.split(args.StringDataFileDSMC)[0]
    if args.IsDust:
        filenames = [path + '/' + filename for filename in os.listdir(path) if (filename.split('.')[-1].lower() == 'dat') and 'Dust' in filename]
    else:
        filenames = [path + '/' + filename for filename in os.listdir(path) if (filename.split('.')[-1].lower() == 'dat') and 'Dust' not in filename]
    print filenames
    print '----------------------------------------------------'

elif args.iModelCase == 1:
    filenames = ['Haser']
elif args.iModelCase == 2:
    filenames = [args.StringUserDataFile]

Triangles = None
for filename in filenames:

    if args.iModelCase == 0:
        x_SC *= -1      # cso reference frame to tenishev reference frame
        y_SC *= -1      # cso reference frame to tenishev reference frame
    ########################################################
    # load data
    ########################################################
    if args.iModelCase == 0:
        print 'dsmc case'
        if args.IsDust:
            print 'dust'
            NumberDensityIndices, allSizeIntervals = getAllDustIntervalIndices(filename, iDim)
            x, y, n = load_dust_data(allSizeIntervals, NumberDensityIndices, iDim, filename, args)
        else:
            print 'gas'
            x, y, n = loadGasData(filename, iDim)
            #x, y, n = load_gas_data(filename, iDim)
    elif args.iModelCase == 1:
        print 'haser case'
        x, n = haserModel(args.QHaser, args.vHaser, args.tpHaser, args.tdHaser)
        y = None
    elif args.iModelCase == 2:
        print 'user case'
        #x, y, n = loadGasData(args.StringUserDataFile, iDim, True, args.DelimiterData, args.nHeaderRowsData)
        x, y, n = load_user_data(args.StringUserDataFile, iDim, args.DelimiterData, args.nHeaderRowsData)

    ##############################################################
    # triangulation and interpolation for 2d case
    if iDim == 1:
        pass
    elif iDim == 2:
        if not Triangles:
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
        x_SC *= -1          # transform back from tenishev to cso frame of reference
        y_SC *= -1          # transform back from thenisev to cso frame of reference
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
plot_result_insitu(args)
print 'Time elapsed:', time.time() - t0
