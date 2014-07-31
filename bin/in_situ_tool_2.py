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
from utils.data_plotting import plot_result_insitu_2
from utils.haser import haserModel
from utils.cmdline_args import cmdline_args

t0 = time.time()
parser = argparse.ArgumentParser()
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
##############################################################
# load data
##############################################################
if args.iModelCase == 0:
    x_SC *= -1      # cso reference frame to tenishev reference frame
    y_SC *= -1      # cso reference frame to tenishev reference frame
    if args.IsDust:
        print 'dsmc case, dust'
        NumberDensityIndices, allSizeIntervals = getAllDustIntervalIndices(args.StringDataFileDSMC, iDim)
        x, y, numberDensities, massDensities = load_dust_data_full(allSizeIntervals, NumberDensityIndices,
                                                                   iDim, args.StringDataFileDSMC, args)
    else:
        print 'dsmc case, gas'
        numberDensities = []
        for filename in filenames:
            x, y, n = loadGasData(filename, iDim)
            numberDensities.append(n)
            
elif args.iModelCase == 1:
    print 'haser case'
    x, n = haserModel(args.QHaser, args.vHaser, args.tpHaser, args.tdHaser)
    y = None
    numberDensities = [n]
elif args.iModelCase == 2:
    print 'user case'
    x, y, n = load_user_data(args.StringUserDataFile, iDim, args.DelimiterData, args.nHeaderRowsData)
    numberDensities = [n]

##############################################################
# triangulation and interpolation for 2d case
if iDim == 1:
    pass
elif iDim == 2:
    Triangles = mtri.Triangulation(x, y)
    if args.IsDust:
        nSpecies = numberDensities.shape[1]
        numberDensities = [numberDensities[:,i] for i in range(nSpecies)]
    Interpolator = [mtri.LinearTriInterpolator(Triangles, n) for n in numberDensities]
elif iDim == 3:
    print '3d not implemented yet'
    
print 'interpolation done'
numberDensities_SC = []
#############################################################
i=0
for n in numberDensities:
    if iDim == 1:
        n_SC = np.interp(np.sqrt(x_SC**2 + y_SC**2 + z_SC**2), x, n)
    elif iDim == 2:
        n_SC = Interpolator[i].__call__(x_SC, np.sqrt(y_SC**2 + z_SC**2))
    elif iDim == 3:
        n_SC = 0
    numberDensities_SC.append(n_SC)
    i+=1

############################################
# write results to file
############################################
if args.iModelCase == 0:
    if args.IsDust:
        species = [size for size in allSizeIntervals if args.DustSizeMin <= size <= args.DustSizeMax]
    else:
        species = [filename.split('.')[-2] for filename in filenames]

    x_SC *= -1          # transform back from tenishev to cso frame of reference
    y_SC *= -1          # transform back from thenisev to cso frame of reference
elif args.iModelCase == 1:
    species = ['Haser']
elif args.iModelCase == 2:
    species = ['User']

#file = open(args.StringOutputDir + '/' + 'in_situ' + '.out', 'w')
with open(args.StringOutputDir + '/' + 'in_situ' + '.out', 'w') as file:
    file.write('Local number densities for the rosetta spacecraft at selected dates. Comet is at (0,0,0) with the sun on the positive x axis.(inf,0,0)\n')
    if args.iModelCase == 0:
        file.write('DSMC case: %s\n' % (os.path.split(args.StringDataFileDSMC)[0].split('/')[-1]))
    elif args.iModelCase == 1:
        file.write("HASER case: Q = %.3e [#/s], v = %f [m/s], tp = %.2e" % (args.QHaser, args.vHaser, args.tpHaser))
        if args.tdHaser == 0:
            file.write('\n')
        else:
            file.write(', %.2e\n' %(args.tdHaser))

    if args.iPointingCase == 0:
        file.write('spice kernel: %s\n' % (args.StringKernelMetaFile.split('/')[-1]))

    file.write('date,x[m],y[m],z[m],distance_from_center[m],')
    for s in species:
        if s == species[-1]:
            file.write('%s [#/m3]' %s)
        else:
            file.write('%s [#/m3],' %s)
    file.write('\n')

    i = 0
    for dd, xx, yy, zz, rr, in zip(dates_SC, x_SC, y_SC, z_SC, r_SC):
        file.write("%s,%e,%e,%e,%e," % (dd, xx, yy, zz, rr))
        for n_SC in numberDensities_SC:
            if (np.sum(n_SC == numberDensities_SC[-1]) == len(n_SC)):
                file.write("%e" % n_SC[i])
            else:
                file.write("%e," % n_SC[i])
        file.write("\n")
        i += 1
    file.close()
print 'done'

#######################################################
# plot results
#######################################################
plot_result_insitu_2(args, numberDensities_SC, species, dates_SC, r_SC)
print 'Wall time:', time.time() - t0
