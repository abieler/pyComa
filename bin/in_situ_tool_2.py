#!/opt/local/anaconda/bin/python
from __future__ import division
import os
import sys
import datetime
import argparse
import time
import numpy as np
import sqlite3
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
from utils.rosettaDefs import *
from utils.batl import *

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

iDim = get_iDim(args)

###############################################
# get coordinates of rosetta in correct frame
###############################################
if args.iPointingCase == 0:
    # DSMC case
    if iDim < 3:
        frame = '67P/C-G_CSO'
    elif iDim == 3:
        frame = '67P/C-G_CK'
        #print '3D in-situ not implemented yet!'
        #sys.exit()
    print 'Rosetta coordinates in reference frame:', frame
    x_SC, y_SC, z_SC, r_SC, dates_SC = spice_functions.get_coordinates(args.StringUtcStartTime, args.StringKernelMetaFile,
                                                             'ROSETTA', frame, "None", "CHURYUMOV-GERASIMENKO",
                                                             args.StringUtcStopTime, args.nDeltaT)
    x_sun, y_sun, z_sun, r_sun, dates_sun = spice_functions.get_coordinates(args.StringUtcStartTime, args.StringKernelMetaFile,
                                                             'SUN', frame, "None", "CHURYUMOV-GERASIMENKO",
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

if args.iModelCase == dsmc_:
    # DSMC case
    path = os.path.split(args.StringDataFileDSMC)[0]
    dsmc_case = os.path.dirname(args.StringDataFileDSMC).split('/')[-1]
    
    if args.IsDust:
        filenames = [path + '/' + filename for filename in os.listdir(path) if (filename.split('.')[-1].lower() == 'dat') and 'Dust' in filename]
    else:
        filenames = [path + '/' + filename for filename in os.listdir(path) if (filename.split('.')[-1].lower() == 'dat') and 'Dust' not in filename]

elif args.iModelCase == haser_:
    filenames = ['Haser']
elif args.iModelCase == userModel_:
    filenames = [args.StringUserDataFile]
##############################################################
# load data (for 3D cases data is loaded from julia program)
##############################################################
if iDim < 3:
    if args.iModelCase == dsmc_:
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
                
    elif args.iModelCase == haser_:
        print 'haser case'
        x, n = haserModel(args.QHaser, args.vHaser, args.tpHaser, args.tdHaser, args.tgHaser)
        y = None
        numberDensities = [n]
    elif args.iModelCase == userModel_:
        print 'user case'
        x, y, n = load_user_data(args.StringUserDataFile, iDim, args.DelimiterData, args.nHeaderRowsData)
        numberDensities = [n]

    if args.SubCase == 'preferred':
        print 'preferred case'
        AU = np.float(args.StringDataFileDSMC.split('/')[-1].split('_')[1])
        print "AU:", AU
        if AU >= 3.0:
            numberDensities = [nn * 0.24 for nn in numberDensities]
        else:
            numberDensities = [nn * 1.25 for nn in numberDensities]
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
        pass 
        
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
if iDim == 3:
    if args.iModelCase == dsmc_:
        #####################################################
        # select calculate lon/lat of Sun for each instance
        # in time and check database for the best suited DSCMC case
        #####################################################
        db = sqlite3.connect(args.StringOutputDir+'/../../ICES.sqlite')
        cur = db.cursor()
        tup = (dsmc_case,)

        # runs is a dictionary with run names as keys,
        # holding at least 3 lists. these lists have
        # dates, coordinates and number densities that belong
        # to this run. if more than 1 species, more number
        # density lists are appended
        runs = {} 
        runNames = []

        # retrieve all possible run names for the selected case
        # that means all dsmc files for all AUs
        cur.execute("SELECT data_prefix FROM select3D WHERE dsmc_case=?", tup)
        queryData = cur.fetchall()
        for qd in queryData:
            runs[qd[0]] =  [[],[]]
            runNames.append(qd[0])

        # get all longitudes and latitudes of the sun for all
        # dates where data has to be interpolated 
        km2AU = 1.0 / 149597871
        lat_sun = []
        lon_sun = []
        r_AU = []
        shapeModel = 'rmoc'
        for xx,yy,zz in zip(x_sun, y_sun, z_sun):
            r, llon, llat = spice.reclat([xx,yy,zz])
            r = r*km2AU
            llon = llon / np.pi * 180
            llat = llat / np.pi * 180
            lat_sun.append(llat)
            lon_sun.append(llon)
            r_AU.append(r)

        print "selecting cases from lat lon of sun"
        kk = 0
        #for llat, llon in zip(lat_sun, lon_sun):
        for tt, llon in zip(dates_SC, lon_sun):
            # new version where case is selected by being closest in time
            sql = "SELECT data_prefix, abs(julianday(date(datetime))- julianday('%s')) AS 'delta_t' FROM select3D WHERE dsmc_case='%s' ORDER BY delta_t, abs( (((longitude-(%.1f)) + 180) %% 360) - 180) LIMIT 1;" %(tt, dsmc_case, llon);

            # old version where selected by latitude and longitude difference.
            #sql = "SELECT data_prefix FROM select3D WHERE dsmc_case='%s' ORDER BY abs( (((latitude-(%.1f)) + 180) %% 360) - 180), abs( (((longitude-(%.1f)) + 180) %% 360) - 180) LIMIT 1;" % (dsmc_case,llat, llon)

            cur.execute(sql)
            runName = cur.fetchone()[0]
            runs[runName][0].append(dates_SC[kk])
            runs[runName][1].append([x_SC[kk], y_SC[kk], z_SC[kk]])
            kk += 1

        species3D = []
        cur.execute("SELECT species from dsmc_species WHERE dcase='%s'"% dsmc_case)
        qData = cur.fetchall()
        for qd in qData:
            species3D.append(qd[0])

        db.close()

        # build path where all fitting DSMC cases are located below
        pathToData = os.path.dirname(args.StringDataFileDSMC)


        # write coordinates of each case to file and then start julia
        # to perform the 3D interpolation
        for key in runNames:
            if len(runs[key][0]) > 0:
                with open("rosettaCoords.txt", 'w') as oFile:
                    for rrr in runs[key][1]:
                        oFile.write("%.5e,%.5e,%.5e\n" %(rrr[0], rrr[1], rrr[2]))
                for spec in species3D:
                    runName = key+"." + spec + ".dat"
                    dsmcFileName = os.path.join(pathToData, runName) 
                    #os.system("su _www -c '/Applications/Julia-0.3.0.app/Contents/Resources/julia/bin/julia /Users/abieler/newLOS/in-situ.jl %s %s'" %(pathToData, args.StringOutputDir))
                    os.system("export JULIA_PKGDIR=/opt/local/share/julia/site ; /opt/local/bin/julia ../../../Models/LoS/pyComa/bin/in-situ.jl %s %s" %(dsmcFileName, args.StringOutputDir))
                    n_SC = np.genfromtxt('interpolation.out', dtype=float)

                    # genfromtxt returns float instead of one element array in case
                    # there is only one entry --> make array out of that
                    if type(n_SC) == float:
                        n_SC = np.array([n_SC])
                    runs[key].append(n_SC)
            else:
                pass

        # combine all cases into one array and sort them according to date
        numberDensities_SC = []
        for i in range(len(species3D)):
            n_SC = []
            dates_SC = []
            r_SC = []
            for key in runNames:
                if len(runs[key][0]) >= 1:
                    n_SC.extend(runs[key][2+i])
                    dates_SC.extend(runs[key][0])
                    r = [np.sqrt(p[0]**2+p[1]**2+p[2]**2) for p in runs[key][1]]
                    r_SC.extend(r)
            n_SC = np.array(n_SC)
            dates_SC = np.array(dates_SC)
            r_SC = np.array(r_SC)
            sort_index = np.argsort(dates_SC)

            dates_SC = dates_SC[sort_index]
            n_SC = n_SC[sort_index]
            r_SC = r_SC[sort_index]
            numberDensities_SC.append(n_SC)

    elif args.iModelCase == batsrus_:
        print "------------------------------------------------"
        print " 3D BATSRUS cases not officially available yet"
        print " please contact kenhan@umich.edu or abieler@umich.edu"
        print " to get early access to these results"
        print "------------------------------------------------"
        sys.exit()
        #batsDataFile  = "/www/ices/Data/Coma/BATSRUS/CG_1.3_au_00/CG_1.3_au_00.idl"
        nPoints = len(x_SC)
        myCoords = np.zeros((nPoints, iDim))
        i = 0
        for xx,yy,zz in zip(x_SC, y_SC, z_SC):
            myCoords[i,0] = xx 
            myCoords[i,1] = yy
            myCoords[i,2] = zz
            i += 1

        S, varNames = interpolate(args.StringDataFileDSMC, myCoords)
        print S.shape
        print len(varNames)
        numberDensities_SC = []
        for iVar in range(S.shape[1]):
            numberDensities_SC.append(S[:,iVar])

############################################
# write results to file
############################################
if args.iModelCase == dsmc_:
    if args.IsDust:
        species = [size for size in allSizeIntervals if args.DustSizeMin <= size <= args.DustSizeMax]
    else:
        species = []
        for filename in filenames:
            spec = filename.split('.')[-2]
            if spec not in species:
                species.append(spec)
        #species = [filename.split('.')[-2] for filename in filenames]

    x_SC *= -1          # transform back from tenishev to cso frame of reference
    y_SC *= -1          # transform back from thenisev to cso frame of reference
elif args.iModelCase == batsrus_:
    species = [varName for varName in varNames]
elif args.iModelCase == haser_:
    species = ['Haser']
elif args.iModelCase == userModel_:
    species = ['User']

#file = open(args.StringOutputDir + '/' + 'in_situ' + '.out', 'w')
with open(args.StringOutputDir + '/' + 'in_situ' + '.out', 'w') as file:
    file.write('Local number densities for the rosetta spacecraft at selected dates. Comet is at (0,0,0) with the sun on the positive x axis.(inf,0,0)\n')
    if args.iModelCase == dsmc_:
        file.write('DSMC case: %s\n' % (os.path.split(args.StringDataFileDSMC)[0].split('/')[-1]))
    elif args.iModelCase == batsrus_:
        file.write('BATSRUS case: %s\n' % (os.path.split(args.StringDataFileDSMC)[0].split('/')[-1]))
        
    elif args.iModelCase == haser_:
        file.write("HASER case: Q = %.3e [#/s], v = %f [m/s], tp = %.2e" % (args.QHaser, args.vHaser, args.tpHaser))
        if args.tdHaser == 0:
            file.write('\n')
        else:
            file.write(', %.2e\n' %(args.tdHaser))

    if args.iPointingCase == spice_:
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
