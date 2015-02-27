#!/opt/local/anaconda/bin/python
'''
filename: pyLOS_vectorized.py

requires pySPICE installed

requires matplotlib version 1.3.1 or higher

***************************************************************
requires mpi4py installed (on osx: anaconda mpi4py 
package is broken for version 1.9.2, create a symlink:
sudo ln -s /yourPathTo/anaconda /opt/anaconda1anaconda2anaconda3
solves the problem
****************************************************************

***************************************************
requires cython installed.
--> execute:

python setupOSX.py build_ext --inplace

in the pyComa directory, this compiles the cython
module to build createRay.so which is then imported
by pyLOS_vectorized.py
****************************************************

'''
# must be first line of program, sets default division to floating point
from __future__ import division  

try:
    # import standard modules
    import sys
    import os
    import time
    import subprocess

    # import 3rd party modules
    import spice
    import numpy as np
    from numpy import dot
    import matplotlib.tri as mtri
    import matplotlib.pyplot as plt
    from mpi4py import MPI
    import argparse

    # import own modules
    from utils.data_loaders import *
    from utils.haser import haserModel
    from utils.dust import dustModel
    from utils.data_plotting import plot_result_LOS, build_plot_title, plot_miro
    import utils.rotations as rotations
    import utils.alice as alice
    import utils.miro as miro
    import utils.sun as sun
    import utils.createRay as createRay
    from utils.cmdline_args import cmdline_args
    from utils.instrument import Instrument

    # import defined constants
    from utils.rosettaDefs import *

except Exception, e:
    print '--' * 20
    print "Error with import of python module:"
    print e
    print '--' * 20
    sys.exit()

startTime = time.time()

#############################################
# setup argparser for cmd line arguments
##############################################
parser = argparse.ArgumentParser()
args = cmdline_args(parser)


iModelCase = args.iModelCase
iPointingCase = args.iPointingCase
iInstrumentSelector = args.iInstrumentSelector
StringOutputDir = args.StringOutputDir

StringDataFileDSMC = args.StringDataFileDSMC
IsDust = args.IsDust
DustSizeMin = args.DustSizeMin
DustSizeMax = args.DustSizeMax

QHaser = args.QHaser
vHaser = args.vHaser
tpHaser = args.tpHaser
tdHaser = args.tdHaser

QDust = args.QDust
vDust = args.vDust

StringKernelMetaFile = args.StringKernelMetaFile
StringUtcStartTime = args.StringUtcStartTime

UserR = args.UserR
UserPhaseAngle = args.UserPhaseAngle
UserLatitude = args.UserLatitude
UserAlpha = args.UserAlpha
UserBeta = args.UserBeta
UserGamma = args.UserGamma

StringUserDataFile = args.StringUserDataFile
DelimiterData = args.DelimiterData
iDimUser = args.iDimUser
nHeaderRowsData = args.nHeaderRowsData

comm = MPI.COMM_WORLD
nMpiSize = comm.Get_size()
iMpiRank = comm.Get_rank()



if iMpiRank == 0:
    print '##########################################'
    print 'modelCase     :', args.iModelCase
    print 'pointing case :', args.iPointingCase
    print 'instrument    :', args.iInstrumentSelector
    print ''
    if args.iPointingCase == spice_:
        print 'SPICE pointing selected:'
        print '   -KernelFile :', args.StringKernelMetaFile.split('/')[-1]
        print '   -Date       :', args.StringUtcStartTime
    elif args.iPointingCase == userPointing_:
        print 'User pointing selected:'
        print '   -R          : %.2e [km]' % args.UserR
        print '   -phase angle: %f [deg]' % args.UserPhaseAngle
        print '   -latitude   : %f [deg]' % args.UserLatitude
        print '   -alpha      : %f [deg]' % args.UserAlpha
        print '   -beta       : %f [deg]' % args.UserBeta
        print '   -gamma      : %f [deg]' % args.UserGamma
    print ''
    if args.iModelCase == dsmc_:
        print "DSMC case selected:"
        print '   -case:      : %s' % args.StringDataFileDSMC.split('/')[-2]
        print '   -species    : %s' % args.StringDataFileDSMC.split('.')[-2]
        if args.IsDust == 1:
            print '   -dust min r : %.3e [m]' %args.DustSizeMin
            print '   -dust max r : %.3e [m]' %args.DustSizeMax
    elif args.iModelCase == haser_:
        print 'HASER case selected:'
        print '   -QHaser     :', QHaser
        print '   -vHaser     :', vHaser
        print '   -tdHaser    :', tdHaser
        print '   -tpHaser    :', tpHaser
    elif args.iModelCase == dust_:
        print 'Analytic DUST case selected:'
        print '   -Qdust      : %f [kg/s]' % QDust 
        print '   -vdust      : %f [m/s]' % vDust
    elif args.iModelCase == userModel_:
        print 'USER coma model uploaded:'
        print '   -filename   : %s' % args.StringUserDataFile
        print '   -delimiter  : "%s"' % args.DelimiterData
        print '   -dimensions : %i' % args.iDimUser
    print '##########################################'


if iModelCase == dsmc_:
    ############################################
    # check if 1d or 2d case
    ############################################

    dataDir = os.path.split(StringDataFileDSMC)[0]
    filenames = os.listdir(dataDir)
    allFilenamesInOneString = ''
    for filename in filenames:
        allFilenamesInOneString += filename

    if '2d' in allFilenamesInOneString:
        iDim = 2
    elif '1d' in allFilenamesInOneString:
        iDim = 1
    elif '3d' in allFilenamesInOneString:
        iDim = 3
    else:
        iDim = 0
        if iModelCase == dsmc_:
            print 'Could not detect number of dimensions of dsmc case. Exiting now.'
        sys.exit()

elif iModelCase in [haser_, dust_]:
    iDim = 1
elif iModelCase == userModel_:
    iDim = args.iDimUser

if iMpiRank == 1:
    print 'iDimensions:', iDim

########################################################
# load data
########################################################
if iModelCase == dsmc_ and iDim < 3:
    if IsDust:
        if iMpiRank == 0:
            print 'dust case'
        NumberDensityIndices, allSizeIntervals = getAllDustIntervalIndices(StringDataFileDSMC, iDim)
        x, y, numberDensities, massDensities = load_dust_data_full(allSizeIntervals, NumberDensityIndices,
                                                                   iDim, StringDataFileDSMC, args)
    else:
        x, y, numberDensities = loadGasData(StringDataFileDSMC, iDim)

elif iModelCase == haser_:
    x, numberDensities = haserModel(QHaser, vHaser, tpHaser, tdHaser)
    y = None

elif iModelCase == dust_:
    x, numberDensities, massDensities, allSizeIntervals = dustModel(QDust, vDust)
    y = None

elif iModelCase == userModel_:
    x, y, numberDensities = load_user_data(StringUserDataFile, iDim, DelimiterData, nHeaderRowsData)

if iDim < 3:
    if numberDensities.ndim == 1:
        numberDensities = np.array([[n] for n in numberDensities])

    nSpecies = numberDensities.shape[1]
    if iMpiRank == 0:
        print 'Nr of species:', nSpecies
elif iDim == 3:
    nSpecies = 1

if args.SubCase == 'preferred':
    print 'preferred case'
    AU = np.float(args.StringDataFileDSMC.split('/')[-1].split('_')[1])
    print "AU:", AU
    if AU >= 3.0:
        numberDensities = numberDensities * 0.24
    else:
        numberDensities = numberDensities * 1.25
        

##############################################################
# triangulation and interpolation for 2d case
if iMpiRank == 0:
    print 'start interpolation'
if iDim == 1:
    pass
elif iDim == 2:
    Triangles = mtri.Triangulation(x, y)
    Interpolator = [mtri.LinearTriInterpolator(Triangles, numberDensities[:, i]) for i in range(nSpecies)]
elif iDim == 3:
    pass

if iMpiRank == 0:
    print 'interpolation done'
#############################################################

#############################################################
# instrument specific definitions
#############################################################

if iInstrumentSelector == osirisw_:       # osiris wac
    nPixelsX = 1024                        # nr of pixels along x axis
    nPixelsY = 1024                        # nr of pixels along y axis
    PhiX = 12 / 2                         # instrument FOV in x (half opening angle) in degrees
    PhiY = 12 / 2                         # instrument FOV in y (half opening angle) in degrees
    iFOV = 0.000993                       # pixel FOV in rad
    PixelSize = 1                         # area of one pixel
    InstrumentFrame = 'ROS_OSIRIS_WAC'

elif iInstrumentSelector == osirisn_:     # osiris nac
    
    nPixelsX = 1024 
    nPixelsY = 1024 
    nPixelsX = 1 
    nPixelsY = 1 
    PhiX = 3 / 2
    PhiY = 3 / 2
    iFOV = 0.0000188
    PixelSize = 1
    InstrumentFrame = 'ROS_OSIRIS_NAC'
    
    '''
    specs = sun.get_specs()
    instrument = Instrument(specs)
    nPixelsX = instrument.nPixelsX
    nPixelsY = instrument.nPixelsY
    PhiX = instrument.PhiX
    PhiY = instrument.PhiY
    InstrumentFrame = instrument.frame
    '''

elif iInstrumentSelector in [alice_, aliceSpec_]:         # alice
    specs = alice.get_specs()
    instrument = Instrument(specs)
    nPixelsX = instrument.nPixelsX
    nPixelsY = instrument.nPixelsY
    PhiX = instrument.PhiX
    PhiY = instrument.PhiY
    InstrumentFrame = instrument.frame

elif iInstrumentSelector in [miro_]:               # miro

    nPixelsX = 1
    nPixelsY = 1
    PhiX = 0.396667 / 2
    PhiY = 0.396667 / 2
    iFOV = 0.00692314
    PixelSize = 1
    InstrumentFrame = 'ROS_MIRO_MM'

elif iInstrumentSelector in [miroDustIR_]:               # miro

    nPixelsX = 1
    nPixelsY = 1
    PhiX = 0.125 / 2
    PhiY = 0.125 / 2
    iFOV = 0.00218166 
    PixelSize = 1
    InstrumentFrame = 'ROS_MIRO_MM'


elif iInstrumentSelector in [miroDustIRSpread_]:               # miro

    nPixelsX = 11
    nPixelsY = 11
    PhiX = 10.0 / 2
    PhiY = 10.0 / 2
    iFOV = 0.00218166
    PixelSize = 1
    InstrumentFrame = 'ROS_MIRO_MM'
    
elif iInstrumentSelector == virtism_:       # virtis m
    nPixelsX = 256
    nPixelsY = 256
    PhiX = 3.6669 / 2
    PhiY = 3.6669 / 2
    iFOV = 0.00025
    PixelSize = 1
    InstrumentFrame = 'ROS_VIRTIS-M'
    
    '''
    # for debugging values
    nPixelsX = 256 
    nPixelsY = 256 
    PhiX = 3.6669
    PhiY = 3.6669
    iFOV = 0.00025
    PixelSize = 1
    '''

elif iInstrumentSelector == virtish_:       # virtis h
    nPixelsX = 1
    nPixelsY = 3
    PhiX = 0.0334
    PhiY = 0.1
    iFOV = 0.000583
    PixelSize = 1
    InstrumentFrame = 'ROS_VIRTIS-H'

Lx = 2 * np.sin(PhiX / 180 * np.pi)
Ly = 2 * np.sin(PhiY / 180 * np.pi)

if nPixelsX > 1:
    Dx = Lx / (nPixelsX - 1)
else:
    Dx = Lx

if nPixelsY > 1:
    Dy = Ly / (nPixelsY - 1)
else:
    Dy = Ly


if iPointingCase == spice_:
    #################################################
    # get rosetta coordinates from spice
    #################################################
    spice.furnsh(StringKernelMetaFile)
    Et = spice.str2et(StringUtcStartTime)
    #R = spice.pxform("ROS_SPACECRAFT", "67P/C-G_CSO", Et)      # create rotation matrix R to go from instrument reference frame to CSO
    if ((iDim == 2) or (iDim == 1)):
        rRosetta, lightTime = spice.spkpos("ROSETTA", Et, "67P/C-G_CSO", "NONE", "CHURYUMOV-GERASIMENKO")        # s/c coordinates in CSO frame of reference
        R = spice.pxform(InstrumentFrame, "67P/C-G_CSO", Et)      # create rotation matrix R to go from instrument reference frame to CSO
    elif iDim == 3:
        print "iDim = 3"
        rRosetta, lightTime = spice.spkpos("ROSETTA", Et, "67P/C-G_CK", "NONE", "CHURYUMOV-GERASIMENKO")        # s/c coordinates in CSO frame of reference
        R = spice.pxform(InstrumentFrame, "67P/C-G_CK", Et)
        # do 3d stuff from here?
        rSun, lt = spice.spkpos("SUN", Et, "67P/C-G_CK", "NONE", "CHURYUMOV-GERASIMENKO")
        rSun = np.array(rSun)
        rSun = rSun / np.sqrt((rSun**2).sum())
        # quick and dirty fix --> andre  needs to define instrument object for all instruments!
        try:
            if instrument.name == 'sun':
                rRosetta = rSun * 150  # set distance from sun at 200 km because of domain size
                zhat = np.array([0,0,1])
                R = rotations.rotMat(zhat, -rSun)
        except:
            pass
        if iMpiRank == 0:
          with open('rSun_hat.dat', 'w') as sFile:
            sFile.write("%.5e,%.5e,%.5e\n" %(rSun[0], rSun[1], rSun[2]))
        print "R:", R
        try:
            print "-rSun", -rSun
            print "rRosetta", rRosetta
            print "pnadir:", dot(R,zhat)
        except:
            pass
	  
    rRosetta = np.array(rRosetta) * 1000            # transform km to m
elif iPointingCase == userPointing_:
    x0 = np.array([-UserR*1000, 0, 0])           # -UserR --> start at subsolar point, in meters
    rRosetta = rotations.rotateVector(x0, UserPhaseAngle, UserLatitude)
    ei, ej, ek = rotations.rotateCoordinateSystem2(UserPhaseAngle, UserLatitude, UserAlpha, UserBeta, UserGamma)
    R = rotations.createRotationMatrix(ei, ej, ek)

    # do 3d stuff from here?
    rSun = np.array([0, 1, 0]) 
    if iMpiRank == 0:
      with open('rSun_hat.dat', 'w') as sFile:
        sFile.write("%.5e,%.5e,%.5e\n" %(rSun[0], rSun[1], rSun[2]))

if iMpiRank == 0:
    print 'Distance from comet  : %.2e [m]' % (np.sqrt(np.sum(rRosetta ** 2)))
    if iDim < 3:
        print 'rRosetta in CSO Frame: (%.2e, %.2e, %.2e)' % (rRosetta[0], rRosetta[1], rRosetta[2])
    else:
        print 'rRosetta in CK Frame: (%.2e, %.2e, %.2e)' % (rRosetta[0], rRosetta[1], rRosetta[2])


##############################################
# create pointing vectors p
##############################################
i = np.arange(nPixelsX)
j = np.arange(nPixelsY)
ii, jj = np.meshgrid(i, j, indexing='ij')
p = np.array([np.zeros((len(i), len(j))) for k in range(3)])

if args.iPointingCase == spice_:
    p[0] = ii*Dx - Lx/2 + Dx/2
    p[1] = jj*Dy - Ly/2 + Dy/2
    p[2] = np.ones((len(i), len(j)))
    p_hat = p / np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)

elif args.iPointingCase == userPointing_:
    p[1] = ii*Dx - Lx/2 + Dx/2
    p[2] = jj*Dy - Ly/2 + Dy/2
    p[0] = np.ones((len(i), len(j)))
    p_hat = p / np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)

ccd = np.zeros((nPixelsX, nPixelsY, nSpecies))
if iDim < 3:
    cso2tenishev = np.array([-1, -1, 1])
else:
    cso2tenishev = np.array([1, 1, 1])
kkk = 0
nnn = 0

percentProgressLast = 0


if iMpiRank == 0:
    if iDim == 3:
        pFile = open('pointing.dat', 'w')
        
    print ''
    print 'Entering pixel loop.  Progress ...'
for i in range(nPixelsX):
    for j in range(nPixelsY):
        if iDim < 3:
            if (kkk == (iMpiRank + nnn * nMpiSize)):
                if iPointingCase == spice_ and iDim < 3:
                    p = np.dot(R, p_hat[:, i, j]) * cso2tenishev
                    rRay = np.array([value for value in rRosetta]) * cso2tenishev
                else:
                    p = np.dot(R, p_hat[:, i, j])
                    rRay = np.array([value for value in rRosetta])
                xTravel = np.array(createRay.createRay(rRay, p))
                dTravel = np.sqrt(np.sum((xTravel[0] - xTravel)**2, axis=1))

                if iDim == 1:
                    xTravel = np.sqrt(np.sum(xTravel ** 2, axis=1))
                elif iDim == 2:
                    xTravel[:, 1] = np.sqrt(xTravel[:, 1]**2 + xTravel[:, 2]**2)
                    #print "xTravel:", xTravel
                    #print "dTravel:", dTravel
                elif iDim == 3:
                    pass

                # loop over species
                for spIndex in range(nSpecies):
                    if iDim == 1:
                        DensityRay = np.interp(xTravel, x, numberDensities[:, spIndex])
                    elif iDim == 2:
                        DensityRay = Interpolator[spIndex].__call__(xTravel[:, 0], xTravel[:, 1]) #interpolated local number density
                    elif iDim == 3:
                        DensityRay = None

                    ColumnDensity = np.trapz(DensityRay, dTravel)
                    if iMpiRank != 0:
                        data = np.array([ColumnDensity, i, j, spIndex])
                        comm.send([ColumnDensity, i, j, spIndex], dest=0, tag=13)
                    else:
                        ccd[i][j][spIndex] = ColumnDensity
                        for jjj in range(1, nMpiSize):
                            ColumnDensity, ii, jj, spIndex = comm.recv(source=jjj, tag=13)
                            ccd[ii][jj][spIndex] = ColumnDensity
                nnn += 1
            kkk += 1
        elif iDim == 3 and iMpiRank == 0:
            if iPointingCase == spice_:
                p = np.dot(R, p_hat[:, i, j])
                rRay = np.array([value for value in rRosetta])
                pFile.write("%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n" %(p[0], p[1], p[2], rRay[0], rRay[1], rRay[2]))
            else:
                #print "Only spice pointing allowed for 3D cases so far!"
                p = np.dot(R, p_hat[:, i, j])
                rRay = np.array([value for value in rRosetta])
                pFile.write("%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n" %(p[0], p[1], p[2], rRay[0], rRay[1], rRay[2]))
                #sys.exit()
    
           

    if iMpiRank == 0:
        percentProgress = np.floor(i / nPixelsX * 10) * 10
        if percentProgress > percentProgressLast:
             percentProgressLast = percentProgress
             print int(percentProgress),'%'
if iDim == 3:
    pFile.close()

    print args.StringOutputDir

    if args.iIlluminationCase == 0:
        if "H2O" not in args.StringDataFileDSMC:
            parts = args.StringDataFileDSMC.split(".")
            args.StringDataFileDSMC = parts[0] + "." + parts[1] + "." + parts[2] +".H2O.dat"
        os.system("su _www -c '/Applications/Julia-0.3.0.app/Contents/Resources/julia/bin/julia /Users/abieler/newLOS/newLOS.jl %s %s'" %(args.StringDataFileDSMC, args.StringOutputDir))
    elif args.iIlluminationCase == 1:
        os.system("su _www -c '/Applications/Julia-0.3.0.app/Contents/Resources/julia/bin/julia /Users/abieler/newLOS/illumination.jl %s'" %(args.StringOutputDir))
    ccdLoaded = np.loadtxt("ccd.dat")

    # for miro ccdLoaded is a 0 sized tuple which cannot be iterated
    # hence we have to replace that with a 1 element array
    if ccdLoaded.size == 1:
        ccdLoaded = np.array([float(ccdLoaded)])
    k = 0
    for i in range(nPixelsX):
      for j in range(nPixelsY):
        ccd[i,j,0] = ccdLoaded[k] 
        k += 1
if iMpiRank == 0:
    print 'pixel loop done'
    if iInstrumentSelector in [miroDustIR_, miroDustIRSpread_]:
        #ccd, aveBright, frequencies = miro.fluxDensity(ccd, allSizeIntervals, iFOV, args)
        ccdTemp, aveBright, frequencies = miro.fluxDensity(ccd, allSizeIntervals, iFOV, args)
        nSpecies = ccdTemp.shape[2]                       
        ccd = np.zeros((nPixelsX, nPixelsY, nSpecies))
        ccd = ccdTemp

    ccd_limits = (ccd.min(), ccd.max())

    ######################################################
    # plot results
    #######################################################
    if args.iInstrumentSelector in [miroDustIR_, miroDustIRSpread_]:
        plot_miro(ccd, aveBright, frequencies, allSizeIntervals, args)

    if args.IsDust:
        UserDustSizes= [size for size in allSizeIntervals
                       if (args.DustSizeMin <= size <= args.DustSizeMax)]

    for spIndex in range(nSpecies):
        if iInstrumentSelector in [alice_, aliceSpec_]:
            ccdFinal, wavelengths = alice.calculateBrightness(instrument.nOversampleX,
                                                              instrument.nOversampleY,
                                                              instrument.PixelFOV,
                                                              ccd[:, :, spIndex],
                                                              args)
        else:
            ccdFinal = ccd[:, :, spIndex]

        ######################################################
        # write results to file
        ######################################################
        pltTitle = build_plot_title(args, 'LOS')
        filename = '/result_%i.txt' % (spIndex)
        with open(StringOutputDir + filename, 'w') as f:
            f.write(pltTitle)
            f.write('\n')
            if args.IsDust:
                f.write('Dust size: %e [m]\n' % UserDustSizes[spIndex])
            if args.iInstrumentSelector not in [alice_, aliceSpec_]:
                f.write("x [pixel], y [pixel], %s \n" % ccdStr[args.iInstrumentSelector])
            f.write("/begin data\n")
            if iInstrumentSelector in [alice_, aliceSpec_]:
                alice.save_results(f, ccdFinal, wavelengths, filename)
            else:
                ix = 0
                iy = 0
                for row in ccdFinal:
                    for value in row:
                        f.write('%5i, %5i, %.5e\n' % (ix, iy, value))
                        iy += 1
                    ix += 1
                    iy = 0

        ######################################################
        # plot results
        #######################################################
        if args.StringPlotting.lower() == 'matplotlib':
            plotName = 'result_%i.png' % (spIndex)
        elif args.StringPlotting.lower() == 'bokeh':
            plotName = 'result_%i.html' % (spIndex)

        plot_result_LOS(ccdFinal, plotName, args, ccd_limits)

if iMpiRank == 0:
    print '**' * 20
    print 'Run completed'
    print 'Time elapsed: %.2f seconds' % (time.time() - startTime)
