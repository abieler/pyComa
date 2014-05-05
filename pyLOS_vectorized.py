#!/opt/local/anaconda/bin/python
'''
filename: pyLOS_vectorized.py

requires pySPICE installed

requires matplotlib version 1.3.1 or higher

***********************************************************************************************************
requires mpi4py installed (on osx: anaconda mpi4py package is broken for version 1.9.2, create a symlink:
sudo ln -s /yourPathTo/anaconda /opt/anaconda1anaconda2anaconda3 solves the problem
***********************************************************************************************************

***************************************************
requires cython installed.
--> execute:

python setupOSX.py build_ext --inplace

in the pyComa directory, this compiles the cython
module to build createRay.so which is then imported
by pyLOS_vectorized.py
****************************************************

'''
from __future__ import division                 # must be first line of program

try:
    # import standard modules
    import sys
    import os
    import time

    # import 3rd party modules
    import spice
    import numpy as np
    import matplotlib.tri as mtri
    import matplotlib.pyplot as plt
    from mpi4py import MPI
    import argparse

    # import own modules
    from utils.data_loaders import *
    from utils.haser import haserModel
    from utils.data_plotting import plot_result_LOS, build_plot_title
    import utils.rotations as rotations
    import utils.alice as alice
    import utils.createRay as createRay
    from cmdline_args import cmdline_args

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
    print 'modelCase     :', iModelCase
    print 'pointing case :', iPointingCase
    print 'instrument    :', iInstrumentSelector
    print 'StringKernelMetaFile:', StringKernelMetaFile
    print 'StringUtcStartTime  :', StringUtcStartTime
    print 'QHaser        :', QHaser
    print 'vHaser        :', vHaser
    print 'tdHaser       :', tdHaser
    print 'tpHaser       :', tpHaser
    print '##########################################'


if iModelCase == 0:
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
    else:
        iDim = 0
        if iModelCase == 0:
            print 'Could not detect number of iDimensions of dsmc case. Exiting now.'
        sys.exit()

elif iModelCase == 1:
    iDim = 1
elif iModelCase == 2:
    iDim = iDimUser

if iMpiRank == 1:
    print 'iDimensions:', iDim

if iPointingCase == 0:
    #################################################
    # get rosetta coordinates from spice
    #################################################
    spice.furnsh(StringKernelMetaFile)
    Et = spice.str2et(StringUtcStartTime)
    rRosetta, lightTime = spice.spkpos("ROSETTA", Et, "67P/C-G_CSO", "NONE", "CHURYUMOV-GERASIMENKO")        # s/c coordinates in CSO frame of reference
    rRosetta = np.array(rRosetta) * 1000            # transform km to m
    R = spice.pxform("ROS_SPACECRAFT", "67P/C-G_CSO", Et)      # create rotation matrix R to go from instrument reference frame to CSO
    if iMpiRank == 0:
        print 'Distance from comet: %.2e' % (np.sqrt(np.sum(rRosetta ** 2)))

elif iPointingCase == 1:
    x0 = np.array([-UserR, 0, 0])           # -R --> start at subsolar point
    rRosetta = rotations.rotateVector(x0, UserPhaseAngle, UserLatitude)
    ei, ej, ek = rotations.rotateCoordinateSystem2(UserPhaseAngle, UserLatitude, UserAlpha, UserBeta, UserGamma)

    R = rotations.createRotationMatrix(ei, ej, ek)

if iMpiRank == 0:
    print 'Distance from comet: %.2e' % (np.sqrt(np.sum(rRosetta ** 2)))
    print 'rRosetta: %.2e, %.2e, %.2e' % (rRosetta[0], rRosetta[1], rRosetta[2])

########################################################
# load data
########################################################
if iModelCase == 0:
    if IsDust:
        if iMpiRank == 0:
            print 'dust case'
        NumberDensityIndices, allSizeIntervals = getAllDustIntervalIndices(StringDataFileDSMC, iDim)
        x, y, numberDensities, massDensities = load_dust_data_miro(allSizeIntervals, NumberDensityIndices, iDim, StringDataFileDSMC, args)
    else:
        x, y, numberDensities = loadGasData(StringDataFileDSMC, iDim)

elif iModelCase == 1:
    x, numberDensities = haserModel(QHaser, vHaser, tpHaser, tdHaser)
    y = None

elif iModelCase == 2:
    #x, y, n = loadGasData(StringUserDataFile, iDim, True, DelimiterData, nHeaderRowsData)
    x, y, numberDensities = load_user_data(StringUserDataFile, iDim, DelimiterData, nHeaderRowsData)


if numberDensities.ndim == 1:
    numberDensities = np.array([[n] for n in numberDensities])

nSpecies = numberDensities.shape[1]
nSpecies = 4
print '+++++++++++++++++++++++++++++'
print numberDensities[:, 0]
print '+++++++++++++++++++++++++++++'

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

if iInstrumentSelector == 1:                 # osiris wac
    nPixelsX = 256                          # nr of pixels along x axis
    nPixelsY = 256                          # nr of pixels along y axis
    PhiX = 12 / 2                            # instrument FOV in x (half opening angle) in degrees
    PhiY = 12 / 2                            # instrument FOV in y (half opening angle) in degrees
    iFOV = 0.000993                         # pixel FOV in rad
    PixelSize = 1                           # area of one pixel

elif iInstrumentSelector == 2:               # osiris nac
    nPixelsX = 512
    nPixelsY = 512
    PhiX = 3 / 2
    PhiY = 3 / 2
    iFOV = 0.0000188
    PixelSize = 1

elif iInstrumentSelector == 3:               # alice
    nOversampleX = 24
    nOversampleY = 20

    nPixelsX = 19 * nOversampleX
    nPixelsY = 1 * nOversampleY
    PhiX = 5.852 / 2
    PhiY = 0.1 / 2

    iFOV = (PhiX * 2 / 180 * np.pi / (19 * nOversampleX)) * (PhiY * 2 / 180 * np.pi / (nOversampleY))

    print iFOV
    PixelSize = 1

    #v_sun = alice.get_v_sun(StringKernelMetaFile, StringUtcStartTime)
    #gFactor = alice.get_gfactor_from_db()

    v_sun = 12
    gFactor = 2.09e-7

elif iInstrumentSelector == 4:               # miro
    nPixelsX = 1
    nPixelsY = 1
    PhiX = 0.33336 / 2
    PhiY = 0.36666 / 2
    iFOV = 0.36666
    PixelSize = 1

elif iInstrumentSelector == 5:               # virtis m
    nPixelsX = 256
    nPixelsY = 256
    PhiX = 3.6669 / 2
    PhiY = 3.6669 / 2
    iFOV = 0.00025
    PixelSize = 1

elif iInstrumentSelector == 6:               # virtis h
    nPixelsX = 1
    nPixelsY = 3
    PhiX = 0.0334
    PhiY = 0.1
    iFOV = 0.000583
    PixelSize = 1

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

##############################################
# create pointing vectors p
##############################################
i = np.arange(nPixelsX)
j = np.arange(nPixelsY)
ii, jj = np.meshgrid(i, j, indexing='ij')
p = np.array([np.zeros((len(i), len(j))) for k in range(3)])

if args.iPointingCase == 0:
    p[0] = ii*Dx - Lx/2 + Dx/2
    p[1] = jj*Dy - Ly/2 + Dy/2
    p[2] = np.ones((len(i), len(j)))
    p_hat = p / np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)

elif args.iPointingCase == 1:
    p[1] = ii*Dx - Lx/2 + Dx/2
    p[2] = jj*Dy - Ly/2 + Dy/2
    p[0] = np.ones((len(i), len(j)))
    p_hat = p / np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)

ccd = np.zeros((nPixelsX, nPixelsY, nSpecies))
ccdFinal = np.zeros((nPixelsX, nPixelsY, nSpecies))
if iInstrumentSelector == 3:
    ccdFinal = np.zeros((19, nSpecies))

cso2tenishev = np.array([-1, -1, 1])

kkk = 0
nnn = 0

if iMpiRank == 0:
    print 'entering pixel loop'
for i in range(nPixelsX):
    for j in range(nPixelsY):
        if (kkk == (iMpiRank + nnn * nMpiSize)):
            if iPointingCase == 0:
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
            elif iDim == 3:
                pass

            # loop over species
            for spIndex in range(nSpecies):
                if iDim == 1:
                    DensityRay = np.interp(xTravel, x, numberDensities[:, spIndex])
                elif iDim == 2:
                    DensityRay = Interpolator[spIndex].__call__(xTravel[:, 0], xTravel[:, 1])        # interpolated local number density
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

    if iMpiRank == 0:
        print i

if iMpiRank == 0:
    print 'pixel loop done'
    ccd_limits = (ccd.min(), ccd.max())
    print 'max ccd:', ccd_limits[1]
    print 'min ccd:', ccd_limits[0]

    for spIndex in range(nSpecies):
        if iInstrumentSelector == 3:
            ccdFinal[:, spIndex] = alice.calculateBrightness(nOversampleX, nOversampleY, ccd[:, :, spIndex], gFactor)
        else:
            ccdFinal[:, :, spIndex] = ccd[:, :, spIndex]

        ######################################################
        # write results to file
        ######################################################
        pltTitle = build_plot_title(args, 'LOS')
        filename = '/result_%i.txt' % (spIndex)
        with open(StringOutputDir + filename, 'w') as f:
            f.write(pltTitle)
            f.write('\n')
            f.write("Each datapoint is the column number density in 1/m2 for an instrument  pixel.\n")
            f.write("Rows correspond to pixels in instruments X axis, starting with the most negative value.\n")
            f.write("Columns correspond to pixels in instrument Y axis, starting with the most negative value.\n")
            f.write("/begin data\n")
            if iInstrumentSelector == 3:
                for value in ccdFinal[:, spIndex]:
                    f.write('%e\n' % value)
            else:
                for row in ccdFinal[:, :, spIndex]:
                    for value in row:
                        f.write('%e,' % value)
                    f.write('\n')

        ######################################################
        # plot results
        #######################################################
        if args.StringPlotting.lower() == 'matplotlib':
            plotName = 'result_%i.png' % (spIndex)
        elif args.StringPlotting.lower() == 'bokeh':
            plotName = 'result_%i.html' % (spIndex)

        if iInstrumentSelector == 3:
            plot_result_LOS(ccdFinal[:, spIndex], plotName, args, ccd_limits)
        else:
            plot_result_LOS(ccdFinal[:, :, spIndex], plotName, args, ccd_limits)
if iMpiRank == 0:
    print '**' * 20
    print 'Time elapsed: %.2f seconds' % (time.time() - startTime)
    print '**' * 20