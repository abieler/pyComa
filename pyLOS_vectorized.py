#!/Users/abieler/anaconda/bin/python
'''
filename: pyLOS.py

requires pySPICE installed

requires matplotlib version 1.3.1 or higher

requires mpi4py installed (on osx: anaconda mpi4py package is broken for version 1.9.2, create a symlink:
sudo ln -s /Users/yourPathTo/anaconda /opt/anaconda1anaconda2anaconda3 solves the problem
'''

from __future__ import division                 # must be first line of program

try:
    import sys
    import os
    import numpy as np
    import time
    import argparse
    import matplotlib.tri as mtri
    import matplotlib.pyplot as plt
    import spice

    from data_loaders import * #loadGasData, loadDustData, getAllDustIntervalIndices, createRay
    from haser import haserModel
    from data_plotting import plot_result
    import rotations
    import alice
    import createRay

except Exception, e:
    print '--' * 20
    print "Error with import of python module:"
    print e
    print '--' * 20
    sys.exit()

try:
    from mpi4py import MPI
except:
    'mpi4py not installed!'
    sys.exit()

startTime = time.time()
#useHaserModel = False

#############################################
# setup argparser for cmd line arguments
##############################################

parser = argparse.ArgumentParser()
parser.add_argument("--iModelCase", type=int, choices=[0, 1, 2], help='0: dsmc model, 1: haser model, 2: user model')
parser.add_argument("--iPointingCase", type=int, choices=[0, 1], help='0: spice pointing, 1: user pointing')
parser.add_argument("--iInstrumentSelector", type=int, choices=[1, 2, 3, 4, 5, 6])
parser.add_argument("--StringOutputDir", type=str)

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

parser.add_argument("--UserR", type=float)                               # Distance in km from nucleus center
parser.add_argument("--UserPhaseAngle", type=float)
parser.add_argument("--UserLatitude", type=float)
parser.add_argument("--UserAlpha", type=float)
parser.add_argument("--UserBeta", type=float)
parser.add_argument("--UserGamma", type=float)

args = parser.parse_args()

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
UserDelimiter = args.UserDelimiter
iUserDim = args.iUserDim
iUserNrOfHeaderRows = args.iUserNrOfHeaderRows

comm = MPI.COMM_WORLD
nMpiSize = comm.Get_size()
iMpiRank = comm.Get_rank()

if iMpiRank == 0:
    print '##########################################'
    print 'modelCase     :', iModelCase
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
    iDim = iUserDim
if iMpiRank == 1:
    print 'iDimensions:', iDim

if iPointingCase == 0:
    #################################################
    # get rosetta coordinates from spice
    #################################################
    spice.furnsh(StringKernelMetaFile)
    Et = spice.str2et(StringUtcStartTime)
    rRosetta, lightTime = spice.spkpos("ROSETTA", Et, "67P/C-G_CSO", "NONE", "CHURYUMOV-GERASIMENKO")        # s/c coordinates in CSO frame of reference
    #rEarth, lightTime = spice.spkpos("EARTH", Et, "J2000", "NONE", "ROSETTA")        # s/c coordinates in CSO frame of reference
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

        x, y, n = loadDustData(allSizeIntervals, NumberDensityIndices, iDim, StringDataFileDSMC)

    else:
        x, y, n = loadGasData(StringDataFileDSMC, iDim)

elif iModelCase == 1:
    x, n = haserModel(QHaser, vHaser, tpHaser, tdHaser)
    y = None

elif iModelCase == 2:
    #x, y, n = loadGasData(StringUserDataFile, iDim, True, UserDelimiter, iUserNrOfHeaderRows)
    x, y, n = load_user_data(StringUserDataFile, iDim, UserDelimiter, iUserNrOfHeaderRows)

##############################################################
# triangulation and interpolation for 2d case
if iDim == 1:
    pass
elif iDim == 2:
    Triangles = mtri.Triangulation(x, y)
    Interpolator = mtri.LinearTriInterpolator(Triangles, n)

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

if iMpiRank == 0:
    print 'entering pixel loop'

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
    p[1] =  ii*Dx - Lx/2 + Dx/2
    p[2] = jj*Dy - Ly/2 + Dy/2
    p[0] = np.ones((len(i), len(j)))
    p_hat = p / np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)

ccd = np.zeros((nPixelsX, nPixelsY))
cso2tenishev = np.array([-1, -1, 1])

kkk = 0
nnn = 0
for i in range(nPixelsX):
    for j in range(nPixelsY):
        if (kkk == (iMpiRank + nnn * nMpiSize)):

            if iPointingCase == 0:
                p = np.dot(R, p_hat[:, i, j]) * cso2tenishev
                rRay = np.array([value for value in rRosetta]) * cso2tenishev
            else:
                p = np.dot(R, p_hat[:, i, j])
                rRay = np.array([value for value in rRosetta])

            xTravel = np.array(createRay.createRay(iDim, rRay, p))
            dTravel = np.sqrt(np.sum((xTravel[0] - xTravel)**2, axis=1))

            if iDim == 1:
                xTravel = np.sqrt(np.sum(xTravel ** 2, axis=1))
            elif iDim == 2:
                xTravel[:, 1] = np.sqrt(xTravel[:, 1]**2 + xTravel[:, 2]**2)

            if iDim == 1:
                DensityRay = np.interp(xTravel, x, n)
            elif iDim == 2:
                DensityRay = Interpolator.__call__(xTravel[:, 0], xTravel[:, 1])        # interpolated local number density
            elif iDim == 3:
                pass
            ColumnDensity = np.trapz(DensityRay, dTravel)

            if iMpiRank != 0:
                data = np.array([ColumnDensity,i,j])
                comm.send([ColumnDensity, i, j], dest=0, tag=13)
            else:
                ccd[i][j] = ColumnDensity
                for jjj in range(1, nMpiSize):
                    ColumnDensity, ii, jj = comm.recv(source=jjj, tag=13)
                    ccd[ii][jj] = ColumnDensity
            nnn += 1
        kkk += 1

    if iMpiRank == 0:
        print i

if iMpiRank == 0:
    print 'pixel loop done'

    ccd = np.array(ccd)

    if iInstrumentSelector == 3:
        ccdFinal = alice.calculateBrightness(nOversampleX, nOversampleY, ccd, gFactor)
    else:
        ccdFinal = ccd

    print '**' * 20
    print '**' * 20
    print time.time() - startTime
    print '**' * 20
    print '**' * 20
    ######################################################
    # write results to file
    ######################################################
    f = open(StringOutputDir + '/result.txt', 'w')
    for row in ccd:
        for value in row:
            f.write('%e,' % value)
        f.write('\n')
    f.close()
    
    ######################################################
    # plot results
    #######################################################
    plot_result(ccdFinal, StringOutputDir, 'result.png', iInstrumentSelector, RunDetails=args, DoShowPlot=True)
