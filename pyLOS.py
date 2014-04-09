'''
filename: pyLOS.py

requires pySPICE installed
requires matplotlib version 1.3.1 or higher
requires mpi4py installed
'''

from __future__ import division                 # must be first line of program

try:
    import sys
    print sys.path
    import os
    import numpy as np
    import time
    import argparse
    import matplotlib.tri as mtri
    import matplotlib.pyplot as plt
    import spice
    #import mpi4py

    from data_loaders import loadGasData, loadDustData, getAllDustIntervalIndices
    from haser import haserModel
    from data_plotting import plot_result
    import rotations
    import alice

except Exception,e:
    print '--' * 20
    print "Error with import of python module:"
    print e
    print '--' * 20
    sys.exit()

useHaserModel = False

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

parser.add_argument("--UserR", type=float)                               # distance in km from nucleus center
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
userNrOfGeaderRows = args.iUserNrOfHeaderRows


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
        print 'Could not detect number of iDimensions of dsmc case. Exiting now.'
        sys.exit()

elif iModelCase == 1:
    iDim = 1

elif iModelCase == 2:
    iDim = iUserDim

print 'iDimensions:', iDim

if iPointingCase == 0:
    #################################################
    # get rosetta coordinates from spice
    #################################################
    spice.furnsh(StringKernelMetaFile)
    et = spice.str2et(StringUtcStartTime)
    rRosetta, lightTime = spice.spkpos("ROSETTA", et, "67P/C-G_CSO", "NONE", "CHURYUMOV-GERASIMENKO")        # s/c coordinates in CSO frame of reference
    rRosetta = np.array(rRosetta) * 1000            # transform km to m
    R = spice.pxform("ROS_SPACECRAFT", "67P/C-G_CSO", et)      # create rotation matrix R to go from instrument reference frame to CSO

    print 'distance from comet:', np.sqrt(np.sum(rRosetta ** 2))

elif iPointingCase == 1:
    x0 = np.array([-UserR, 0, 0])           # -R --> start at subsolar point
    rRosetta = rotations.rotateVector(x0, UserPhaseAngle, UserLatitude)
    ei, ej, ek = rotations.rotateCoordinateSystem2(UserPhaseAngle, UserLatitude, UserAlpha, UserBeta, UserGamma)

    R = rotations.createRotationMatrix(ei,ej,ek)

print 'rRosetta:', rRosetta

########################################################
# load data
########################################################
if iModelCase == 0:
    if IsDust:
        print 'dust case'
        numberDensityIndices, allSizeIntervals = getAllDustIntervalIndices(StringDataFileDSMC, iDim)

        x, y, n = loadDustData(allSizeIntervals, numberDensityIndices, iDim, StringDataFileDSMC)

    else:
        x, y, n = loadGasData(StringDataFileDSMC, iDim)

elif iModelCase == 1:
    x, n = haserModel(QHaser, vHaser, tpHaser, tdHaser)
    y = None

elif iModelCase == 2:
    x, y, n = loadGasData(StringUserDataFile, iDim, True, UserDelimiter, userNrOfGeaderRows)

##############################################################
# triangulation and interpolation for 2d case
if iDim == 1:
    pass
elif iDim == 2:
    triangles = mtri.Triangulation(x, y)
    Interpolator = mtri.LinearTriInterpolator(triangles, n)
print 'interpolation done'
#############################################################

#############################################################
# instrument specific definitions
#############################################################

if iInstrumentSelector == 1:                 # osiris wac
    pixelsX = 256                           # nr of pixels along x axis
    pixelsY = 256                           # nr of pixels along y axis
    phi_x = 12 / 2                            # instrument FOV in x (half opening angle) in degrees
    phi_y = 12 / 2                            # instrument FOV in y (half opening angle) in degrees
    iFOV = 0.000993                         # pixel FOV in rad
    pixelSize = 1                           # area of one pixel

elif iInstrumentSelector == 2:               # osiris nac
    pixelsX = 512
    pixelsY = 512
    phi_x = 3 / 2
    phi_y = 3 / 2
    iFOV = 0.0000188
    pixelSize = 1

elif iInstrumentSelector == 3:               # alice
    N_oversampleX = 24
    N_oversampleY = 20

    pixelsX = 19 * N_oversampleX
    pixelsY = 1 * N_oversampleY
    phi_x = 5.852 / 2
    phi_y = 0.1 / 2

    iFOV = (phi_x * 2 / 180 * np.pi / (19 * N_oversampleX)) * (phi_y * 2 / 180 * np.pi / (N_oversampleY))

    print iFOV
    pixelSize = 1

    #v_sun = alice.get_v_sun(StringKernelMetaFile, StringUtcStartTime)
    #gFactor = alice.get_gfactor_from_db()

    v_sun = 12
    gFactor = 2.09e-7

elif iInstrumentSelector == 4:               # miro
    pixelsX = 1
    pixelsY = 1
    phi_x = 0.33336 / 2
    phi_y = 0.36666 / 2
    iFOV = 0.36666
    pixelSize = 1

elif iInstrumentSelector == 5:               # virtis m
    pixelsX = 256
    pixelsY = 256
    phi_x = 3.6669 / 2
    phi_y = 3.6669 / 2
    iFOV = 0.00025
    pixelSize = 1

elif iInstrumentSelector == 6:               # virtis h
    pixelsX = 1
    pixelsY = 3
    phi_x = 0.0334
    phi_y = 0.1
    iFOV = 0.000583
    pixelSize = 1

lx = 2 * np.sin(phi_x / 180 * np.pi)
ly = 2 * np.sin(phi_y / 180 * np.pi)

if pixelsX > 1:
    dx = lx / (pixelsX - 1)
else:
    dx = lx

if pixelsY > 1:
    dy = ly / (pixelsY - 1)
else:
    dy = ly

print 'entering pixel loop'

p = np.zeros(3)
ccd = np.zeros((pixelsX, pixelsY))
for i in range(pixelsX):
    for j in range(pixelsY):
        if iPointingCase == 0:
            p[0] = i*dx - lx/2 + dx/2
            p[1] = j*dy - ly/2 + dy/2
            p[2] = 1
            p /= np.sqrt(np.sum(p**2))          # p = pointing vector in instrument coordinates
            pSpice = spice.mxv(R, p)            # pointing vector in CSO frame of reference (sun on +x axis)
            p[0] = -pSpice[0]                   # p = pointing vector in simulation coordinates (sun on -x axis)
            p[1] = -pSpice[1]
            p[2] = pSpice[2]

            rRay = np.array([value for value in rRosetta])  # ray for line of sight, starting at rosetta position (transform to coords with sun on -x)
            rRay[0] *= -1
            rRay[1] *= -1
        else:
            p[1] = i*dx - lx/2 + dx/2
            p[2] = j*dy - ly/2 + dy/2
            p[0] = 1
            p /= np.sqrt(np.sum(p**2))          # p = pointing vector in instrument coordinates
            p = R.dot(p)
            rRay = np.array([value for value in rRosetta])

        #########################################
        # calculate ray for line of sight
        #########################################
        distance = np.sqrt(np.sum(rRay**2))
        xTravel = []
        dTravel = []                                           # distance from s/c
        columnDensity = 0
        drTotal = 0
        while (distance <= 1e8 and distance >= 2e3):            # adapt dr according to distance from the nucleus, the closer to it, the smaller dr becomes.
            distance = np.sqrt(np.sum(rRay**2))
            dTravel.append(drTotal)

            if iDim == 1:
                xTravel.append(distance)
            elif iDim == 2:
                xTravel.append((rRay[0], np.sqrt(np.sum(rRay[1]**2 + rRay[2] ** 2))))

            if distance < 1e4:
                dr = distance / 40
                if distance < 4e3:
                    dr = distance / 100
                    if distance < 2.5e3:
                        dr = distance / 250
                        if distance < 2030:
                            dr = 1

            else:
                dr = distance / 10
            rRay = rRay + p * dr
            drTotal += dr

        xTravel = np.array(xTravel)

        if iDim == 1:
            nRay = np.interp(xTravel, x, n)
        elif iDim == 2:
            nRay = Interpolator.__call__(xTravel[:, 0], xTravel[:, 1])        # interpolated local number density
        elif iDim == 3:
            pass

        columnDensity = np.trapz(nRay, dTravel)
        ccd[i][j] = columnDensity

    print i
print 'pixel loop done'

ccd = np.array(ccd)

if iInstrumentSelector == 3:
    ccdFinal = alice.calculateBrightness(N_oversampleX, N_oversampleY, ccd, gFactor)
else:
    ccdFinal = ccd

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
plot_result(ccdFinal, StringOutputDir, 'result.png', iInstrumentSelector, RunDetails=args, showPlot=True)
