#!/Users/abieler/Library/Enthought/Canopy_64bit/User/bin/python
'''
filename: pyLOS.py

performs line of sight integration on 2d amps output considering spice kernels for
the rosetta spacecraft. includes 6 different instruments and their field of views.

requires pySPICE installed

ToDo:
-----
* ALICE slit geometry
* MIRO dust implementation
* test integration dr if sufficient

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
    
    from dataLoaders import loadGasData, loadDustData, getAllDustIntervalIndices
    from haser import haserModel
    from dataPlotting import plotResult
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
parser.add_argument("--modelCase", type=int, choices=[0,1,2],help='0: dsmc model, 1: haser model, 2: user model')
parser.add_argument("--pointingCase", type=int, choices=[0,1],help='0: spice pointing, 1: user pointing')
parser.add_argument("--instrumentSelector",type=int,choices=[1,2,3,4,5,6])
parser.add_argument("--outputDir",type=str)

parser.add_argument("--dataFile", type=str)
parser.add_argument("--isDust",type=int,choices=[0,1],help='1 for dust case, 0 for gas case')
parser.add_argument("--minSize",type=float)
parser.add_argument('--maxSize',type=float)


parser.add_argument("--QHaser",type=float)
parser.add_argument("--vHaser",type=float)
parser.add_argument("--tpHaser",type=float)
parser.add_argument("--tdHaser",type=float)

parser.add_argument("--userDataFile",type=str)                          # file to upload from user which contains user coma model
parser.add_argument("--userDelimiter",type=str)                         # delimiter used in datafile
parser.add_argument("--userNrOfHeaderRows",type=int)                   # number of header lines in datafile
parser.add_argument("--userDim",type=int)                               # number of dimenisons of user coma model

parser.add_argument("--kernelMetaFile",type=str)
parser.add_argument("--utcStartTime",type=str)

parser.add_argument("--userR",type=float)                               # distance in km from nucleus center
parser.add_argument("--userPhaseAngle",type=float)
parser.add_argument("--userLatitude",type=float)
parser.add_argument("--userAlpha",type=float)
parser.add_argument("--userBeta",type=float)
parser.add_argument("--userGamma",type=float)

args = parser.parse_args()

modelCase = args.modelCase
pointingCase = args.pointingCase
instrumentSelector = args.instrumentSelector
outputDir = args.outputDir

dataFile = args.dataFile
isDust = args.isDust
minSize = args.minSize
maxSize =args.maxSize

QHaser = args.QHaser
vHaser = args.vHaser
tpHaser = args.tpHaser
tdHaser = args.tdHaser
kernelMetaFile = args.kernelMetaFile
utcStartTime = args.utcStartTime

userR = args.userR
userPhaseAngle = args.userPhaseAngle
userLatitude = args.userLatitude
userAlpha = args.userAlpha
userBeta = args.userBeta
userGamma = args.userGamma

userDataFile = args.userDataFile
userDelimiter = args.userDelimiter
userDim = args.userDim
userNrOfGeaderRows= args.userNrOfHeaderRows



print '##########################################'
print 'modelCase     :', modelCase
print 'instrument    :', instrumentSelector
print 'kernelMetaFile:',kernelMetaFile
print 'utcStartTime  :', utcStartTime
print 'QHaser        :', QHaser
print 'vHaser        :',vHaser
print 'tdHaser       :', tdHaser
print 'tpHaser       :',tpHaser
print '##########################################'


if modelCase == 0:
    ############################################
    # check if 1d or 2d case
    ############################################
    
    dataDir = os.path.split(dataFile)[0]
    filenames = os.listdir(dataDir)
    allFilenamesInOneString = ''
    for filename in filenames:
        allFilenamesInOneString += filename
        
    if '2d' in allFilenamesInOneString:
        dim = 2
    elif '1d' in allFilenamesInOneString:
        dim = 1
    else:
        dim = 0
        print 'Could not detect number of dimensions of dsmc case. Exiting now.'
        sys.exit()

elif modelCase == 1:
    dim = 1

elif modelCase == 2:
    dim = userDim
    
print 'dimensions:', dim

if pointingCase == 0:
    #################################################
    # get rosetta coordinates from spice
    #################################################
    spice.furnsh(kernelMetaFile)                    # load spice kernels
    et = spice.str2et(utcStartTime)                 # ephemeris time calculated from utc
    rRosetta,lightTime = spice.spkpos("ROSETTA", et, "67P/C-G_CSO", "NONE", "CHURYUMOV-GERASIMENKO")        # s/c coordinates in CSO frame of reference
    rRosetta = np.array(rRosetta) * 1000            # transform km to m
    R = spice.pxform("ROS_SPACECRAFT", "67P/C-G_CSO", et);      # create rotation matrix R to go from instrument reference frame to CSO
    
    print 'distance from comet:', np.sqrt( np.sum( rRosetta**2 ) )
    
elif pointingCase == 1:
    x0 = np.array([-userR, 0, 0])
    rRosetta = rotations.rotateVector(x0,userPhaseAngle, userLatitude)
    ei, ej, ek = rotations.rotateCoordinateSystem2(userPhaseAngle, userLatitude, userAlpha, userBeta, userGamma)
    
    R = rotations.createRotationMatrix(ei,ej,ek)
    
    
print 'rRosetta:', rRosetta

########################################################
# load data
########################################################
if modelCase == 0:
    if isDust:
        print 'dust case'
        numberDensityIndices,allSizeIntervals = getAllDustIntervalIndices(dataFile,dim)
        
        x,y,n = loadDustData(allSizeIntervals,numberDensityIndices,dim,dataFile)
        
    else:
        x,y,n = loadGasData(dataFile, dim)
            
elif modelCase == 1:
    x,n = haserModel(QHaser,vHaser,tpHaser,tdHaser)
    y = None
    
elif modelCase == 2:
    x,y,n = loadGasData(userDataFile,dim,True,userDelimiter,userNrOfGeaderRows)
        
        
##############################################################
# triangulation and interpolation for 2d case
if dim == 1:
    pass
elif dim == 2:
    triangles = mtri.Triangulation(x,y)
    Interpolator = mtri.LinearTriInterpolator(triangles, n)
    
print 'interpolation done'
#############################################################


#############################################################
# instrument specific definitions
#############################################################

if instrumentSelector == 1:                 # osiris wac
    pixelsX = 256                           # nr of pixels along x axis
    pixelsY = 256                           # nr of pixels along y axis
    phi_x = 12/2                            # instrument FOV in x (half opening angle) in degrees
    phi_y = 12/2                            # instrument FOV in y (half opening angle) in degrees
    iFOV = 0.000993                         # pixel FOV in rad
    pixelSize = 1                           # area of one pixel
    
elif instrumentSelector == 2:               # osiris nac
    pixelsX = 512
    pixelsY = 512 
    phi_x = 3/2
    phi_y = 3/2
    iFOV = 0.0000188
    pixelSize = 1

elif instrumentSelector == 3:               # alice
    N_oversampleX = 24
    N_oversampleY = 20
    
    pixelsX = 19 * N_oversampleX
    pixelsY = 1 * N_oversampleY
    phi_x = 5.852/2
    phi_y = 0.1/2
    
    iFOV = (phi_x*2 / 180 * np.pi / (19*N_oversampleX)) * ( phi_y*2 / 180 * np.pi / (N_oversampleY))
    
    print iFOV
    pixelSize = 1
    
    v_sun = alice.get_v_sun(kernelMetaFile, utcStartTime)
    gFactor = alice.getGFactorFromDB()
    
    
elif instrumentSelector == 4:               # miro
    pixelsX = 1
    pixelsY = 1
    phi_x = 0.33336/2
    phi_y = 0.36666/2
    iFOV = 0.36666
    pixelSize = 1
    
elif instrumentSelector == 5:               # virtis m
    pixelsX = 256
    pixelsY = 256
    phi_x = 3.6669/2
    phi_y = 3.6669/2
    iFOV = 0.00025
    pixelSize = 1

elif instrumentSelector == 6:               # virtis h
    pixelsX = 1
    pixelsY = 3
    phi_x = 0.0334
    phi_y = 0.1
    iFOV = 0.000583
    pixelSize = 1
    
lx = 2 * np.sin(phi_x / 180 * np.pi)
ly = 2 * np.sin(phi_y / 180 * np.pi)

if pixelsX > 1:
    dx = lx / (pixelsX-1)
else:
    dx = lx

if pixelsY > 1:
    dy = ly / (pixelsY-1)
else:
    dy = ly
    
print 'entering pixel loop'

p = np.zeros(3)
ccd = np.zeros((pixelsX,pixelsY))
for i in range(pixelsX):
    for j in range(pixelsY):
        if pointingCase == 0:
            p[0] = i*dx - lx/2 + dx/2
            p[1] = j*dy - ly/2 + dy/2
            p[2] = 1
            p /= np.sqrt(np.sum(p**2))          # p = pointing vector in instrument coordinates
            pSpice = spice.mxv(R, p)            # pointing vector in CSO frame of reference (sun on +x axis)
            p[0] = -pSpice[0]                   # p = pointing vector in simulation coordinates (sun on -x axis)
            p[1] = -pSpice[1]
            p[2] = pSpice[2]
            
            rRay = 5 * np.array([value for value in rRosetta])  # ray for line of sight, starting at rosetta position (transform to coords with sun on -x)
            rRay[0] *= -1
            rRay[1] *= -1
        else:
            p[1] = i*dx - lx/2 + dx/2
            p[2] = j*dy - ly/2 + dy/2
            p[0] = 1
            p /= np.sqrt(np.sum(p**2))          # p = pointing vector in instrument coordinates
            #p = spice.vrotv(p,(0,-1,0),np.pi/2)
            p = R.dot(p)
            #print 'p:',p
            rRay = np.array([value for value in rRosetta])
            #print 'rRosetta:', rRosetta
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
            
            if dim == 1:
                xTravel.append(distance)
            elif dim == 2:
                xTravel.append((rRay[0], np.sqrt(np.sum( rRay[1]**2 + rRay[2]**2) ) ))
            
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
        
        if dim == 1:
            nRay = np.interp(xTravel,x,n)
        elif dim == 2:
            nRay = Interpolator.__call__(xTravel[:,0], xTravel[:,1])        # interpolated local number density
        
        if instrumentSelector == 3:
            brightness = alice.calculateColumn(nRay,dTravel,0,iFOV)
            ccd[i][j] = brightness
        else:
            columnDensity = np.trapz(nRay,dTravel)                          # integration along xRay
            ccd[i][j] = columnDensity
            
    print i
print 'pixel loop done'

ccd = np.array(ccd)

if instrumentSelector == 3:
    ccdFinal = alice.calculateBrightness(pixelsX,N_oversampleX, N_oversampleY,ccd,gFactor)
else:
    ccdFinal = ccd

######################################################
# write results to file
######################################################
f = open(outputDir + '/result.txt','w')
for row in ccd:
    for value in row:
        f.write('%e,' % value)
    f.write('\n')
f.close()


######################################################
# plot results
#######################################################
plotResult(ccdFinal, outputDir, 'result.png', instrumentSelector,runDetails = args,showPlot=True)