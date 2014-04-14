#!/opt/local/bin/python2.7
#!/usr/bin/env python2.7

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

import spice
from data_loaders import *
import spice_functions


parser = argparse.ArgumentParser()
parser.add_argument("--iModelCase", type=int, choices=[0, 1, 2], help='0: dsmc model, 1: haser model, 2: user model')
parser.add_argument("--iPointingCase", type=int, choices=[0, 1], help='0: spice pointing, 1: user pointing')
parser.add_argument("--iInstrumentSelector", type=int, choices=[1, 2, 3, 4, 5, 6])
parser.add_argument("--StringOutputDir", type=str)

parser.add_argument("--StringDSMCdir", type=str)
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
print 'dsmcDur:', args.StringDSMCdir
print '--' * 20

x_SC, y_SC, z_SC = spice_functions.get_coordinates(args.StringUtcStartTime, 'ROSETTA', 'J2000', "None",
                                                   "CHURYUMOV-GERASIMENKO", args.UtcStopTime, args.nDeltaT)

n_SC_all = []
species_all = []

firstSpecies = True
filenames = os.listdir(dsmcDir)
allFilenamesInOneString = ''
for filename in filenames:
    allFilenamesInOneString += filename


if '2d' in allFilenamesInOneString:
    dim = 2
elif '1d' in allFilenamesInOneString:
    dim = 1
else:
    dim = 0

print "dimension: ", dim
for filename in filenames:
    if ("Dust" not in filename) and (filename.split('.')[-1] == 'dat'):
        #####################################################
        # load data
        #####################################################
        species = filename.split('.')[-2]
        print "load dsmc data", species

        if dim == 1:
            r, n = np.genfromtxt(dsmcDir + '/' + filename, skip_header=2, autostrip=True, usecols=(0, 3), unpack=True)
            n_SC = np.interp(r_SC, r, n)
        elif dim == 2:
            if triImport:
                file = open(dsmcDir + '/' + filename, 'r')
                for line in file:
                    if 'VARIABLES' in line:
                        if 'Gx' in line:
                            densityIndex = 7
                        else:
                            densityIndex = 5

                x, y, n = np.genfromtxt(dsmcDir + '/' + filename, skip_header=3, skip_footer=155236,
                                        autostrip=True, usecols=(0, 1, densityIndex), unpack=True)
                triangles = mtri.Triangulation(x, y)
                Interpolator = mtri.LinearTriInterpolator(triangles, n)
                n_SC = Interpolator.__call__(x_SC, np.sqrt(y_SC**2 + z_SC**2))
                print "successfully imported 2d data"
            else:
                print "2d cases not supported. install matplotlib version 0.1.3.1 or newer"

        n_SC_all.append(n_SC)
        species_all.append(species)
        #####################################################
        # create dates for the figure
        #####################################################
        if firstSpecies:
            yy = np.int(utcStartTime.split('-')[0])
            mm = np.int(utcStartTime.split('-')[1])
            dd = np.int(utcStartTime.split('-')[2][0:2])
            HH = np.int(utcStartTime.split('T')[1].split(':')[0])
            MM = np.int(utcStartTime.split('T')[1].split(':')[1])
            SS = np.int(utcStartTime.split('T')[1].split(':')[2])

            dates = []

            t = datetime.datetime(yy,mm,dd,HH,MM,SS)
            dt = datetime.timedelta(seconds = deltaT)

            while len(dates) < len(n_SC):
                dates.append(t)
                t += dt
            firstSpecies = False
        
        ############################################
        # write results to file
        ###############################################

        file = open(outputDir + '/' + species + '.out','w')
        file.write('Local number densities for the rosetta spacecraft at selected dates. Comet is at (0,0,0) with the sun on the positive x axis.(inf,0,0)\n')
        file.write('DSMC case: %s\n' %(dsmcDir.split('/')[-1]))
        file.write('spice kernel: %s\n' %(mkFile.split('/')[-1]))
        file.write('date,x[m],y[m],z[m],numberDensity [1/m3]\n')
        #print "outfile:", outputDir + '/' + species + '.out'
        for dd,xx,yy,zz,nn in zip(dates,x_SC,y_SC,z_SC,n_SC):
            file.write("%s,%e,%e,%e,%e\n" %(dd,xx,yy,zz,nn))
        file.close()

        ####################################################
        # plot results
        ####################################################

        plt.semilogy(dates,n_SC,label=species)


n_SC_all = np.array(n_SC_all)
file = open(outputDir + '/' + 'allSpecies.out','w')
file.write('Local number densities [1/m3] for the rosetta spacecraft at selected dates. Comet is at (0,0,0) with the sun on the positive x axis.(inf,0,0)\n')
file.write('DSMC case: %s\n' %(dsmcDir.split('/')[-1]))
file.write('spice kernel: %s\n' %(mkFile.split('/')[-1]))
file.write('date,x[m],y[m],z[m],')
for spec in species_all:
    file.write("nrDensity_"+ spec + ",")
file.write("\n")
for dd,xx,yy,zz,i in zip(dates,x_SC,y_SC,z_SC,range(len(x_SC))):
    file.write("%s,%e,%e,%e," %(dd,xx,yy,zz))
    for value in n_SC_all[:,i]:
        file.write("%e," %value)
    file.write("\n")
file.close()

plt.title('Neutral densitiy along trajectory')
plt.ylabel('Number density [1/m^3]')
plt.grid(True)
plt.gcf().autofmt_xdate()
plt.legend()
#plt.show()
plt.savefig(outputDir + '/' + 'result.png')
