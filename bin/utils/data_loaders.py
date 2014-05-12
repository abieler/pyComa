#!/usr/bin/env python
from __future__ import division
import numpy as np
import os
import datetime
import time
from pandas import read_csv


def load_hybrid2_data(filename):
    x, y, z, n = np.genfromtxt(filename, usecols=(0, 1, 2, 13), unpack=True)
    return np.array([x, y, z]), n


def get_iDim(args):

    if args.iModelCase == 0:
        ############################################
        # check if 1d or 2d case
        ############################################

        filenames = os.listdir(os.path.split(args.StringDataFileDSMC)[0])
        allFilenamesInOneString = ''
        for filename in filenames:
            allFilenamesInOneString += filename

        if '2d' in allFilenamesInOneString:
            iDim = 2
        elif '1d' in allFilenamesInOneString:
            iDim = 1
        else:
            iDim = 0
            if args.iModelCase == 0:
                print 'Could not detect iDim of dsmc case. Exiting now.'
            sys.exit()

    elif args.iModelCase == 1:
        iDim = 1

    elif args.iModelCase == 2:
        iDim = args.iDimUser

    print 'iDim:', iDim
    return iDim


def getAllDustIntervalIndices(filename, dim):
    '''
    returna a list containing the indices of
    number density variables in the amps dust output file.
    also returns a list of all size intervals contained in the data file.
    '''

    print 'getting dust size intervals from', filename

    allSizeIntervals = []
    allIndices = []

    f = open(filename, 'r')
    for line in f:
        if 'VARIABLES' in line:
            variables = line.split(',')
            for element, j in zip(variables, range(len(variables))):
                condition1 = ('NUMBERDENSITY' in element)
                condition2 = ('R=' in element)
                condition3 = '"X"' in element
                condition4 = '"Y"' in element
                if condition1 and condition2 or condition3 or condition4:
                    if '"X"' in element or '"Y"' in element:
                        pass
                    else:
                        allSizeIntervals.append(np.float(element.split('=')[1][0:-2]))
                    allIndices.append(j)
            break
    f.close()
    #print len(allSizeIntervals), "dust size intervals found."
    #for size, index in zip(allSizeIntervals, allIndices):
    #    print 'size: %.2e m, index: %i' % (size, index)
    return allIndices, allSizeIntervals


def load_gas_data(dataFile, dim, userData=False, userDelimiter=',',
                  userNrOfHeaderRows=0):

    if not userData:
        '''
        load gas data files from amps in 1d or 2d, and return values
        for x,y and number density n.
        '''
        print 'opening data file', dataFile
        f = open(dataFile, 'r')
        for line in f:
            if 'VARIABLES' in line:
                variables = line.split(',')

                for element, j in zip(variables, range(len(variables))):
                    if '"n"' in element:
                        densityIndex = j
                        print 'found density index!!', j
                        break
        f.close()

        if dim == 1:
            dataIndices = [0, densityIndex]
            headerRows = 2
            footerRows = 0
            varNames = ['x', 'n']
        elif dim == 2:
            dataIndices = [0, 1, densityIndex]
            headerRows = 3
            footerRows = 155236
            varNames = ['x', 'y', 'n']

        # load data from amps
        print 'usecols: ', dataIndices
        data = read_csv(dataFile, dtype=float, skiprows=headerRows,
                        skip_footer=footerRows, usecols=(dataIndices),
                        names=varNames, sep=r"\s+")
        print 'gas data loaded'

        if dim == 1:
            x = data['x']
            y = None
            n = data['n']

        elif dim == 2:
            x = data['x']
            y = data['y']
            n = data['n']                           # number density
    else:

        print 'userDim:', dim

        if dim == 1:
            dataIndices = [0, 1]
        elif dim == 2:
            dataIndices = [0, 1, 2]

        data = np.genfromtxt(dataFile, dtype=float,
                             skip_header=userNrOfHeaderRows,
                             delimiter=userDelimiter,
                             usecols=(0, 1))

        if dim == 1:
            x = data[:, 0]
            y = None
            n = data[:, 1]
        elif dim == 2:
            x = data[:, 0]
            y = data[:, 1]
            n = data[:, 2]

    return x, y, n


def loadGasData(dataFile, dim, userData=False, userDelimiter=',',
                userNrOfHeaderRows=0):

    if not userData:
        '''
        load gas data files from amps in 1d or 2d, and return values
        for x,y and number density n.
        '''
        print 'opening data file', dataFile
        f = open(dataFile, 'r')
        for line in f:
            if 'VARIABLES' in line:
                variables = line.split(',')

                for element, j in zip(variables, range(len(variables))):
                    if '"n"' in element:
                        densityIndex = j
                        print 'found density index!!', j
                        break
        f.close()

        if dim == 1:
            dataIndices = [0, densityIndex]
            headerRows = 2
            footerRows = 0
        elif dim == 2:
            dataIndices = [0, 1, densityIndex]
            headerRows = 3
            footerRows = 155236

        # load data from amps
        print 'usecols: ', dataIndices
        data = np.genfromtxt(dataFile, dtype=float, skip_header=headerRows,
                             skip_footer=footerRows, usecols=(dataIndices))
        print 'gas data loaded'

        if dim == 1:
            x = data[:, 0]
            y = None
            n = data[:, 1]

        elif dim == 2:
            x = data[:, 0]
            y = data[:, 1]
            n = data[:, 2]                           # number density
    else:

        print 'userDim:', dim

        if dim == 1:
            dataIndices = [0, 1]
        elif dim == 2:
            dataIndices = [0, 1, 2]

        data = np.genfromtxt(dataFile, dtype=float,
                             skip_header=userNrOfHeaderRows,
                             delimiter=userDelimiter,
                             usecols=dataIndices)

        if dim == 1:
            x = data[:, 0]
            y = None
            n = data[:, 1]
        elif dim == 2:
            x = data[:, 0]
            y = data[:, 1]
            n = data[:, 2]

    return x, y, n


def load_user_trajectory(args):

    print 'loading trajectory from: ', args.StringUserTrajectoryFile
    if args.DelimiterTraj.lower() == 'space':
        args.DelimiterTraj = " "

    print 'trajectory file delimiter:"%s"' % args.DelimiterTraj
    x, y, z = np.genfromtxt(args.StringUserTrajectoryFile, dtype=float,
                            skip_header=args.nHeaderRowsData,
                            delimiter=args.DelimiterTraj, unpack=True)
    return x, y, z


def load_user_data(DataFile, iDim, Delimiter, nHeaderRows):

    if iDim == 1:
        DataIndices = [0, 1]
    elif iDim == 2:
        DataIndices = [0, 1, 2]
    elif iDim == 3:
        DataIndices = [0, 1, 2, 3]

    if Delimiter == 'SPACE':
        Delimiter = ' '

    data = np.genfromtxt(DataFile, dtype=float,
                         skip_header=nHeaderRows,
                         delimiter=Delimiter,
                         usecols=DataIndices)

    if iDim == 1:
        x = data[:, 0]
        y = None
        n = data[:, 1]
    elif iDim == 2:
        x = data[:, 0]
        y = data[:, 1]
        n = data[:, 2]

    return x, y, n


def load_dust_data_full(allSizeIntervals, numberDensityIndices, iDim, dataFile, args):

    if iDim == 1:
        allIndices = [0]
    elif iDim == 2:
        allIndices = [0, 1]

    UserIndices = [i for i, size in zip(numberDensityIndices[iDim:], allSizeIntervals)
                   if (args.DustSizeMin <= size <= args.DustSizeMax)]

    for i in UserIndices:
        allIndices.append(i)        # index of number density
        allIndices.append(i+1)      # index of mass density

    data = np.genfromtxt(dataFile, dtype=float, skip_header=3,
                         skip_footer=155236, usecols=allIndices)

    #print 'user indices', UserIndices
    #print 'all indices', allIndices

    if iDim == 1:
        x = data[:, 0]
        y = None
    elif iDim == 2:
        x = data[:, 0]
        y = data[:, 1]

    numberDensities = data[:, iDim::2]
    massDensities = data[:, iDim+1::2]

    return x, y, numberDensities, massDensities


def load_dust_data(allSizeIntervals, numberDensityIndices, dim, dataFile, args):

    print len(allSizeIntervals)
    print len(numberDensityIndices)

    userIndices = []
    for size, index in zip(allSizeIntervals, numberDensityIndices):
        if args.DustSizeMin <= size <= args.DustSizeMax:
            if dim == 1:
                userIndices.append(index + 1)
                # add +1 becaus x is also contained in data
            elif dim == 2:
                userIndices.append(index+2)
                # add +2 because x and y are also contained in data

    data = np.genfromtxt(dataFile, dtype=float, skip_header=3,
                         skip_footer=155236, usecols=numberDensityIndices)
    print 'dust data loaded'

    x = data[:, 0]
    if dim == 1:
        y = None
    elif dim == 2:
        y = data[:, 1]
    elif dim == 3:
        y = None         # 3d not implemented yet!!!

    n = np.zeros(len(data[:, 0]))
    for i in userIndices:
        n += data[:, i]

    return x, y, n


def load_in_situ_output(filename):
    i = 0
    dates_SC, r_SC, n_SC = [], [], []
    file = open(filename, 'r')
    for line in file:
        if i >= 0:
            try:
                dd, xx, yy, zz, rr, nn = line.split(',')
                r_SC.append(np.float(rr))
                n_SC.append(np.float(nn))
                try:
                    dates_SC.append(datetime.datetime.strptime(dd, '%Y-%m-%d %H:%M:%S'))
                except:
                    dates_SC.append(dd)
            except Exception, e:
                pass
        i += 1
    file.close()
    return dates_SC, np.array(r_SC), np.array(n_SC)
