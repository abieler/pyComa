#!/usr/bin/env python
from __future__ import division
import numpy as np


def getAllDustIntervalIndices(filename, dim):
    '''
    returna a list containing x,y indices and the indices of
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
                        allSizeIntervals.append(
                                    np.float(element.split('=')[1][0:-2]))
                    allIndices.append(j)
            break
    f.close()
    print len(allSizeIntervals), "size intervals found."
    return allIndices, allSizeIntervals


def loadGasData(dataFile, dim, userData=False, userDelimiter=',',
                userNrOfHeaderRows=0):

    if not userData:
        '''
        load gas data files from amps in 1d or 2d, and return values
        for x,y and number density n.
        '''

        f = open(dataFile, 'r')
        for line in f:
            if 'VARIABLES' in line:
                variables = line.split(',')

                for element, j in zip(variables, range(len(variables))):
                    if '"n"' in element:
                        densityIndex = j
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


def loadDustData(allSizeIntervals, numberDensityIndices, dim, dataFile):

    userIndices = []
    for size, index in zip(allSizeIntervals, numberDensityIndices):
        if minSize <= size <= maxSize:
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