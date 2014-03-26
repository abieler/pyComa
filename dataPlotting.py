#!/usr/bin/env python
from __future__ import division
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})


def plotResult(ccd,outputDir,outFileName,instrumentSelector,runDetails,showPlot=False):
    
    
    if instrumentSelector == 1:
        pltInstrument = 'OSIRIS WAC'
    elif instrumentSelector == 2:
        pltInstrument = 'OSIRIS NAC'
    elif instrumentSelector == 3:
        pltInstrument = 'ALICE'
    elif instrumentSelector == 4:
        pltInstrument = 'MIRO'
    elif instrumentSelector == 5:
        pltInstrument = 'VIRTIS M'
    elif instrumentSelector == 6:
        pltInstrument = 'VIRTIS H'
    
    
    if runDetails.modelCase == 0:
        pltModelCase = 'dsmc model'
    elif runDetails.modelCase == 1:
        pltModelCase = 'haser model'
        pltQHaser = str(runDetails.QHaser)
        pltvHaser = str(runDetails.vHaser)
        plttdHaser = str(runDetails.tdHaser)
        plttpHaser = str(runDetails.tpHaser)
    elif runDetails.modelCase == 2:
        pltModelCase = 'user model'
        
    if runDetails.pointingCase == 0:
        pltPointing = 'spice pointing'
        pltSpiceCase = runDetails.kernelMetaFile
        pltUtcStartTime = runDetails.utcStartTime
    elif runDetails.pointingCase == 1:
        pltPointing = 'manual pointing'
        pltR = ''
        pltPhaseAngle = ''
        pltLatitude = ''
        pltEl = ''
        pltAz = ''
        
    
    pltTitle = '%s, %s, %s\n' %(pltInstrument,pltModelCase, pltPointing)
    if runDetails.pointingCase == 0:
        pltTitle += '%s, %s\n' %( os.path.split(runDetails.kernelMetaFile)[1], pltUtcStartTime )
    else:
        pltTitle += 'R:%i, PA:%i, LAT:%i, a:%i, b:%i, c:%i\n' %(runDetails.userR,runDetails.userPhaseAngle,runDetails.userLatitude,runDetails.userAlpha,runDetails.userBeta, runDetails.userGamma)
    
    if runDetails.modelCase == 0:
        pltTitle += '%s\n' % (os.path.split(runDetails.dataFile)[1])
    elif runDetails.modelCase == 1:
        pltTitle += 'Q: %.2e [#/s],  v: %i [m/s],  Tp: %.0e [s]' %(runDetails.QHaser, runDetails.vHaser, runDetails.tpHaser)
        
        if runDetails.tdHaser == None:
            pltTitle += '\n'
        else:
            pltTitle += ', Td: %0.e [s]\n' % ( runDetails.tdHaser )
    elif runDetails.modelCase == 2:
        pass
    
    
    plt.figure()
    if instrumentSelector == 3:           # alice
        plt.semilogy(range(5,24),ccd/10**10,'-ok',linewidth=2)
        plt.grid(True)
        plt.xlabel('Pixel Number')
        plt.ylabel('Rayleigh')
        plt.xticks(range(5,24))
        plt.xlim((5,23))
    else:
        plt.contourf(ccd,200)
        plt.colorbar()
        
    plt.title(pltTitle)
    if showPlot:
        plt.show()
