#!/usr/bin/env python
from __future__ import division
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})


def plot_result(ccd, StringOutputDir, StringOutFileName, iInstrumentSelector, RunDetails, DoShowPlot=False):

    if iInstrumentSelector == 1:
        pltInstrument = 'OSIRIS WAC'
    elif iInstrumentSelector == 2:
        pltInstrument = 'OSIRIS NAC'
    elif iInstrumentSelector == 3:
        pltInstrument = 'ALICE'
    elif iInstrumentSelector == 4:
        pltInstrument = 'MIRO'
    elif iInstrumentSelector == 5:
        pltInstrument = 'VIRTIS M'
    elif iInstrumentSelector == 6:
        pltInstrument = 'VIRTIS H'

    if RunDetails.iModelCase == 0:
        pltModelCase = 'dsmc model'
    elif RunDetails.iModelCase == 1:
        pltModelCase = 'haser model'
        pltQHaser = str(RunDetails.QHaser)
        pltvHaser = str(RunDetails.vHaser)
        plttdHaser = str(RunDetails.tdHaser)
        plttpHaser = str(RunDetails.tpHaser)
    elif RunDetails.iModelCase == 2:
        pltModelCase = 'user model'

    if RunDetails.iPointingCase == 0:
        pltPointing = 'spice pointing'
        pltSpiceCase = RunDetails.StringKernelMetaFile
        pltUtcStartTime = RunDetails.StringUtcStartTime
    elif RunDetails.iPointingCase == 1:
        pltPointing = 'manual pointing'
        pltR = ''
        pltPhaseAngle = ''
        pltLatitude = ''
        pltEl = ''
        pltAz = ''

    pltTitle = '%s, %s, %s\n' % (pltInstrument, pltModelCase, pltPointing)
    if RunDetails.iPointingCase == 0:
        pltTitle += '%s, %s\n' % (os.path.split(RunDetails.StringKernelMetaFile)[1],
                                  pltUtcStartTime)
    else:
        pltTitle += 'R: %i km, PA: %i, LAT: %i, a: %i, b: %i, c: %i\n' %\
                    (RunDetails.UserR / 1000, RunDetails.UserPhaseAngle,
                     RunDetails.UserLatitude, RunDetails.UserAlpha,
                     RunDetails.UserBeta, RunDetails.UserGamma)

    if RunDetails.iModelCase == 0:
        pltTitle += '%s\n' % (os.path.split(RunDetails.StringDataFileDSMC)[1])
    elif RunDetails.iModelCase == 1:
        pltTitle += 'Q: %.2e [#/s],  v: %i [m/s],  Tp: %.0e [s]' % \
                    (RunDetails.QHaser, RunDetails.vHaser, RunDetails.tpHaser)

        if RunDetails.tdHaser is None:
            pltTitle += '\n'
        else:
            pltTitle += ', Td: %0.e [s]\n' % (RunDetails.tdHaser)
    elif RunDetails.iModelCase == 2:
        pass

    plt.figure()
    if iInstrumentSelector == 3:           # alice
        plt.plot(range(5, 24), ccd, '-ok', linewidth=2)
        plt.grid(True)
        plt.xlabel('Pixel Number')
        plt.ylabel('Flux [photons / m2 / s]')
        plt.xticks(range(5, 24))
        plt.xlim((5, 23))
    else:
        plt.contourf(ccd, 200)
        plt.colorbar()

    plt.title(pltTitle)
    if DoShowPlot:
        plt.show()
