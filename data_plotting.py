#!/usr/bin/env python
from __future__ import division
import os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
import datetime

from data_loaders import load_in_situ_output

rcParams.update({'figure.autolayout': True})


def plot_in_situ(args):

    path = args.StringOutputDir
    pltTitle = 'ICES in-situ tool\n'
    plt.xkcd()
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    if args.iModelCase == 0:
        pltTitle += 'model: DSMC, case: %s\n' % (os.path.split(args.StringDataFileDSMC)[0].split('/')[-1])
        pltTitle += 'spice: %s' % (os.path.split(args.StringKernelMetaFile)[1].split('.')[0])
    elif args.iModelCase == 1:
        pltTitle += 'model: HASER, spice: %s\n' % (os.path.split(args.StringKernelMetaFile)[1].split('.')[0])
        pltTitle += 'Q = %.2e [#/s], v = %.0f [m/s], tp = %.2e [s]' % (args.QHaser, args.vHaser, args.tpHaser)

    filenames = [path + '/' + filename for filename in os.listdir(path) if filename.split('.')[-1] == 'out']
    for filename in filenames:
        species = os.path.split(filename)[1].split('.')[0]
        dates_SC, r_SC, n_SC = load_in_situ_output(filename)
        ax1.semilogy(dates_SC, n_SC, label=species, lw=2)

    ax2.plot(dates_SC, r_SC / 1000, '-k', lw=2)
    ax1.set_title(pltTitle)
    ax1.set_ylabel('Number density [#/m3]')
    ax1.grid(True)
    ax1.legend()

    ax2.set_ylabel('Distance from comet center [km]')
    ax2.grid(True)
    plt.gcf().autofmt_xdate()
    plt.savefig(path + '/' + 'result.png')
    plt.show()


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
