#!/usr/bin/env python
from __future__ import division
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams

import bokeh.plotting as bplt
from bokeh.objects import Range1d

from data_loaders import load_in_situ_output

rcParams.update({'figure.autolayout': True})


def plot_result_insitu(args):

    print 'plotting results...'
    path = args.StringOutputDir
    nLinesPerFig = 4
    pltTitle = build_plot_title(args, 'insitu')

    filenames = [path + '/' + filename for filename in os.listdir(path) if filename.split('.')[-1] == 'out']
    all_n_SC = []
    all_species = []
    for filename, i in zip(filenames, range(len(filenames))):
        species = os.path.split(filename)[1].split('.')[0]
        dates_SC, r_SC, n_SC = load_in_situ_output(filename)
        all_n_SC.append(n_SC)
        all_species.append(species)

    if args.StringPlotting.lower() == 'bokeh':
        create_plot_insitu_bokeh(args, all_n_SC, all_species, dates_SC, r_SC, nLinesPerFig, pltTitle)
    elif args.StringPlotting.lower() == 'matplotlib':
        create_plot_insitu_matplotlib(args, all_n_SC, all_species, dates_SC, r_SC, nLinesPerFig, pltTitle)


def create_plot_insitu_bokeh(args, all_n_SC, all_species, dates_SC, r_SC, nLinesPerFig, pltTitle):

    pltColors = ['black', 'red', '#006400', 'blue', '#8B008B', '#66CDAA', '#FF8C00', 'cyan']
    for species, i in zip(all_species, range(len(all_species))):
        if (i % nLinesPerFig == 0):
            bplt.output_file(args.StringOutputDir + '/' + 'result_%i.html' % (i // nLinesPerFig))
            bplt.figure(x_axis_type='datetime')
            bplt.hold()
        bplt.line(dates_SC, np.log10(all_n_SC[i]), legend=all_species[i], line_width=2, line_color=pltColors[i % nLinesPerFig])
        bplt.xaxis().axis_label = 'Time'
        bplt.yaxis().axis_label = 'log10(n) [#/m3]'
        if (i % nLinesPerFig == (nLinesPerFig - 1)) or (i == len(all_species) - 1):
            bplt.grid().grid_line_alpha = 0.4
            bplt.figure(x_axis_type='datetime')
            bplt.line(dates_SC, r_SC/1000, line_width=2, line_color='black')
            bplt.grid().grid_line_alpha = 0.4
            bplt.xaxis().axis_label = 'Time'
            bplt.yaxis().axis_label = 'Distance from comet center [km]'
            bplt.save()

    if args.DoShowPlots:
        bplt.show()
    else:
        print 'not showing results on screen'


def create_plot_insitu_matplotlib(args, all_n_SC, all_species, dates_SC, r_SC, nLinesPerFig, pltTitle):

    for species, i in zip(all_species, range(len(all_species))):
            if (i % nLinesPerFig == 0):   # make new figure for every 4 species
                fig = plt.figure(figsize=(14, 14))
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
            ax1.semilogy(dates_SC, all_n_SC[i], label=all_species[i], lw=2)

            # plot distance from comet and save figure
            if (i % nLinesPerFig == (nLinesPerFig - 1)) or (i == len(all_species) - 1):
                ax1.set_title(pltTitle)
                ax1.set_ylabel('Number density [#/m3]')
                ax1.grid(True)
                ax1.legend(framealpha=1)

                ax2.plot(dates_SC, r_SC / 1000, '-k', lw=2)
                ax2.set_ylabel('Distance from comet center [km]')
                ax2.grid(True)
                plt.gcf().autofmt_xdate()
                plt.savefig(args.StringOutputDir + '/' + 'result_%i.png' % ((i-(nLinesPerFig-1)) // nLinesPerFig))

    if args.DoShowPlots:
        plt.show()
    else:
        print 'not showing results on screen'


def plot_result_LOS(ccd, StringOutputDir, StringOutFileName, iInstrumentSelector, args, DoShowPlot=False):

    pltTitle = build_plot_title(args, 'LOS')

    if args.StringPlotting == 'matplotlib':
        if args.iInstrumentSelector in [1, 2, 4]:
            create_plot_LOS_2d(args, ccd, pltTitle)
        else:
            create_plot_LOS_1d(args, ccd, pltTitle)

    elif args.StringPlotting == 'bokeh':
        print 'bokeh plotting not supported yet for LOS cases'
        create_plot_LOS_2d_bokeh(args, ccd, pltTitle)


def create_plot_LOS_2d_bokeh(args, ccd, pltTitle):

    nPixels = len(ccd[:, 0])
    bplt.output_file(args.StringOutputDir + '/' + 'result.html')
    bplt.image(image=[ccd], x=[0], y=[0], dw=[nPixels], dh=[nPixels], palette=["Spectral-11"],
               x_range=Range1d(start=0, end=nPixels), y_range=Range1d(start=0, end=nPixels))
    print pltTitle
    bplt.curplot().title = pltTitle
    bplt.save()

    if args.DoShowPlots:
        bplt.show()
    else:
        print 'not showing results on screen'



def build_plot_title(args, measurement='LOS'):

    if measurement == 'LOS':
        if args.iInstrumentSelector == 1:
            pltInstrument = 'OSIRIS WAC'
        elif args.iInstrumentSelector == 2:
            pltInstrument = 'OSIRIS NAC'
        elif args.iInstrumentSelector == 3:
            pltInstrument = 'ALICE'
        elif args.iInstrumentSelector == 4:
            pltInstrument = 'MIRO'
        elif args.iInstrumentSelector == 5:
            pltInstrument = 'VIRTIS M'
        elif args.iInstrumentSelector == 6:
            pltInstrument = 'VIRTIS H'

        if args.iModelCase == 0:
            pltModelCase = 'dsmc model'
        elif args.iModelCase == 1:
            pltModelCase = 'haser model'
            pltQHaser = str(args.QHaser)
            pltvHaser = str(args.vHaser)
            plttdHaser = str(args.tdHaser)
            plttpHaser = str(args.tpHaser)
        elif args.iModelCase == 2:
            pltModelCase = 'user model'

        if args.iPointingCase == 0:
            pltPointing = 'spice pointing'
            pltSpiceCase = args.StringKernelMetaFile
            pltUtcStartTime = args.StringUtcStartTime
        elif args.iPointingCase == 1:
            pltPointing = 'manual pointing'
            pltR = ''
            pltPhaseAngle = ''
            pltLatitude = ''
            pltEl = ''
            pltAz = ''

        pltTitle = '%s, %s, %s\n' % (pltInstrument, pltModelCase, pltPointing)
        if args.iPointingCase == 0:
            pltTitle += '%s, %s\n' % (os.path.split(args.StringKernelMetaFile)[1],
                                      pltUtcStartTime)
        else:
            pltTitle += 'R: %i km, PA: %i, LAT: %i, a: %i, b: %i, c: %i\n' %\
                        (args.UserR / 1000, args.UserPhaseAngle,
                         args.UserLatitude, args.UserAlpha,
                         args.UserBeta, args.UserGamma)

        if args.iModelCase == 0:
            pltTitle += '%s\n' % (os.path.split(args.StringDataFileDSMC)[1])
        elif args.iModelCase == 1:
            pltTitle += 'Q: %.2e [#/s],  v: %i [m/s],  Tp: %.0e [s]' % \
                        (args.QHaser, args.vHaser, args.tpHaser)

            if args.tdHaser is None:
                pltTitle += '\n'
            else:
                pltTitle += ', Td: %0.e [s]\n' % (args.tdHaser)
        elif args.iModelCase == 2:
            pass

    elif measurement == 'insitu':
        pltTitle = 'ICES in-situ tool\n'
        if args.iModelCase == 0:
            pltTitle += 'model: DSMC, case: %s\n' % (os.path.split(args.StringDataFileDSMC)[0].split('/')[-1])
            pltTitle += 'spice: %s' % (os.path.split(args.StringKernelMetaFile)[1].split('.')[0])
        elif args.iModelCase == 1:
            pltTitle += 'model: HASER, spice: %s\n' % (os.path.split(args.StringKernelMetaFile)[1].split('.')[0])
            pltTitle += 'Q = %.2e [#/s], v = %.0f [m/s], tp = %.2e [s]' % (args.QHaser, args.vHaser, args.tpHaser)

    return pltTitle


def create_plot_LOS_1d(args, ccd, pltTitle):

    plt.figure()
    if args.iInstrumentSelector == 3:           # alice
        plt.plot(range(5, 24), ccd, '-ok', linewidth=2)
        plt.grid(True)
        plt.xlabel('Pixel Number')
        plt.ylabel('Flux [photons / m2 / s]')
        plt.xticks(range(5, 24))
        plt.xlim((5, 23))

    plt.title(pltTitle)

    if args.DoShowPlots:
        plt.show()


def create_plot_LOS_2d(args, ccd, pltTitle):

    plt.figure()
    if args.iInstrumentSelector == 1:
        N = 8
        phi = 12
    elif args.iInstrumentSelector == 2:
        N = 8
        phi = 6
    elif args.iInstrumentSelector == 5:
        N = 6
        phi = 3.7

    dphi = phi/N

    xxticks = np.arange(-phi/2, phi/2 + dphi, dphi)
    yyticks = np.arange(-phi/2, phi/2 + dphi, dphi)

    plt.contourf(np.log10(ccd+0.1), 200)
    plt.colorbar(label='Column density [#/m2]')
    plt.xlabel("Instrument y axis [deg]")
    plt.ylabel("Instrument x axis [deg]")
    plt.xticks(np.arange(N+1)/N * (len(ccd[0])-1), np.round(xxticks,2))
    plt.yticks(np.arange(N+1)/N*(len(ccd[0])-1), np.round(yyticks,2))

    plt.title(pltTitle)
    if args.DoShowPlots:
        plt.show()

