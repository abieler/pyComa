from __future__ import division
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib import rcParams

import bokeh.plotting as bplt
from bokeh.objects import Range1d

from data_loaders import load_in_situ_output, load_in_situ_output_hybrid

rcParams.update({'figure.autolayout': True})
font = {
        'family':'serif',
        'weight':'bold',
        'size':15}

matplotlib.rc('font', **font)

def plot_result_insitu_2(args, numberDensities_SC, species, dates_SC, r_SC):
    print 'plotting in-situ results...'
    nLinesPerFig = 6
    if args.StringHybridCase is None:
        pltTitle = build_plot_title(args, 'insitu')
    else:
        pltTitle = build_plot_title(args, 'hybrid')

    if args.StringPlotting.lower() == 'bokeh':
        create_plot_insitu_bokeh(args, numberDensities_SC, species, dates_SC, r_SC, nLinesPerFig, pltTitle)
    elif args.StringPlotting.lower() == 'matplotlib':
        create_plot_insitu_matplotlib(args, numberDensities_SC, species, dates_SC, r_SC, nLinesPerFig, pltTitle)

def plot_result_insitu(args):

    print 'plotting results...'
    path = args.StringOutputDir
    nLinesPerFig = 4
    if args.StringHybridCase is None:
        pltTitle = build_plot_title(args, 'insitu')
    else:
        pltTitle = build_plot_title(args, 'hybrid')

    filenames = [path + '/' + filename for filename in os.listdir(path) if filename.split('.')[-1] == 'out']
    all_n_SC = []
    all_species = []
    for filename, i in zip(filenames, range(len(filenames))):
        species = os.path.split(filename)[1].split('.')[0]
        if args.StringHybridCase is None:
            dates_SC, r_SC, n_SC = load_in_situ_output(filename)
        else:
            dates_SC, r_SC, n_SC = load_in_situ_output_hybrid(filename)
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
            bplt.figure(x_axis_type='datetime', plot_width=650, plot_height=400, title='in-situ tool')
            bplt.hold()
        bplt.line(dates_SC, np.log10(all_n_SC[i]), legend=all_species[i], line_width=2,
                  line_color=pltColors[i % nLinesPerFig])

        bplt.xaxis().axis_label = 'Time'
        bplt.yaxis().axis_label = 'log10(n) [#/m3]'
        if (i % nLinesPerFig == (nLinesPerFig - 1)) or (i == len(all_species) - 1):
            bplt.grid().grid_line_alpha = 0.4
            bplt.figure(x_axis_type='datetime', plot_width=650, plot_height=400, title='in-situ tool')
            bplt.line(dates_SC, r_SC/1000, line_width=2, line_color='black')
            bplt.grid().grid_line_alpha = 0.4
            bplt.xaxis().axis_label = 'Date'
            bplt.yaxis().axis_label = 'r from comet center [km]'
            bplt.save()

    if args.DoShowPlots:
        bplt.show()
    else:
        print 'not showing plots on screen'


def create_plot_insitu_matplotlib(args, all_n_SC, all_species, dates_SC, r_SC, nLinesPerFig, pltTitle):

    for species, i in zip(all_species, range(len(all_species))):
            if (i % nLinesPerFig == 0):   # make new figure for every 4 species
                fig = plt.figure(figsize=(14, 14))
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
            ax1.semilogy(dates_SC, all_n_SC[i], label=all_species[i], lw=2)

            # plot distance from comet and save figure
            if (i % nLinesPerFig == (nLinesPerFig - 1)) or (i == len(all_species) - 1):
                nFig = (i-(nLinesPerFig-1)) // nLinesPerFig
                if nFig < 0:
                    nFig = 0
                ax1.set_title(pltTitle)
                ax1.set_ylabel('Number density [#/m3]')
                ax1.grid(True)
                ax1.legend(framealpha=1, loc=0)

                ax2.plot(dates_SC, r_SC / 1000, '-k', lw=2)
                ax2.set_ylabel('Distance from comet center [km]')
                ax2.grid(True)
                plt.gcf().autofmt_xdate()
                plt.savefig(args.StringOutputDir + '/' + 'result_%i.png' % nFig)

    if args.DoShowPlots:
        plt.show()
    else:
        print 'not showing plots on screen'


def plot_result_LOS(ccd, StringOutFileName, args, ccd_limits):

    pltTitle = build_plot_title(args, 'LOS')

    if args.StringPlotting == 'matplotlib':
        if args.iInstrumentSelector in [1, 2, 5, 9, 10]:
            create_plot_LOS_2d_matplotlib(args, ccd, pltTitle, StringOutFileName, ccd_limits)
        elif args.iInstrumentSelector in [3, 6, 7]:
            create_plot_LOS_1d_matplotlib(args, ccd, pltTitle, StringOutFileName)
        elif args.iInstrumentSelector in [4, 8]:
            print 'not generating plot for MIRO. Single pixel only.'

    elif args.StringPlotting == 'bokeh':
        if args.iInstrumentSelector in [1, 2, 5]:
            create_plot_LOS_2d_bokeh(args, ccd, pltTitle, StringOutFileName)
        elif args.iInstrumentSelector in [3, 6, 7]:
            create_plot_LOS_1d_bokeh(args, ccd, pltTitle, StringOutFileName)


def create_plot_LOS_1d_bokeh(args, ccd, pltTitle, StringOutFileName):
    print 'setting plot title to:'
    print pltTitle
    bplt.output_file(args.StringOutputDir + '/' + StringOutFileName, title=pltTitle)
    if args.iInstrumentSelector in [3, 7]:
        bplt.line(np.arange(5, len(ccd[0, :])+5), ccd[0, :], line_width=2, line_color='black')
        bplt.xaxis().axis_label='Pixel number'
        if args.iInstrumentSelector == 3:
            bplt.yaxis().axis_label='Column Density [#/m2]'
        else:
            bplt.yaxis().axis_label='photons / m2 / s'
    else:
        bplt.line(range(0, len(ccd[0, :])), ccd[0, :], line_width=2, line_color='black')
    bplt.save()

    if args.DoShowPlots:
        bplt.show()
    else:
        print 'bokeh: not showing plots on screen'


def create_plot_LOS_2d_bokeh(args, ccd, pltTitle, StringOutFileName):

    nPixels = len(ccd[:, 0])
    bplt.output_file(args.StringOutputDir + '/' + StringOutFileName)
    bplt.image(image=[np.log10(ccd)], x=[0], y=[0],
               dw=[nPixels], dh=[nPixels], palette=["Spectral-11"],
               x_range=Range1d(start=0, end=nPixels),
               y_range=Range1d(start=0, end=nPixels))
    print pltTitle
    bplt.curplot().title = pltTitle
    bplt.save()

    if args.DoShowPlots:
        bplt.show()
    else:
        print 'bokeh: not showing plots on screen'


def create_plot_LOS_1d_matplotlib(args, ccd, pltTitle, figName):

    plt.figure(figsize=(15, 12))
    plt.grid(True)
    plt.xlabel('Pixel Number')
    if args.iInstrumentSelector in [3, 7]:           # alice
        plt.plot(range(5, 24), ccd[0, :], '-ok', linewidth=2)
        if args.iInstrumentSelector == 3:
            plt.ylabel('Column Density [#/m2]')
        elif args.iInstrumentSelector == 7:
            plt.ylabel('Flux [photons / m2 / s]')
        plt.xticks(range(5, 24))
        plt.xlim((5, 23))
    else:
        plt.plot(ccd[0, :])
        plt.ylabel('Column Density [#/m2]')

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.title(pltTitle)
    plt.savefig(args.StringOutputDir + '/' + figName)
    if args.DoShowPlots:
        plt.show()


def create_plot_LOS_2d_matplotlib(args, ccd, pltTitle, figName, ccd_limits):

    plt.figure(figsize=(14, 12))
    if args.iInstrumentSelector == 1:
        N = 8
        phi = 12
    elif args.iInstrumentSelector == 2:
        N = 8
        phi = 6
    elif args.iInstrumentSelector == 5:
        N = 6
        phi = 3.7
    elif args.iInstrumentSelector == 10:
        N = 5
        phi = 10.0

    dphi = phi/N

    xxticks = np.arange(-phi/2, phi/2 + dphi, dphi)
    yyticks = np.arange(-phi/2, phi/2 + dphi, dphi)

    nLevels = 200
    dLevels = (np.log10(ccd_limits[1]) - np.log10(ccd_limits[0])) / nLevels
    pltLevels = np.arange(np.log10(ccd_limits[0]), np.log10(ccd_limits[1])+dLevels, dLevels)
    try:
        plt.contourf(np.log10(ccd+1.0e-200), levels=pltLevels)
    except:
        print '################################################'
        print 'plotting error, contour levels not well defined'
        print 'trying with auto-defined contours...'
        print '################################################'
        plt.contourf(np.log10(ccd+1.0e-200))

    if args.iInstrumentSelector == 10:
       plt.colorbar(label=r'log10 $S_{total} [W m^{-2} Hz^{-1}$]')
    else:
       plt.colorbar(label='log10 column density [#/m2]')

    plt.xlabel("Instrument y axis [deg]")
    plt.ylabel("Instrument x axis [deg]")
    plt.xticks(np.arange(N+1)/N * (len(ccd[0])-1), np.round(xxticks,2))
    plt.yticks(np.arange(N+1)/N*(len(ccd[0])-1), np.round(yyticks,2))

    plt.title(pltTitle)
    plt.savefig(args.StringOutputDir + '/' + figName)
    print 'saved result as', args.StringOutputDir + '/' + figName
    if args.DoShowPlots:
        plt.show()


def plot_miro(sTotal, aveBright, F, aDust, args):

    pltTitle = build_plot_title(args, 'LOS')

    # Make the first figure.  Fixed frequency.  Brightness vs. Dust Radius)
    figName = "MIRO_brightness_v_radius.png"
    plt.figure(figsize=(15, 12))
    plt.grid(True)
    plt.xlabel(r'Dust Radius [$\mu m$]')
    plt.ylabel(r'Average Brightness [$W m^{-2} Hz^{-1} sr^{-1}$]')
    plt.loglog() 
    plt.plot(aDust/1e-6, aveBright[0,0,:,0], 'or', markersize=15)
    plt.title(pltTitle)

    plt.savefig(args.StringOutputDir + '/' + figName)
    if args.DoShowPlots:
        plt.show()

    # # Make the second figure.  Fixed Radius.  Brightness vs. Frequency)
    # figName = "MIRO_brightness_v_freq.png"
    # plt.figure(figsize=(15, 12))
    # plt.grid(True)
    # plt.xlabel('Frequency')
    # plt.ylabel('Average Brightness [??]')
    # #plt.loglog() 
    # plt.plot(F, aveBright[0,0,5,:], 'or', markersize=20)
    # plt.title(pltTitle)
    # 
    # plt.savefig(args.StringOutputDir + '/' + figName)
    # if args.DoShowPlots:
    #     plt.show()

    # # Make the third figure.  Fixed Radius and freq.  varied pointing
    # cmap = cm.jet
    # figName = "MIRO_brightness_v_FOV.png"
    # plt.figure(figsize=(15, 12))
    # plt.grid(True)
    # plt.xlabel('Angle')
    # plt.ylabel('Angle')
    # x = np.arange(sTotal.shape[0])
    # y = np.arange(sTotal.shape[1])
    # plt.scatter(x, y, c=aveBright[:,:,5,1], cmap=cmap, marker='+', s=3, linewidths=30)
    # plt.title(pltTitle)
    # 
    # plt.savefig(args.StringOutputDir + '/' + figName)
    # if args.DoShowPlots:
    #     plt.show()


def build_plot_title(args, measurement='LOS'):

    if measurement == 'LOS':
        if args.iInstrumentSelector == 1:
            pltInstrument = 'OSIRIS WAC '
        elif args.iInstrumentSelector == 2:
            pltInstrument = 'OSIRIS NAC '
        elif args.iInstrumentSelector in [3, 7]:
            pltInstrument = 'ALICE '
        elif args.iInstrumentSelector in [4, 8, 10]:
            pltInstrument = 'MIRO '
        elif args.iInstrumentSelector == 5:
            pltInstrument = 'VIRTIS M '
        elif args.iInstrumentSelector == 6:
            pltInstrument = 'VIRTIS H '
        elif args.iInstrumentSelector == 9:
            pltInstrument = 'VIRTIS'

        pltTitle = 'ICES line of sight tool for '
        pltTitle += pltInstrument
        pltTitle += '\n'
        if args.iModelCase == 0:
            pltTitle += 'Coma model: DSMC (%s)\n' % (os.path.split(args.StringDataFileDSMC)[1])
        elif args.iModelCase == 1:
            pltTitle += 'Coma model: Haser (Q: %.2e[#/s],  v: %i[m/s],  Tp: %.0e[s]' % \
                        (args.QHaser, args.vHaser, args.tpHaser)
            if args.tdHaser == 0:
                pltTitle += ')\n'
            else:
                pltTitle += ', Td: %0.e[s])\n' % (args.tdHaser)
        elif args.iModelCase == 2:
            pltTitle += 'Coma model: user defined (%s)\n' % (os.path.split(args.StringUserDataFile)[1])
        elif args.iModelCase == 3:
            pltTitle += 'Dust analytic model:  (Q: %.2e [kg/s],  v: %i [m/s])\n' % \
                        (args.QDust, args.vDust)


        if args.iPointingCase == 0:
            pltTitle += 'Pointing: SPICE (%s, %s)\n' % (os.path.split(args.StringKernelMetaFile)[1],
                                                        args.StringUtcStartTime)
        else:
            pltTitle += 'Pointing: user defined (R: %i km, PA: %i, LAT: %i, a: %i, b: %i, c: %i)\n' %\
                        (args.UserR, args.UserPhaseAngle,
                         args.UserLatitude, args.UserAlpha,
                         args.UserBeta, args.UserGamma)

    elif measurement == 'insitu':
        pltTitle = 'ICES in-situ tool\n'
        if args.iModelCase == 0:
            pltTitle += 'Coma model: DSMC '
            if args.IsDust:
                pltTitle += 'Dust, rmin=%.2e, rmax=%.2e ' %(args.DustSizeMin, args.DustSizeMax)
            pltTitle += '(%s)\n' % (os.path.split(args.StringDataFileDSMC)[0].split('/')[-1])
        elif args.iModelCase == 1:
            pltTitle += 'Coma model: HASER '
            pltTitle += 'Q = %.2e [#/s], v = %.0f [m/s], tp = %.2e [s]\n' % (args.QHaser,
                                                                             args.vHaser,
                                                                             args.tpHaser)
        elif args.iModelCase == 2:
            pltTitle += 'Coma model: user defined (%s)\n' % (os.path.split(args.StringUserDataFile)[1])

        if args.iPointingCase == 0:
            pltTitle += 'Trajectory: spice (%s)\n' % (os.path.split(args.StringKernelMetaFile)[1].split('.')[0])
        elif args.iPointingCase == 2:
            pltTitle += 'Trajectory: user defined (%s)\n' % (os.path.split(args.StringUserTrajectoryFile)[1])
        else:
            pltTitle += '\n'
    elif measurement == 'hybrid':
        pltTitle = 'ICES in-situ tool\n'
        pltTitle += 'Coma model: Hybrid AIKEF, (%s)\n' % args.StringHybridCase
        if args.iPointingCase == 0:
            pltTitle += 'Trajectory: spice (%s)' % (os.path.split(args.StringKernelMetaFile)[1].split('.')[0])
        elif args.iPointingCase == 2:
            pltTitle += 'Trajectory: user defined (%s)' % (os.path.split(args.StringUserTrajectoryFile)[1])
    return pltTitle
