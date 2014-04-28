# -*- coding: utf-8 -*-
import os

iModelCase = 2
iPointingCase = 0
iInstrumentSelector = 1
StringOutputDir = '/Users/abieler/pyComa/results'
StringDataFileDSMC = '/Users/abieler/data/Coma/DSMC/CG_1.3_au_03/CG_1.3_au_03.Dust.dat'
isDust = 1
DustSizeMin = 0
DustSizeMax = 1e-7

QHaser = 1e27
vHaser = 500
tpHaser = 1e5
tdHaser = 0

StringKernelMetaFile = '/Users/abieler/PySPICE-master/cspice/kernels/metafiles/kernelsLTP_002H.tm'
StringUtcStartTime = '2014-11-24T00:00:00'
StringUtcStopTime = '2015-02-24T01:00:00'
nDeltaT = 600

UserR = 100e3
UserPhaseAngle = 0
UserLatitude = 0
UserAlpha = 2
UserBeta = 0
UserGamma = 0

StringUserDataFile = 'test_data/testUserFile_1d.txt'
StringUserTrajectoryFile = '/Users/abieler/pyComa/test_data/testUserTrajFile.txt'
UserDelimiter = ','
iUserNrOfHeaderRows = 0
iUserDim = 1

DelimiterTraj = 'space'
nHeaderRowsTraj = 0


cmd_string = '/Users/abieler/anaconda/bin/python in_situ_tool.py '

cmd_string += "--iModelCase=%i --iPointingCase=%i --iInstrumentSelector=%i --StringOutputDir=%s --StringDataFileDSMC=%s " %(iModelCase,
                                                                                                                          iPointingCase,
                                                                                                                          iInstrumentSelector,
                                                                                                                          StringOutputDir,
                                                                                                                          StringDataFileDSMC)

cmd_string += '--QHaser=%e --vHaser=%e --tpHaser=%e --tdHaser=%e ' % (QHaser, vHaser, tpHaser, tdHaser)
cmd_string += "--StringKernelMetaFile=%s --StringUtcStartTime=%s --StringUtcStopTime=%s --nDeltaT=%i " % (StringKernelMetaFile,
                                                                                                        StringUtcStartTime,
                                                                                                        StringUtcStopTime,
                                                                                                        nDeltaT)
cmd_string += "--UserR=%f --UserPhaseAngle=%f --UserLatitude=%f --UserAlpha=%f --UserBeta=%f --UserGamma=%f " %(UserR,
                                                                                                                UserPhaseAngle,
                                                                                                                UserLatitude,
                                                                                                                UserAlpha,
                                                                                                                UserBeta,
                                                                                                                UserGamma)
cmd_string += "--StringUserDataFile=%s --UserDelimiter=%s --iUserNrOfHeaderRows=%i --iUserDim=%i " % (StringUserDataFile,
                                                                                                   UserDelimiter,
                                                                                                   iUserNrOfHeaderRows,
                                                                                                   iUserDim)
#cmd_string += "--StringPlotting=bokeh"
cmd_string += "--DoShowPlots='y' --StringPlotting=matplotlib --StringUserTrajectoryFile=%s " % (StringUserTrajectoryFile)
cmd_string += "--DelimiterTraj=%s --nHeaderRowsTraj=%i --IsDust=%i --DustSizeMin=%e --DustSizeMax=%e" % (DelimiterTraj, nHeaderRowsTraj, isDust, DustSizeMin, DustSizeMax)

print cmd_string
os.system(cmd_string)
