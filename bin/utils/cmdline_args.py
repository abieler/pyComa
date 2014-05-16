import argparse

parser = argparse.ArgumentParser()


def cmdline_args(parser):

    parser.add_argument("--iModelCase", type=int, choices=[0, 1, 2],
                        default=1, help='0: dsmc model, 1: haser model, 2: user model')
    parser.add_argument("--iPointingCase", type=int, choices=[0, 1, 2],
                        default=1, help='0: spice pointing, 1: user pointing, 2: user trajectory')
    parser.add_argument("--iInstrumentSelector", type=int, choices=[1, 2, 3, 4, 5, 6], default=3)
    parser.add_argument("--StringOutputDir", type=str, default='results')

    parser.add_argument("--StringDSMCdir", type=str)
    parser.add_argument("--StringDataFileDSMC", type=str)
    parser.add_argument("--IsDust", type=int, choices=[0, 1], help='1 for dust case, 0 for gas case')
    parser.add_argument("--DustSizeMin", type=float, default=0.0)   # dust particle radius in meters
    parser.add_argument('--DustSizeMax', type=float, default=2e-7)  # dust particle radius in meters

    parser.add_argument("--QHaser", type=float, default=1e27)    # production rate in molecules per second
    parser.add_argument("--vHaser", type=float, default=500.0)   # outflow velocity in m/s
    parser.add_argument("--tpHaser", type=float, default=1e5)    # exponential life time parent
    parser.add_argument("--tdHaser", type=float, default=0)      # exponential life time daughter

    parser.add_argument("--StringUserDataFile", type=str)        # file to upload from user which contains user coma model
    parser.add_argument("--DelimiterData", type=str)             # delimiter used in datafile
    parser.add_argument("--nHeaderRowsData", type=int)           # number of header lines in datafile
    parser.add_argument("--iDimUser", type=int)                  # number of dimenisons of user coma model
    parser.add_argument("--StringUserTrajectoryFile", type=str)  # file to upload from user which contains user coma model
    parser.add_argument("--DelimiterTraj", type=str)             # number of dimenisons of user coma model
    parser.add_argument("--nHeaderRowsTraj", type=int)           # number of dimenisons of user coma model

    parser.add_argument("--StringKernelMetaFile", type=str)      # spice .tm file that contains kernel info
    parser.add_argument("--StringUtcStartTime", type=str)        # format = yyyy-mm-ddTHH:MM:SS, e.g. 2014-07-11T00:00:00
    parser.add_argument("--StringUtcStopTime", type=str)
    parser.add_argument("--nDeltaT", type=int)                   # delta T in seconds for sampling along in situ trajectory

    parser.add_argument("--UserR", type=float, default=1e2)      # Distance in km from nucleus center
    parser.add_argument("--UserPhaseAngle", type=float, default=0.0)
    parser.add_argument("--UserLatitude", type=float, default=0.0)
    parser.add_argument("--UserAlpha", type=float, default=0.0)
    parser.add_argument("--UserBeta", type=float, default=0.0)
    parser.add_argument("--UserGamma", type=float, default=0.0)

    parser.add_argument("--StringHybridCase", type=str)
    parser.add_argument("--StringRuntimeDir", type=str)
    parser.add_argument("--StringMeasurement", type=str, default='LOS')

    parser.add_argument('--StringPlotting', type=str, default='matplotlib')  # matplotlib or bokeh selectable
    parser.add_argument('--DoShowPlots', type=str)

    parser.add_argument('--gFactor', type=float)
    parser.add_argument('--species', type=str, default='CO_')
    parser.add_argument('--gasTemp', type=float, default=100.0)
    parser.add_argument('--aliceDate', type=str, default='2014-07-11T23:23:00')

    args = parser.parse_args()

    return args
