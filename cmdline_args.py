import argparse

parser = argparse.ArgumentParser()

def cmdline_args(parser):

    parser.add_argument("--iModelCase", type=int, choices=[0, 1, 2], default=1, help='0: dsmc model, 1: haser model, 2: user model')
    parser.add_argument("--iPointingCase", type=int, choices=[0, 1], default=1, help='0: spice pointing, 1: user pointing')
    parser.add_argument("--iInstrumentSelector", type=int, choices=[1, 2, 3, 4, 5, 6], default=1)
    parser.add_argument("--StringOutputDir", type=str, default='.')

    parser.add_argument("--StringDSMCdir", type=str)
    parser.add_argument("--StringDataFileDSMC", type=str)
    parser.add_argument("--IsDust", type=int, choices=[0, 1], help='1 for dust case, 0 for gas case')
    parser.add_argument("--DustSizeMin", type=float)
    parser.add_argument('--DustSizeMax', type=float)

    parser.add_argument("--QHaser", type=float, default=1e27)
    parser.add_argument("--vHaser", type=float, default=500.0)
    parser.add_argument("--tpHaser", type=float, default=1e5)
    parser.add_argument("--tdHaser", type=float, default=0)

    parser.add_argument("--StringUserDataFile", type=str)                          # file to upload from user which contains user coma model
    parser.add_argument("--UserDelimiter", type=str)                         # delimiter used in datafile
    parser.add_argument("--iUserNrOfHeaderRows", type=int)                   # number of header lines in datafile
    parser.add_argument("--iUserDim", type=int)                               # number of dimenisons of user coma model

    parser.add_argument("--StringKernelMetaFile", type=str)
    parser.add_argument("--StringUtcStartTime", type=str)
    parser.add_argument("--StringUtcStopTime", type=str)
    parser.add_argument("--nDeltaT", type=int)

    parser.add_argument("--UserR", type=float, default=1e5)                               # Distance in km from nucleus center
    parser.add_argument("--UserPhaseAngle", type=float, default=0.0)
    parser.add_argument("--UserLatitude", type=float, default=0.0)
    parser.add_argument("--UserAlpha", type=float, default=0.0)
    parser.add_argument("--UserBeta", type=float, default=0.0)
    parser.add_argument("--UserGamma", type=float, default=0.0)

    parser.add_argument("--StringHybridCase", type=str)
    parser.add_argument("--StringRuntimeDir", type=str)

    parser.add_argument('--StringPlotting', type=str, default='matplotlib')
    parser.add_argument('--DoShowPlots', type=str)

    args = parser.parse_args()

    return args
