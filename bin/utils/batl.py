import sys
import numpy as np

sys.path.append("../../../Models/LoS/BATL/srcReadAmr/")
import read_amr_wrapper

def interpolate(batsFileName, coords):
    batl = read_amr_wrapper.ReadBATL("../../../Models/LoS/BATL/lib/libWRAPAMR.so")
    batl.load_file(batsFileName)
    minDomain, maxDomain = batl.domain_limits()
    print " - BATS: minDomain: ", minDomain / 1000, " km"
    print " - BATS: maxDomain: ", maxDomain / 1000, " km"

    # limit coordinates to be well within the domain boundaries
    xMin = minDomain[0] * 0.95
    xMax = maxDomain[0] * 0.95

    # get number of variables, variable names and units from datafile.
    nVars = batl.get_nVar()
    varNames = batl.varnames()
     
    # do interpolation at coordinates
    S, lFound = batl.get_data_array(coords)

    batl.clean()

    return S, varNames[:nVars]
