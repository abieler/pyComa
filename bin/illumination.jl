#using HDF5
include("types2.jl")

pFileName = ARGS[1] * "/pointing.dat"

br = calcBrightnessFromNucleus(pFileName)
writedlm(ARGS[1]* "/ccd.dat", br)
