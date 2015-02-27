using HDF5
using JLD
include("typesCleaned.jl")


pFileName = ARGS[1] * "/pointing.dat"

br = calcBrightnessFromNucleus(pFileName)
writedlm(ARGS[1]* "/ccd.dat", br)
