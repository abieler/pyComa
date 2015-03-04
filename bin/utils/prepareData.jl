using HDF5, JLD
include("../types2.jl")
#=
# load amps output files and save details thereof
# as binary .jld files (jld files are hdf5 compatible)

fileName = "data.dat"
fileName = ARGS[1]
nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices, numberDensity = load_AMPS_data(fileName)
nodes = build_nodes(nCells, nodeCoordinates, cubeIndices)

@save "nodes.jld" nodes
@save "nodeCoordinates.jld" nodeCoordinates
@save "cubeIndices.jld" cubeIndices
@save "numberDensity.jld" numberDensity
@save "nCells.jld" nCells
@save "nCellsPerBlock.jld" nCellsPerBlock
@save "nBlocks.jld" nBlocks
=#
