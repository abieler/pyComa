using HDF5
using JLD
include("typesCleaned.jl")


nDim = 3
#fileName = ARGS[1] * "/Data_dsmc_15digits_update.dat"
#fileName = "/www/ices/Data/Coma/DSMC/CG_3.5_au_103/AMPS_1800.dat"
fileName = ARGS[1]
pFileName = ARGS[2] * "/pointing.dat"

halfSize = [200000., 200000., 200000.]
root = [0.0, 0.0, 0.0]

nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices, numberDensity = load_AMPS_data_full(fileName)
cellList, nodes = build_cells(nodeCoordinates, cubeIndices, nCells, numberDensity)
blocks = build_blocks(nBlocks, cellList, nodes)

octree = Block(root, halfSize,1, initChildren, initCells)
println("populating octree")
populate_octree(octree, blocks, nBlocks)
n = doIntegration(octree, pFileName)
writedlm(ARGS[2]* "/ccd.dat", n)
