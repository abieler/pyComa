using HDF5, JLD
include("types2.jl")

# this script has to be applied for all new 3D cases with amps.
# it reads the original amps output file and saves some content
# of it into binary hdf5 files for quicker loading when doing
# the actual line of sight calculations with newLOS2.jl

# This script has to be run only once and needs the filename
# of the AMPS output as only input parameter.

fileNameExtension = split(basename(ARGS[1]), ".")[end]

fileName = ARGS[1]
filePath = dirname(ARGS[1])

# subtract extension and last dot from filename
fileNameBase = basename(ARGS[1])[1:end-(length(fileNameExtension)+1)]

#fileName = "data.dat"
nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices, numberDensity = load_AMPS_data(fileName)
nodes = build_nodes(nCells, nodeCoordinates, cubeIndices)

h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nCells", nCells)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nBlocks", nBlocks)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nCellsPerBlock", nCellsPerBlock)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nodes", nodes)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nodeCoordinates", nodeCoordinates)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/cubeIndices", cubeIndices)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/numberDensity", numberDensity)

