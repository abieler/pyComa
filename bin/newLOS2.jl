using HDF5, JLD
include("types2.jl")

const fileName = ARGS[1]
#const pFileName = ARGS[2]
const pFileName = ARGS[2] * "/pointing.dat"
const filePath = dirname(ARGS[1])

# subtract extension and last dot from filename
fileNameExtension = split(basename(fileName), ".")[end]
fileNameBase = basename(fileName)[1:end-(length(fileNameExtension)+1)]

nCells = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nCells")
nBlocks = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nBlocks")
nCellsPerBlock = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nCellsPerBlock")
nodes = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nodes")
nodeCoordinates = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nodeCoordinates")
cubeIndices = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/cubeIndices")
numberDensity= h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/numberDensity")


const xMax = maximum(nodes[:,:,1])
const yMax = maximum(nodes[:,:,2])
const zMax = maximum(nodes[:,:,3])

const halfSize = [xMax, yMax, zMax]
const root = [0.0, 0.0, 0.0]

cellList = build_cells(nodes, cubeIndices, numberDensity, nCells)
blocks = build_blocks(nBlocks, cellList, nodes)
octree = Block(root, halfSize,1, initChildren, initCells)
println(" - populating octree")
populate_octree(octree, blocks, nBlocks)
@time n = doIntegration(octree, pFileName)
writedlm(ARGS[2]* "/ccd.dat", n)
