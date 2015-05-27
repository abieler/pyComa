using HDF5, JLD
include("types2.jl")

fileName = ARGS[1]
pFileName = ARGS[2] * "/rosettaCoords.txt"
filePath = dirname(ARGS[1])

# subtract extension and last dot from filename
fileNameExtension = split(basename(fileName), ".")[end]
fileNameBase = basename(fileName)[1:end-(length(fileNameExtension)+1)]
println("- (julia) dsmc-run: ", joinpath(filePath, fileNameBase * ".h5"))

nCells = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nCells")
nBlocks = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nBlocks")
nCellsPerBlock = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nCellsPerBlock")
println(" - (julia) loading data")
@time nodes = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nodes")
@time nodeCoordinates = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nodeCoordinates")
@time cubeIndices = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/cubeIndices")
@time numberDensity= h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/numberDensity")


const xMax = maximum(nodes[:,:,1])
const yMax = maximum(nodes[:,:,2])
const zMax = maximum(nodes[:,:,3])

const halfSize = [xMax, yMax, zMax]
const root = [0.0, 0.0, 0.0]

println(" - (julia) building octree")
@time cellList = build_cells(nodes, cubeIndices, numberDensity, nCells)
@time blocks = build_blocks(nBlocks, cellList, nodes)
octree = Block(root, halfSize,1, initChildren, initCells)
@time populate_octree(octree, blocks, nBlocks)


@time n = doInSituCalculation(octree, pFileName)
println(" - (julia) write result to file")
println(" - (julia) size n: ", size(n))
writedlm(ARGS[2] * "/interpolation.out", n)
println(" - (julia) write result to file ok")
println(" - (julia) done julia interpolation")
println("+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +")
