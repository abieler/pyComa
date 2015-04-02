using HDF5, JLD
include("types2.jl")

const fileName = ARGS[1]
#const pFileName = ARGS[2]
const pFileName = ARGS[2] * "/pointing.dat"
const filePath = dirname(ARGS[1])

# subtract extension and last dot from filename
fileNameExtension = split(basename(fileName), ".")[end]
fileNameBase = basename(fileName)[1:end-(length(fileNameExtension)+1)]
println("fileNameBase: ", fileNameBase)
@load joinpath(filePath,  fileNameBase * "_nodes.jld")
@load joinpath(filePath, fileNameBase * "_nodeCoordinates.jld")
@load joinpath(filePath, fileNameBase * "_cubeIndices.jld") 
@load joinpath(filePath, fileNameBase * "_numberDensity.jld")
@load joinpath(filePath, fileNameBase * "_nCells.jld")
@load joinpath(filePath, fileNameBase * "_nCellsPerBlock.jld")
@load joinpath(filePath, fileNameBase * "_nBlocks.jld")

const xMax = maximum(nodes[:,:,1])
const yMax = maximum(nodes[:,:,2])
const zMax = maximum(nodes[:,:,3])

const halfSize = [xMax, yMax, zMax]
println("halfSize: ", halfSize)
const root = [0.0, 0.0, 0.0]

@time cellList = build_cells(nodes, cubeIndices, numberDensity, nCells)
blocks = build_blocks(nBlocks, cellList, nodes)
octree = Block(root, halfSize,1, initChildren, initCells)
println("populating octree")
populate_octree(octree, blocks, nBlocks)
n = doIntegration(octree, pFileName)
writedlm(ARGS[2]* "/ccd.dat", n)
