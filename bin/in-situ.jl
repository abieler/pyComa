println("start julia in-situ calculation. not finished yet")
using HDF5
using JLD
include("types2.jl")


nDim = 3
path = ARGS[1]
pFileName = ARGS[2] * "/rosettaCoords.txt"

println("dsmc path: ", path)
DSMCfileNames = readdir(path)

fileName = ""
for f in DSMCfileNames
    println(f)
    if contains(f, "H2O")
        fileName = path * f
    end
end

halfSize = [200000., 200000., 200000.]
root = [0.0, 0.0, 0.0]

nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices, numberDensity = load_AMPS_data_full(fileName)
cellList, nodes = build_cells(nodeCoordinates, cubeIndices, nCells, numberDensity)
blocks = build_blocks(nBlocks, cellList, nodes)

octree = Block(root, halfSize,1, initChildren, initCells)
println("populating octree")
populate_octree(octree, blocks, nBlocks)
n = doInSituCalculation(octree, pFileName)
println(" - write result to file")
println(size(n))
writedlm(ARGS[2] * "/interpolation.out", n)
println(" - write result to file ok")
println(" - done julia interpolation")
println("+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +")
