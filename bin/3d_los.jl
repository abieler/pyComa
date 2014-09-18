include("types.jl")

ampsFileName = ARGS[1]
interpFileName = ARGS[1]

fileName = "/home/abieler/julia/octree/input/VolumeMesh_15digits.dat"
fileName2 = "/home/abieler/julia/octree/input/pic.H2O.s=0.out=6.dat"

numberDensity = load_AMPS_densities(fileName2)
nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices = load_AMPS_data(fileName)
cellList, nodes = build_cells(nodeCoordinates, cubeIndices, nCells, numberDensity)
blocks = build_blocks(nBlocks, cellList, nodes)

root = [0.0, 0.0, 0.0]
halfSize = [200000., 200000., 200000.]
isLeaf = 1
oct = Octree(root, halfSize, isLeaf, None)
populate_octree(oct, blocks, nBlocks)

if interpFileName == "input/sc_points.dat"
  sc_points = loadPointsToIterate(interpFileName)
  nSC_points = size(sc_points)[1]

  for i=1:nSC_points
    myBlock = findNodeContainingPoint(vec(sc_points[i,:]), oct)
    cellIndex = findCellInBlock(myBlock, vec(sc_points[i,:]))
    nIter = triLinearInterpolation(myBlock.cells[cellIndex], vec(sc_points[i,:]))
    println(nIter)
  end
end