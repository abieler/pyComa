using DataFrames

type Octree
  origin::Array{Float64,1}
  halfSize::Array{Float64,1}
  isLeaf::Int64
  children::Vector{Any}
  block::Any
  
  function Octree(origin::Array{Float64,1}, halfSize::Array{Float64,1},
		  isLeaf::Int64, block::Any)
		  
	  children = [None, None, None, None, None, None, None, None, None]
	  new(origin, halfSize, isLeaf, children, block)
  end
end

type Cell
  origin::Array{Float64,1}
  halfSize::Array{Float64,1}
  nodes::Array{Float64,2}
  volume::Float64
  densities::Array{Float64,1}
end

function create_ray(rRay::Array{Float64,1}, p::Array{Float64,1})
  
  distance = sqrt(sum(rRay.^2))
  rTravel = Array{Float64, 1}
  dr = 0.0
  
  rMax = 4.0*10^5
  rMin = 500.0
  
  while (distance < rMAx && distance > rMin)
    dr = distance / 10.0
  
  end


end
function loadPointsToIterate(fileName::String)
  df = readtable(fileName, skipstart=0, separator=',', header=false)
  nCoordinates = size(df)[1]
  sc_points = zeros(nCoordinates,3)
  for i=1:nCoordinates
    for j=1:3
      sc_points[i,j] = df[i,j]
    end
  end
  
  println("scp: ", sc_points[1,:])
  return sc_points
end

function findCellInBlock(block::Any, point::Array{Float64, 1})
  
  
  x = point[1] - block.cells[1].nodes[1,1]
  y = point[2] - block.cells[1].nodes[1,2]
  z = point[3] - block.cells[1].nodes[1,3]
  lx, ly, lz = block.halfSize / 2.5
  
  #println("x, y, z : ", x, " ", y, " ", z)
  #println("lx,ly,lz: ", lx, " ", ly, " ", lz)
  
  cellIndex = 1
  cellIndex += div(x, lx)
  cellIndex += div(y, ly) * 5
  cellIndex += div(z, lz) * 25
  
  #println("CellIndex: ", int(cellIndex))
  return int(cellIndex)

end


function triLinearInterpolation(cell::Any, point::Array{Float64,1})
  
  xd = (point[1] - cell.nodes[1,1]) / (cell.nodes[2,1] - cell.nodes[1,1])
  yd = (point[2] - cell.nodes[1,2]) / (cell.nodes[3,2] - cell.nodes[1,2])
  zd = (point[3] - cell.nodes[1,3]) / (cell.nodes[5,3] - cell.nodes[1,3])
  
  c00 = cell.densities[1] * (1-xd) + cell.densities[2] * xd
  c10 = cell.densities[5] * (1-xd) + cell.densities[6] * xd
  c01 = cell.densities[4] * (1-xd) + cell.densities[3] * xd
  c11 = cell.densities[8] * (1-xd) + cell.densities[7] * xd
  
  c0 = c00*(1-zd) + c10*zd
  c1 = c01*(1-zd) + c11*zd
  
  c = c0*(1-yd) + c1*yd
  
  #println("xd: ", xd)
  #println("yd: ", yd)
  #println("zd: ", zd)
  
  return c
end


function build_cells(nodeCoordinates::DataFrame, cubeIndices::DataFrame, nCells::Int64, numberDensity::DataArray{Float64, 1})

  nodes = Array(Float64,(nCells,8,3))
  nodeDensities = Array(Float64, (nCells, 8))
  origin = Array(Float64, (nCells,3))
  lCell = Array(Float64, 3)
  halfSize = Array(Float64, 3)
  allCells = Array(Cell, nCells)
  lx = 0.0
  ly = 0.0
  lz = 0.0

  print(" - building cells...              ")
  for i = 1:nCells
    for j = 1:8
      for k in 1:3
	nodes[i,j,k] = nodeCoordinates[cubeIndices[i,j],k]
      end
      nodeDensities[i,j] = numberDensity[cubeIndices[i,j]]
    end
    volume = lx * ly * lz
    for k=1:3
      lCell[k] = maximum(nodes[i,:,k]) - minimum(nodes[i,:,k])
      halfSize[k] = lCell[k] / 2.0 / 5.0
      origin[i,k] = nodes[i,1,k] + halfSize[k]
    end
    volume = lCell[1] * lCell[2] * lCell[3]
    allCells[i] = Cell(vec(origin[i,:]),halfSize,reshape(nodes[i,:,:], (8,3)),volume, vec(nodeDensities[i,:]))
  end

  println("OK")
  return allCells, nodes
 
end



type Block
  origin::Array{Float64, 1}
  halfSize::Array{Float64, 1}
  cells::Array{Cell, 1}
end

function getOctantContainingPoint(point::Array{Float64,1}, node::Octree)
    oct::Int64 = 1
    if (point[1] >= node.origin[1])
      oct += 4
    end
    if (point[2] >= node.origin[2])
      oct += 2
    end
    if (point[3] >= node.origin[3])
      oct += 1
    end
  return oct
end

function insertNode(point::Array{Float64,1}, node::Octree, block::Block)
  if (node.isLeaf == 1 && point != node.origin)
    node.isLeaf = 0
    splitNode(node)
    oct = getOctantContainingPoint(point, node)
    insertNode(point, node.children[oct], block)
  elseif (node.isLeaf == 0)
    oct = getOctantContainingPoint(point,  node)
    insertNode(point, node.children[oct], block)
  else
    node.block = block
  end
end

function findNodeContainingPoint(point::Array{Float64,1}, node::Octree)
  if (node.isLeaf == 0)
    oct = getOctantContainingPoint(point, node)
    findNodeContainingPoint(point, node.children[oct])
  elseif (node.isLeaf == 1)
    return node.block
  end
end

function splitNode(node::Octree)
  
  xc1 = [node.origin[1] - node.halfSize[1]/2, node.origin[2] - node.halfSize[2]/2,  node.origin[3] - node.halfSize[3]/2] 
  xc2 = [node.origin[1] - node.halfSize[1]/2, node.origin[2] - node.halfSize[2]/2,  node.origin[3] + node.halfSize[3]/2]
  xc3 = [node.origin[1] - node.halfSize[1]/2, node.origin[2] + node.halfSize[2]/2,  node.origin[3] - node.halfSize[3]/2] 
  xc4 = [node.origin[1] - node.halfSize[1]/2, node.origin[2] + node.halfSize[2]/2,  node.origin[3] + node.halfSize[3]/2]
  xc5 = [node.origin[1] + node.halfSize[1]/2, node.origin[2] - node.halfSize[2]/2,  node.origin[3] - node.halfSize[3]/2]
  xc6 = [node.origin[1] + node.halfSize[1]/2, node.origin[2] - node.halfSize[2]/2,  node.origin[3] + node.halfSize[3]/2]
  xc7 = [node.origin[1] + node.halfSize[1]/2, node.origin[2] + node.halfSize[2]/2,  node.origin[3] - node.halfSize[3]/2]
  xc8 = [node.origin[1] + node.halfSize[1]/2, node.origin[2] + node.halfSize[2]/2,  node.origin[3] + node.halfSize[3]/2]
  
  node.children[1] = Octree(xc1, node.halfSize/2, 1, None)
  node.children[2] = Octree(xc2, node.halfSize/2, 1, None)
  node.children[3] = Octree(xc3, node.halfSize/2, 1, None)
  node.children[4] = Octree(xc4, node.halfSize/2, 1, None)
  node.children[5] = Octree(xc5, node.halfSize/2, 1, None)
  node.children[6] = Octree(xc6, node.halfSize/2, 1, None)
  node.children[7] = Octree(xc7, node.halfSize/2, 1, None)
  node.children[8] = Octree(xc8, node.halfSize/2, 1, None)
end


function load_AMPS_densities(fileName::String)
  df = readtable(fileName, skipstart=2, separator=' ', nrows=1265221, header=false)
  numberDensity = df[:,7]
  return numberDensity
end


function load_AMPS_data(fileName::String)
  print(" - loading AMPS data              ")
  nodeCoordinates = readtable(fileName, skipstart=2, separator=' ', nrows=1265221, header=false)
  cubeIndices = readtable(fileName, skipstart=(2+1265221), separator=' ', header=false)
  println("OK")
  nCells = length(cubeIndices[:,1])
  nCellsPerBlock = 125
  nBlocks = int(nCells / nCellsPerBlock)
  println(" - nCells          : ", nCells)
  println(" - cells per block : ", nCellsPerBlock)
  println(" - nBlocks         : ", nBlocks)

  return nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices
end




function build_blocks(nBlocks::Int64, allCells::Array{Cell,1}, nodes::Array{Float64,3})

  print(" - assigning cells to blocks...   ")
  lBlock = Array(Float64, (nBlocks,3))
  origin = Array(Float64, (nBlocks,3))
  blocks = Array(Block, nBlocks)
  
  for i=1:nBlocks
    blockCells = allCells[(i-1)*nCellsPerBlock+1:i*nCellsPerBlock]
    blockCellCoords = nodes[(i-1)*nCellsPerBlock+1: i*nCellsPerBlock, :, :]
    for j=1:3
      lBlock[i,j] = maximum(blockCellCoords[:,:,j]) - minimum(blockCellCoords[:,:,j])  
      origin[i,j] = minimum(blockCellCoords[:,:,j]) + lBlock[i,j] / 2
    end
    blocks[i] = Block(vec(origin[i,:]), vec(lBlock[i,:]).*0.5, blockCells)
  end
  println("OK")
  return blocks
end

function populate_octree(oct::Octree, blocks::Array{Block, 1}, nBlocks::Int64)
  print(" - build octree...                ")
  for i=1:nBlocks
    insertNode(blocks[i].origin, oct, blocks[i])
  end
  println("OK")
end