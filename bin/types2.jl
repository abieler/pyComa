using DataFrames
#using HDF5

type Cell
  origin::Array{Float64,1}
  halfSize::Array{Float64,1}
  nodes::Array{Float64,2}
  volume::Float64
  densities::Array{Float64,1}
end

type Block
  origin::Array{Float64, 1}
  halfSize::Array{Float64, 1}
  isLeaf::Int64
  children::Array{Block, 1}
  cells::Array{Cell,1}
end

initChildren = Array(Block, 8) 
initCells = Array(Cell,1)


function load_AMPS_data(fileName::String)
  f = open(fileName, "r")
  nNodes = 0
  nCells = 0
  nCellsPerBlock = 125
  nBlocks = 0
  nHeaderRows = 2
  
  while !eof(f)
    line = readline(f)
    if ismatch(r"ZONE ", line)
      nNodes, nCells = [int(value) for value in matchall(r"(\d+)", line)]
      nBlocks = int(nCells / nCellsPerBlock)
      println("nNodes: ", nNodes)
      println("nCells: ", nCells)
      break
    end
  end
  close(f)
  nodeCoordinates = Array(Float64, nNodes, 3)
  numberDensity = Array(Float64, nNodes)
  cubeIndices = Array(Int64, nCells, 8)
  
  f = open(fileName, "r")
  i = 1
  while !eof(f)
    line = readline(f)
    if (i <= (nNodes+nHeaderRows)) && (i>nHeaderRows)
      x,y,z,n = [float(value) for value in matchall(r"(-?\d.\d+[eE][+-]\d+)", line)[[1,2,3,4]]]
      nodeCoordinates[i-nHeaderRows,1] = x
      nodeCoordinates[i-nHeaderRows,2] = y
      nodeCoordinates[i-nHeaderRows,3] = z
      numberDensity[i-nHeaderRows] = n
    elseif (i > (nNodes+nHeaderRows))
      ijk = [int(value) for value in matchall(r"(\d+)", line)]
      for index in 1:8
        cubeIndices[i-nNodes-nHeaderRows, index] = ijk[index]
      end
    end
    i+=1
  end
  close(f)
  return nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices, numberDensity

end


function build_nodes(nCells::Int64, nodeCoordinates::Array{Float64,2}, cubeIndices::Array{Int64,2})
  nodes = Array(Float64, (nCells,8,3))
  
  for i=1:nCells
    for j=1:8
       for k=1:3
         @inbounds nodes[i,j,k] = nodeCoordinates[cubeIndices[i,j],k]
       end
    end
  end
  return nodes
end

function build_cells(nodes::Array{Float64,3}, cubeIndices::Array{Int64,2}, numberDensity::Array{Float64,1}, nCells::Int64)
  nodeDensities = Array(Float64, (nCells, 8))
  origin = Array(Float64, (nCells,3))
  lCell = Array(Float64, 3)
  halfSize = Array(Float64, 3)
  allCells = Array(Cell, nCells)
  
  const lx = 0.0
  const ly = 0.0
  const lz = 0.0
  for i = 1:nCells
    for j = 1:8
      @inbounds nodeDensities[i,j] = numberDensity[cubeIndices[i,j]]
    end
    volume = lx * ly * lz
    for k=1:3
      @inbounds lCell[k] = maximum(nodes[i,1:8,k]) - minimum(nodes[i,1:8,k])
      @inbounds halfSize[k] = lCell[k] / 2.0 / 5.0
      @inbounds origin[i,k] = nodes[i,1,k] + halfSize[k]
    end
    volume = lCell[1] * lCell[2] * lCell[3]
    @inbounds allCells[i] = Cell(vec(origin[i,1:3]),halfSize,reshape(nodes[i,1:8,1:3], (8,3)),volume, vec(nodeDensities[i,1:8]))
  end
  return allCells
end

function build_blocks(nBlocks::Int64, allCells::Array{Cell,1}, nodes::Array{Float64,3})
  if myid() == 1
      print(" - assigning cells to blocks...   ")
  end
  lBlock = Array(Float64, (nBlocks,3))
  origin = Array(Float64, (nBlocks,3))
  blocks = Array(Block, nBlocks)
  
  for i=1:nBlocks
    blockCells = allCells[(i-1)*nCellsPerBlock+1:i*nCellsPerBlock]
    blockCellCoords = nodes[(i-1)*nCellsPerBlock+1: i*nCellsPerBlock, 1:8, 1:3]
    for j=1:3
      lBlock[i,j] = maximum(blockCellCoords[1:nCellsPerBlock,1:8,j]) - minimum(blockCellCoords[1:nCellsPerBlock,1:8,j])  
      origin[i,j] = minimum(blockCellCoords[1:nCellsPerBlock,1:8,j]) + lBlock[i,j] / 2
    end
    @inbounds blocks[i] = Block(vec(origin[i,1:3]), vec(lBlock[i,1:3]).*0.5, 1, initChildren, blockCells)
  end
  if myid() == 1
      println("OK")
  end
  return blocks
end

function splitBlock(node::Block)
  
  xc1 = [node.origin[1] - node.halfSize[1]/2, node.origin[2] - node.halfSize[2]/2,  node.origin[3] - node.halfSize[3]/2] 
  xc2 = [node.origin[1] - node.halfSize[1]/2, node.origin[2] - node.halfSize[2]/2,  node.origin[3] + node.halfSize[3]/2]
  xc3 = [node.origin[1] - node.halfSize[1]/2, node.origin[2] + node.halfSize[2]/2,  node.origin[3] - node.halfSize[3]/2] 
  xc4 = [node.origin[1] - node.halfSize[1]/2, node.origin[2] + node.halfSize[2]/2,  node.origin[3] + node.halfSize[3]/2]
  xc5 = [node.origin[1] + node.halfSize[1]/2, node.origin[2] - node.halfSize[2]/2,  node.origin[3] - node.halfSize[3]/2]
  xc6 = [node.origin[1] + node.halfSize[1]/2, node.origin[2] - node.halfSize[2]/2,  node.origin[3] + node.halfSize[3]/2]
  xc7 = [node.origin[1] + node.halfSize[1]/2, node.origin[2] + node.halfSize[2]/2,  node.origin[3] - node.halfSize[3]/2]
  xc8 = [node.origin[1] + node.halfSize[1]/2, node.origin[2] + node.halfSize[2]/2,  node.origin[3] + node.halfSize[3]/2]
  
  node.children[1] = Block(xc1, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1))
  node.children[2] = Block(xc2, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1))
  node.children[3] = Block(xc3, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1))
  node.children[4] = Block(xc4, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1))
  node.children[5] = Block(xc5, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1))
  node.children[6] = Block(xc6, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1))
  node.children[7] = Block(xc7, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1))
  node.children[8] = Block(xc8, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1))
end

function getOctantContainingPoint(point::Array{Float64,1}, block::Block)
    octant::Int64 = 1
    if (point[1] >= block.origin[1])
      octant += 4
    end
    if (point[2] >= block.origin[2])
      octant += 2
    end
    if (point[3] >= block.origin[3])
      octant += 1
    end
  return octant
end

function insertChild(point::Array{Float64,1}, parent::Block, child::Block)
  if (parent.isLeaf == 1 && point != parent.origin)
    parent.isLeaf = 0
    splitBlock(parent)
    octant = getOctantContainingPoint(point, parent)
    insertChild(point, parent.children[octant], child)
  elseif (parent.isLeaf == 0)
    octant = getOctantContainingPoint(point,  parent)
    insertChild(point, parent.children[octant], child)
  else
    parent.cells = child.cells
  end
end


function populate_octree(oct::Block, blocks::Array{Block, 1}, nBlocks::Int64)
  for i=1:nBlocks
    insertChild(blocks[i].origin, oct, blocks[i])
  end
  return oct
end


function findBlockContainingPoint(point::Array{Float64,1}, block::Block)
  if (block.isLeaf == 0)
    oct = getOctantContainingPoint(point, block)
    findBlockContainingPoint(point, block.children[oct])
  elseif (block.isLeaf == 1)
    return block
  end
end

function findCellInBlock(block::Block, point::Array{Float64, 1})
  
  x = point[1] - block.cells[1].nodes[1,1]
  y = point[2] - block.cells[1].nodes[1,2]
  z = point[3] - block.cells[1].nodes[1,3]
  lx = block.halfSize[1] / 2.5
  ly = block.halfSize[2] / 2.5
  lz = block.halfSize[3] / 2.5

  cellIndex = 1.0 + div(x, lx) + div(y, ly) * 5.0 + div(z, lz) * 25.0
  
  return int(cellIndex)

end


function load_pointing_vectors(fileName::String)
  df = readtable(fileName, skipstart=0, separator=',', header=false)
  const nVectors = size(df)[1]
  p = zeros(nVectors,3)
  rRay = zeros(nVectors,3)
  
  for i=1:nVectors
    for j=1:3
      p[i,j] = df[i,j]
      rRay[i,j] = df[i,j+3]
    end
  end

  return p, rRay, nVectors
end

function loadPointsToIterate(fileName::String)
  df = readtable(fileName, skipstart=0, separator=',', header=false)
  const nCoordinates = size(df)[1]
  sc_points = zeros(nCoordinates,3)
  for i=1:nCoordinates
    for j=1:3
      sc_points[i,j] = df[i,j]
    end
  end
  return sc_points
end

function triLinearInterpolation(cell::Cell, point::Array{Float64,1})
  
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
  
  return c
end


function doIntegrationDEBUG(oct::Block, pFileName::String)
    pVectors, rRay, nRays = load_pointing_vectors(pFileName)
    columnDensity::Float64 = 0.0
    distance::Float64 = 0.0
    const rMax::Float64 = 2.0*10^5
    n_old::Float64=0.0
    n = Array(Float64, nRays)
    dr = Array(Float64, 3)
    p = Array(Float64, 3)
    r = Array(Float64, 3)

    # 1 == true, 0 == false
    doCheckShadowing = 0
    if (doCheckShadowing == 1)
        println("shadowing of gas on")
        #rSun_hat = vec(readdlm("rSun_hat.dat", ','))
        #nTriangles, nodeCoords, triIndices, triangles, n_hat, p, triArea = utils.load_shape_model_vtk("/Users/abieler/meshes/CSHP_87_14kCells.ply")
    end

    for i=1:nRays
      columnDensity = 0.0
      n_old = 0.0
      for k=1:3
        p[k] = pVectors[i,k]
        r[k] = rRay[i,k]
      end
      distance = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
      while distance < rMax
        myBlock = findBlockContainingPoint(r, oct)
        l = myBlock.halfSize[1]/3
        cellIndex = findCellInBlock(myBlock, r)
        nIter = triLinearInterpolation(myBlock.cells[cellIndex], r)
        columnDensity += (n_old + nIter)/2.0 * l
        if nIter == 0.0
          break
        end
        for k=1:3
          dr[k] = p[k] * l
          r[k] = r[k] + dr[k]
        end
        distance = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
      end
      n[i] = columnDensity
    end
    return n
end

function doIntegration(oct::Block, pFileName::String)
    println(" - start 3d los calculation")
    # pVector = pointing of vector, rRay = coordinates of origin of ray
    pVectors, rRay, nRays = load_pointing_vectors(pFileName)
    columnDensity::Float64 = 0.0
    distance::Float64 = 0.0
    const rMax::Float64 = 2.0*10^5
    n_old::Float64 = 0.0
    n = Array(Float64,nRays)
    dr = Array(Float64, 3)
    p = Array(Float64, 3)
    r = Array(Float64, 3)
    
    for i=1:nRays
      columnDensity = 0.0
      n_old = 0.0
      for k=1:3
        p[k] = pVectors[i,k]
        r[k] = rRay[i,k]
      end
      distance = sqrt(r[1]^2 + r[2]^2 + r[3]^2)

      while distance < rMax
        myBlock = findBlockContainingPoint(r, oct)
        l = myBlock.halfSize[1]/3
        cellIndex = findCellInBlock(myBlock, r)
        nIter = triLinearInterpolation(myBlock.cells[cellIndex], r)
        if nIter == 0.0
          break
        end
        columnDensity += (n_old + nIter)/2.0 * l
        n_old = nIter * 1.0
        for k=1:3
          dr[k] = p[k] * l
          r[k] = r[k] + dr[k]
        end
        distance = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
      end
      n[i] = columnDensity
      #println("columnDensity: ", columnDensity)
    end
    return n
end

function doIntegrationShadow(oct::Block, pFileName::String)
    println(" - start 3d los calculation")
    # pVector = pointing of vector, rRay = coordinates of origin of ray
    pVectors, rRay, nRays = load_pointing_vectors(pFileName)
    columnDensity::Float64 = 0.0
    distance::Float64 = 0.0
    const rMax::Float64 = 2.0*10^5
    n_old::Float64 = 0.0
    n = Array(Float64,nRays)
    dr = Array(Float64, 3)
    p = Array(Float64, 3)
    r = Array(Float64, 3)
    isVisible = 1
    println("shadowing of gas on")
    rSun_hat = vec(readdlm("rSun_hat.dat", ','))
    nTriangles, nodeCoords, triIndices, triangles, n_hat, p, triArea = utils.load_shape_model_vtk("/Users/abieler/meshes/CSHP_87_14kCells.ply")
    for i=1:nRays
      columnDensity = 0.0
      n_old = 0.0
      for k=1:3
        p[k] = pVectors[i,k]
        r[k] = rRay[i,k]
      end
      distance = sqrt(r[1]^2 + r[2]^2 + r[3]^2)

      while distance < rMax
        myBlock = findBlockContainingPoint(r, oct)
        l = myBlock.halfSize[1]/3
        cellIndex = findCellInBlock(myBlock, r)
        nIter = triLinearInterpolation(myBlock.cells[cellIndex], r)
        if nIter == 0.0
          break
        end
        if (abs(r[1]) > 3000 && abs(r[2]) > 3000 && abs(r[3]) > 3000)
          isVisible = 1
        else
          isVisible = isSunlit(triangles, n_hat, rSun_hat, r, nTriangles, -1)
        end
        for k=1:3
          dr[k] = p[k] * l
          r[k] = r[k] + dr[k]
        end
        if n_old > 0.0
          if isVisible == 0
            columnDensity += 0.0
          end
          if isVisible == 1
            columnDensity += (n_old + nIter)/2.0 * l 
          end
        end
        n_old = nIter * 1.0
      end
      n[i] = columnDensity
    end
    return n
end


function doIntegrationOld(oct::Block, pFileName::String)
    println(" - start 3d los calculation")
    # pVector = pointing of vector, rRay = coordinates of origin of ray
    pVectors, rRay, nRays = load_pointing_vectors(pFileName)
    columnDensity::Float64 = 0.0
    distance::Float64 = 0.0
    rMax::Float64 = 2.0*10^5
    n_old::Float64=0.0
    n = zeros(nRays)
    
    # 1 == true, 0 == false
    doCheckShadowing = 0
    if doCheckShadowing == 1
        println("shadowing of gas on")
        rSun_hat = vec(readdlm("rSun_hat.dat", ','))
        nTriangles, nodeCoords, triIndices, triangles, n_hat, p, triArea = utils.load_shape_model_vtk("/Users/abieler/meshes/CSHP_87_14kCells.ply")
    end
    for i=1:nRays
      columnDensity = 0.0
      n_old = 0.0
      p = vec(pVectors[i,:])
      r = vec(rRay[i,:])
      distance = sqrt(sum(r.^2))

      while distance < rMax
        myBlock = findBlockContainingPoint(r, oct)
        l = myBlock.halfSize[1]/3
        cellIndex = findCellInBlock(myBlock, r)
        nIter = triLinearInterpolation(myBlock.cells[cellIndex], r)
        if nIter == 0.0
          break
        end
        # check if in shadow
        isVisible = 1 
        if doCheckShadowing == 1
          if (abs(r[1]) > 3000 && abs(r[2]) > 3000 && abs(r[3]) > 3000)
           isVisible = 1
          else
              isVisible = isSunlit(triangles, n_hat, rSun_hat, r, nTriangles, -1)
          end
        end
        dr = p * l
        r += dr
        if doCheckShadowing == 0
          #println("no shadow calculation of gas")
          if n_old > 0.0 && isVisible == 1
            columnDensity += (n_old + nIter)/2.0 * l
          end
        end
        if doCheckShadowing == 1
          if n_old > 0.0
            if isVisible == 0
              columnDensity += 0.0
            end
            if isVisible == 1
              columnDensity += (n_old + nIter)/2.0 * l 
            end
          end
        end
        n_old = nIter * 1.0
        #println(r, nIter)
        distance = sqrt(sum(r.^2))
      end
      n[i] = columnDensity
      #println("columnDensity: ", columnDensity)
    end
    return n
end


function doIntegrationParallel(I)

  halfSize = [200000., 200000., 200000.]
  root = [0.0, 0.0, 0.0]

  nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices = load_AMPS_data(fileName)
  numberDensity = load_AMPS_densities(fileName2)
  cellList, nodes = build_cells(nodeCoordinates, cubeIndices, nCells, numberDensity)
  blocks = build_blocks(nBlocks, cellList, nodes)
  
  
   xrange = I[1]
    if myid() == 2
        println(" - start 3d los calculation")
    end    
    println("worker: ", typeof(oct))
    pVectors,rRay, nRays = load_pointing_vectors(pFileName)
    nR = size(xrange)[1]
    columnDensity::Float64 = 0.0
    distance::Float64 = 0.0
    rMax::Float64 = 2.0*10^5
    n_old::Float64=0.0
    n = zeros(nR)
    n0 = xrange[1]
    for i=xrange
      columnDensity = 0.0
      n_old = 0.0
      p = vec(pVectors[i,:])
      r = vec(rRay[i,:])
      distance = sqrt(sum(r*r))
      
      while distance < rMax
	myBlock = findNodeContainingPoint(r, oct)
	l = myBlock.halfSize[1]/3.0
	cellIndex = findCellInBlock(myBlock, r)
	nIter = triLinearInterpolation(myBlock.cells[cellIndex], r)
	if nIter == 0.0
	  break
	end
	dr = p * l 
	r += dr 
	if n_old > 0.0
	  columnDensity += (n_old + nIter)/2.0 * l
	end
	n_old = nIter * 1.0
	distance = sqrt(sum(r*r))
      end
      n[i-(n0-1)] = columnDensity
    end
    println("finished 3d los calculation")
    return n

end