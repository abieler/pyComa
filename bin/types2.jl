using DataFrames
using PyCall
@pyimport abieler.shapeUtils as utils

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


function load_AMPS_data(fileName::AbstractString)
  f = open(fileName, "r")
  nNodes = 0
  nCells = 0
  nCellsPerBlock = 125
  nBlocks = 0
  nHeaderRows = 2
  
  while !eof(f)
    line = readline(f)
    if ismatch(r"ZONE ", line)
      nNodes, nCells = [parse(Int, value) for value in matchall(r"(\d+)", line)]
      nBlocks = round(Int, nCells / nCellsPerBlock)
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
      x,y,z,n = [parse(Float64, value) for value in matchall(r"(-?\d.\d+[eE][+-]\d+)", line)[[1,2,3,4]]]
      nodeCoordinates[i-nHeaderRows,1] = x
      nodeCoordinates[i-nHeaderRows,2] = y
      nodeCoordinates[i-nHeaderRows,3] = z
      numberDensity[i-nHeaderRows] = n
    elseif (i > (nNodes+nHeaderRows))
      ijk = Int[parse(Int, value) for value in matchall(r"(\d+)", line)]
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
    @inbounds allCells[i] = Cell(vec(origin[i,1:3]), halfSize,reshape(nodes[i,1:8,1:3], (8,3)),volume, vec(nodeDensities[i,1:8]))
  end
  return allCells
end

function build_blocks(nBlocks::Int64, allCells::Array{Cell,1}, nodes::Array{Float64,3})
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
    @inbounds blocks[i] = Block(vec(origin[i,1:3]), vec(lBlock[i,1:3]).*0.5, 1, Array(Block, 8), blockCells)
  end
  if myid() == 1
      println("OK")
  end
  return blocks
end

function splitBlock(node::Block)
  
  xc1 = [node.origin[1] - node.halfSize[1]/2.0, node.origin[2] - node.halfSize[2]/2.0,  node.origin[3] - node.halfSize[3]/2.0] 
  xc2 = [node.origin[1] - node.halfSize[1]/2.0, node.origin[2] - node.halfSize[2]/2.0,  node.origin[3] + node.halfSize[3]/2.0]
  xc3 = [node.origin[1] - node.halfSize[1]/2.0, node.origin[2] + node.halfSize[2]/2.0,  node.origin[3] - node.halfSize[3]/2.0] 
  xc4 = [node.origin[1] - node.halfSize[1]/2.0, node.origin[2] + node.halfSize[2]/2.0,  node.origin[3] + node.halfSize[3]/2.0]
  xc5 = [node.origin[1] + node.halfSize[1]/2.0, node.origin[2] - node.halfSize[2]/2.0,  node.origin[3] - node.halfSize[3]/2.0]
  xc6 = [node.origin[1] + node.halfSize[1]/2.0, node.origin[2] - node.halfSize[2]/2.0,  node.origin[3] + node.halfSize[3]/2.0]
  xc7 = [node.origin[1] + node.halfSize[1]/2.0, node.origin[2] + node.halfSize[2]/2.0,  node.origin[3] - node.halfSize[3]/2.0]
  xc8 = [node.origin[1] + node.halfSize[1]/2.0, node.origin[2] + node.halfSize[2]/2.0,  node.origin[3] + node.halfSize[3]/2.0]
  
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
  nx = 5.0
  ny = 5.0
  nz = 5.0
  x = point[1] - block.cells[1].nodes[1,1]
  y = point[2] - block.cells[1].nodes[1,2]
  z = point[3] - block.cells[1].nodes[1,3]
  lx = block.halfSize[1] * 2.0 / nx 
  ly = block.halfSize[2] * 2.0 / ny
  lz = block.halfSize[3] * 2.0 / nz
  fx = fld(x, lx)
  fy = fld(y, ly)
  fz = fld(z, lz)

  if fx > (nx-1.0)
      fx = nx - 1.0
  end
  if fy > (ny-1.0)
      fy = ny-1.0
  end
  if fz > (nz-1.0)
      fz = nz-1.0
  end

  cellIndex = 1 + fx + fy*nx + fz*nx*ny
  return round(Int64, cellIndex)
end

function load_pointing_vectors(fileName::AbstractString)
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

function loadPointsToIterate(fileName::AbstractString)
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

function doIntegration(oct::Block, pFileName::AbstractString)
    println(" - start 3d los calculation")
    # pVector = pointing of vector, rRay = coordinates of origin of ray
    pVectors, rRay, nRays = load_pointing_vectors(pFileName)
    println(" - pVectors loaded")
    columnDensity::Float64 = 0.0
    distance::Float64 = 0.0
    # rMax::Float64 = 1.99*10^5
    rMax = minimum(oct.halfSize)*0.99
    println(" - new rMax determined: ", rMax)
    n_old::Float64 = 0.0
    l::Float64 = 0.0
    n = Array(Float64,nRays)
    dr = Array(Float64, 3)
    p = Array(Float64, 3)
    r = Array(Float64, 3)
    println("nRays: ", nRays)
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
        l = myBlock.halfSize[1]/3.0
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
        distance = sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])
      end

      n[i] = columnDensity
    end
    return n
end

function doIntegrationShadow(oct::Block, pFileName::AbstractString)
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


function doIntegrationOld(oct::Block, pFileName::AbstractString)
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

function doInSituCalculation(oct::Block, pFileName::AbstractString)
  rRay = loadPointsToIterate(pFileName)
  nPoints = size(rRay)[1]
  colDensity = zeros(Float64, nPoints)
  for i=1:nPoints
      r = vec(rRay[i,:])
      myBlock = findBlockContainingPoint(r, oct)
      cellIndex = findCellInBlock(myBlock, r)
      colDensity[i] = triLinearInterpolation(myBlock.cells[cellIndex], r)
  end

  return colDensity

end

function calcBrightnessFromNucleus(pFileName::AbstractString)

    path = dirname(pFileName)
    pVectors, rRay, nRays = load_pointing_vectors(pFileName)
    println(" - pointing vectors loaded")
    println("load shape model from new location")
    nTriangles, nodeCoords, triIndices, triangles, n_hat, p, triArea = utils.load_shape_model_vtk("/www/ices/Data/Nucleus/shapeModels/CSHP_87_14kCells.ply")
    #nTriangles, nodeCoords, triIndices, triangles, n_hat, p, triArea = utils.load_shape_model_vtk("/www/ices/Data/Nucleus/shapeModels/CSHP_87_35004Cells.ply")
    #nTriangles, nodeCoords, triIndices, triangles, n_hat, p, triArea = utils.load_shape_model_vtk("/www/ices/Data/Nucleus/shapeModels/SHAP5_stefano.ply")
    
    triangles2 = ones(3,3,nTriangles) 
    n_hat2 = ones(3, nTriangles)

    for i=1:nTriangles
      for j=1:3
        n_hat2[j,i] = n_hat[i,j]
        for k=1:3
          triangles2[k,j,i] = triangles[i,j,k]
        end
      end
    end


    f = open("rSun_hat.dat", "r")
    line = readline(f)
    println(line)
    d = matchall(r"(-?\d+.?\d+[eE]?[+-]?\d+)", line)
    println(d)
    rSun_hat = zeros(Float64, 3)    
    rRosetta = zeros(Float64,3)
    for i=1:3
      rSun_hat[i]= parse(Float64, d[i])
      rRosetta[i] = rRay[1,i]
    end

    br = Float64[]
    for i=1:nRays
        p = vec(pVectors[i,:])
        iFacet = visibleFacetIndex(triangles2, n_hat2, p, rRosetta, nTriangles)
        if iFacet > 0
           alpha = angle_btw_vec(vec(rSun_hat), vec(n_hat2[:,iFacet])) 
           push!(br, alpha)
        else
            push!(br, pi)
        end
    end
    return br
end

function angle_btw_vec(v1::Array{Float64,1}, v2::Array{Float64,1})
  if (v1 == v2)
      return 0.0
  end
  acos(dot(v1,v2))
end

function visibleFacetIndex(triangles::Array{Float64,3}, n_hat::Array{Float64, 2}, r::Array{Float64,1}, p0::Array{Float64,1}, nTriangles::Int64)
  i = 0
  j = 0
  k = 0
  
  a = 0.0
  b = 0.0
  rI = 0.0
  
  dot_uv = 0.0
  dot_uu = 0.0
  dot_vv = 0.0
  dot_wv = 0.0
  dot_wu = 0.0
  
  divisor = 0.0
  sI = 0.0
  tI = 0.0
  
  pI = [0.,0.,0.]
  u = [0.,0.,0.]
  v = [0.,0.,0.]
  w = [0.,0.,0.]
  
  
  for i=1:nTriangles
    a = 0.0
    b = 0.0
    @simd for k=1:3
      @inbounds a = a + n_hat[k,i] * (triangles[k,1,i] - p0[k]) 
      @inbounds b = b + n_hat[k,i] * r[k]
    end
    if ((a != 0.0) && (b != 0.0))
      rI = a / b
      if rI >= 0.0
        dot_uv = 0.0
        dot_uu = 0.0
        dot_vv = 0.0
        dot_wu = 0.0
        dot_wv = 0.0
        
        for k=1:3
	  @inbounds pI[k] = p0[k] + rI * r[k]
	  @inbounds u[k] = triangles[k,2,i] - triangles[k,1,i]
	  @inbounds v[k] = triangles[k,3,i] - triangles[k,1,i]
	  @inbounds w[k] = pI[k] - triangles[k,1,i]
	  
	  @inbounds dot_uv = dot_uv + u[k]*v[k]
	  @inbounds dot_uu = dot_uu + u[k]*u[k]
	  @inbounds dot_vv = dot_vv + v[k]*v[k]
	  @inbounds dot_wu = dot_wu + w[k]*u[k]
	  @inbounds dot_wv = dot_wv + w[k]*v[k]
        end
        
	divisor = dot_uv*dot_uv - dot_uu * dot_vv
	sI = (dot_uv*dot_wv - dot_vv*dot_wu) / divisor
	tI = (dot_uv*dot_wu - dot_uu*dot_wv) / divisor
	
	if ((tI >= 0.0) && (sI >= 0.0) && (sI + tI < 1.0))
	  return i
	end
    
      end
      
    end
  end
  return -1

end
