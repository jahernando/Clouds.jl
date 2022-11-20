
#using Base
import StatsBase     as SB
import LinearAlgebra as LA
import Graphs        as GG

export moves, clouds

#-----------------
# Data Types
#-----------------

# sort-cuts type definitions
N   = Vararg
T2I = Tuple{Int64, Int64}
TNI = Tuple{N{Int64}}
TNV = Tuple{N{Float64}}
TNN = Tuple{N{<:Number}}
VI  = Vector{Int64}
VB  = Vector{Bool}
VF  = Vector{Float64}
VN  = Vector{<:Number}
VAN = Vector{T} where T <: Array{<:Number}
AN  = T where T <: Array{<:Number}
AI  = T where T <: Array{Int64}


# Helper Data Type to hold one step movements in a 2D or 3D grid
# one steps 1D movements are : [1, 0], [0, 1], [1, 1], [-1, 0], ...
struct Moves
	moves  ::Tuple{N{VI}}   # list with the vector of unitary movements in 2D or 3D
	i0     ::Int64          # index of the null movement (0., 0) or (0, 0, 0)
	isym   ::Tuple{N{T2I}}  # list of the pair of indices of symmetric movements i.e (1, 0), (-1, 0)
	iortho ::Dict{T2I, VI}  # dictionary with the indices of the ortogonal movements
	                        # i.e movements (0, 1) (0, -1) are orthogornal to (1, 0), (-1, 0)
end

# Data from cells (2d pixels or 3d voxels) of Clouds
struct DataCells
	coors    ::Tuple{N{VN}}   # (x,y) or (x,y, z) tuple with the cell cordinates
	contents ::VN             # cell contents
	cells    ::Vector{TNI}    # Vector with the cell indices, (i, j) or (i, j, k)
	grad     ::VN             # cell gradient
	igrad    ::VI             # cell index of the gradient move
	lap      ::VN             # cell laplacian
	curmax   ::VN             # cell maximum curvature
	icurmax  ::VI             # cell index of the move with the maximum curvature
	curmin   ::VN             # cell minimum curvature
	icurmin  ::VI             # cell index of the move with the minimum curvature
	node     ::VI             # node index which this cell belongs to
	cloud    ::VI             # cloud index which this cell belongs to
	nborders ::VI             # number of border cells whose belong to another node
end

# Data of Nodes (the elements of clouds, also the vertices of the graph)
struct DataNodes
	size        ::VI             # node size (number of cells into this node)
	contents    ::VN             # node contents (sum of the cells contents of this node)
	maxcontent  ::VN             # node maximum content (content of the cell with the maximum content)
	maxgrad     ::VN             # node maximum gradient
	maxlap      ::VN             # node maximum laplacian
	minlap      ::VN             # node minimum laplacian
	curmax      ::VN             # node maximum curvature
	curmin      ::VN             # node minimum curvature
	nedges      ::VI             # number of nodes connected to this node
	cloud       ::VI             # cloud index which this node belongs to
	coors       ::Tuple{N{VN}}   # coordenates of the node barycenter
	coors_std   ::Tuple{N{VN}}   # std of the node coordinates
	coors_cell  ::Tuple{N{VN}}   # coordenates of the cell with null gradient (top of the node)
	ecc         ::VI             # distance in nodes of this node to the end of the graph (eccenticity)
	extremes    ::VB             # true if the node is the pair extreme of the graph with the largest contents
end

#----------
# Moves
#----------

"""
	moves(ndim)

Returns a `Moves` DataType for 2 or 3 dimentions.

`Moves` holds data for the one step movements in one unit in any coordinate.

# Examples

```julia-repr
julia> mm = Moves(2)
julia> mm.moves
([-1, -1], [-1, 0], [-1, 1], [0, -1], [0, 0], [0, 1], [1, -1], [1, 0], [1, 1])
```
"""
function moves(ndim::Int64)

	if (ndim < 1) || (ndim > 3)
		throw(ArgumentError("moves dimensions valid are 2, 3"))
	end

	# moves
	moves = ndim == 2 ? Tuple([i, j] for i in -1:1:1 for j in -1:1:1) :
	     Tuple([i, j, k] for i in -1:1:1 for j in -1:1:1 for k in -1:1:1)

	# null move
	move0 = ndim == 2 ? [0, 0] : [0, 0, 0]
	kmove0 = [i for (i, move) in enumerate(moves) if (move == move0)][1]

	# pair of indices of symmetric moves
	smoves = Tuple((i, j) for (i, movei) in enumerate(moves)
	          for (j, movej) in enumerate(moves)
		      if (movei == -1 .* movej ) &  (i > j))

	#list of orthogonal indeces of moves orthogonal to a pair os symmetric indices of moves
	omoves = Dict{T2I, VI}()
	for ii in smoves
		movei = moves[ii[1]]
		omoves[ii] = [j for (j, movej) in enumerate(moves)
		    if ((sum(movej .* movei) == 0.0) & (sum(movej .* movej) != 0.0))]
	end

	#return (moves = moves, i0 = kmove0, isym = smoves, iortho = omoves)
	return Moves(moves, kmove0, smoves, omoves)
end


#-----
# Clouds
#-----

"""
	clouds(coors, energy, steps; cellmode = 1, threshold = 0.)

From a collection of points (`coors`) with weights (`energy`),
the Clouds algorithm creates a pixel or voxelixed space.
Every pixel or voxel is called a *cell*.
The energy (content) on every cell is compute. Only cells above a threshold
are considered. For every cell the algoritm computes the local gradient, laplacian,
the maximum and minimum curvature.

Cells are grouped into nodes. A *node* can be constructed from one cell if `cellnode = true`,
otherwise cells whose gradient ends in the same cell (which is the cell with
the highest contents or energy) are grouped into a node. The maximum gradient,
laplaciant, curvarure, coordiantes and other information are computed for each node.

Nodes which share adjacent cells are connected. With the list of Nodes and their
connections is created a *graph*.

Parameters:

	`coors     :: Tuple{N{Vector{Number}}}` Tuple with the vector of the coordinates of the points

	`energy    :: Vector{Number}` Vector with the values of the contents or energy in each point

	`steps     :: Tuple{N{Number}}`` tuple the size of the pizel/voxels in each coordinate

	`threshold :: float = 0.0` minimum value of the contents to be create a cell, default = 0.

	`cellnode  :: Bool = false` create a node for each cell (true), otherwise create a node from
		the cells which are connected via the gradient

Returns:

	`DataCells`: DataType with information of the cells, include the coordinates,
				contents, the gradient, laplacian, etc...

	`DataNodes`: DataType with the information of the nodes, include coordinates,
				content, max and min laplacian, etc..

	`Spine`    : DataType with the graph that makes the spine of the cloud,
				the extremes nodes of the spine and the distance matrix between nodes

	`edges`    : Tuple with the ranges of the different coordinates

# Examples

```julia-repr
julia > coors  = ([1, 2, 3, 1, 2, 3, 1, 2, 3], [1, 1, 1, 2, 2, 2, 3, 3, 3])
julia > energy = [1, 2, 1, 2, 4, 2, 1, 2, 1]
julia > steps  = (1.0, 1.0))
julia > datacells, datanodes, spine, edges = clouds(coors, energy, steps)
julia > Tuple(datacells.grad)
(2.1213203435596424, 2.0, 2.1213203435596424, 2.0, 0.0, 2.0, 2.1213203435596424, 2.0, 2.1213203435596424)
julia > Tuple(datacells.lap)
(4.121320343559642, 0.0, 4.121320343559642, 0.0, -16.485281374238568, 0.0, 4.121320343559642, 0.0, 4.121320343559642)
```

"""
function clouds(coors     ::Tuple{N{VN}},  # Tuple with the point coordinates
	 	        energy    ::VN,            # energy or content of the point coordinates
				steps     ::TNN;           # Tuple with the step size in each coordinate
				threshold ::Float64 = 0.0, # energy threshold, default = 0.
				cellnode  ::Bool = false)  # if cellnode == true create a node for each cell
						                   # otherwise the node is composed by the cells whose gradient ends
										   # into a cell with null gradient! (that is, the higuest energy cell)

    ndim  = length(coors)
    nsize = length(coors[1])

	# assert inputs
	_assert(coors, energy, steps)

    # define the extended edges
    edges = (minimum(x) - 1.5*step : step : maximum(x) + 1.5*step for (x, step) in zip(coors, steps))
	edges = Tuple(Vector(ei) for ei in edges)

    # alias
	mm  = moves(ndim)
	m0  = mm.moves[mm.i0]
	ucoors = reduce(hcat, coors)

    # main histogram
    histo    = SB.fit(SB.Histogram, _hcoors(ucoors, m0), SB.weights(energy), edges)
    contents = deepcopy(histo.weights)
	cells    = findall(x -> x .> threshold, contents)
	coors_cells  = _xcoors(cells, edges)
	ucoors_cells = reduce(hcat, coors_cells)

    # deltas
    deltas = _deltas(ucoors, energy, edges, steps, mm, threshold)

    # gradient
    grad, igrad = _gradient(deltas, mm)

    # curvatures
    curves = _curvatures(deltas, mm)
    lap    = reduce(.+, deltas)

    # maximum and monimum curvatures
    curmax, icurmax, curmin, icurmin = _maxmin_curvatures(Base.size(lap), curves, mm)

	# nodes
    xnodes  = cellnode ?  Vector{Int64}(1:length(cells)) : _nodes(igrad, cells, mm)
    xborders, xneigh = _neighbour_node(ucoors_cells, xnodes, edges, steps, mm)

	xnodes_edges  = _nodes_edges(xnodes, xneigh)
	xcloud_nodes  = _cloud_nodes(xnodes_edges)
	xcloudid      = _cloudid(xnodes, xcloud_nodes)

	# cells data
	dfcells = DataCells(coors_cells, contents[cells], Tuple.(cells),
	 			grad[cells], igrad[cells], lap[cells],
                curmax[cells], icurmax[cells], curmin[cells], icurmin[cells],
			    xnodes, xcloudid, xborders[cells])

	# node data
	dfn = _dfnodes(dfcells, xnodes_edges, cellnode = cellnode)

	# create graph
	#graph, ecc = _graph(xnodes_edges)
	sp  = spine(xnodes_edges, dfn.contents)
	graph = sp.spine
	ecc   = sp.eccentricity
	extrs = zeros(Bool, length(ecc))
	for (i, j) in sp.extremes_maxcontents
		extrs[i] = true
		extrs[j] = true
	end
	#print(sp.extremes_maxcontents)
	dfn     = merge(dfn, (ecc = ecc, extremes = extrs))

	# create DataNodes
	dfnodes = DataNodes(dfn.size, dfn.contents, dfn.maxcontent,
			dfn.maxgrad, dfn.maxlap, dfn.minlap, dfn.curmax, dfn.curmin,
			dfn.nedges, dfn.cloud, dfn.coors, dfn.coors_std, dfn.coors_cell,
			dfn.ecc, dfn.extremes)

	return dfcells, dfnodes, sp, edges

end


#-----------------------------
# Nodes
#-----------------------------

function _dfnodes(xcl           ::DataCells,
	 			 nodes_edges    ::Dict{Int64, VI};
	 			 cellnode       ::Bool = false)
	# create extra information of the nodes

	nnodes   = maximum(xcl.node)
	nsize    = [sum(xcl.node .== inode) for inode in 1:nnodes]
	contents = [sum(xcl.contents[xcl.node .== inode]) for inode in 1:nnodes]
	maxcontent  = [maximum(xcl.contents[xcl.node .== inode]) for inode in 1:nnodes]
	maxgrad  = [maximum(xcl.grad[xcl.node .== inode])    for inode in 1:nnodes]
	maxlap   = [maximum(xcl.lap[xcl.node .== inode])     for inode in 1:nnodes]
	minlap   = [minimum(xcl.lap[xcl.node .== inode])     for inode in 1:nnodes]
	curmax   = [maximum(xcl.curmax[xcl.node .== inode])  for inode in 1:nnodes]
	curmin   = [minimum(xcl.curmin[xcl.node .== inode])  for inode in 1:nnodes]
	cloudid  = [maximum(xcl.cloud[xcl.node .== inode])   for inode in 1:nnodes]
	nedges   = [length(nodes_edges[inode]) for inode in 1:nnodes]

	function _stats(coor, inode)
		sel   = xcl.node .== inode
		norma = sum(xcl.contents[sel])
		mean  = sum(xcl.contents[sel] .* coor[sel])/norma
		std   = sqrt(sum(xcl.contents[sel] .* (coor[sel] .- mean).^2)/norma)
		return  mean, std
	end

	coors = Tuple([_stats(coor, inode)[1] for inode in 1:nnodes] for coor in xcl.coors)
	stds  = Tuple([_stats(coor, inode)[2] for inode in 1:nnodes] for coor in xcl.coors)

	_sel = inode -> cellnode ? xcl.node .== inode : (xcl.node .== inode) .&& (xcl.grad .== 0.0)
	coors_cell = Tuple([coor[_sel(inode)][1] for inode in 1:nnodes] for coor in xcl.coors)

	df = (size    = nsize  , contents = contents, maxcontent = maxcontent,
	      maxgrad = maxgrad, maxlap   = maxlap  , minlap     = minlap,
		  curmax  = curmax , curmin   = curmin,
		  nedges  = nedges,  cloud    = cloudid,
		  coors   = coors,   coors_std  = stds,
		  coors_cell = coors_cell)
	return df
end

#-----------------------------
# Graph
#-----------------------------

function _graph(nodes_edges ::Dict{Int64, VI})
	# create the graph of the Cloud
	nnodes = length(keys(nodes_edges))
	graphs = []
	g = GG.Graph(nnodes)
	for inode in keys(nodes_edges)
		for knode in nodes_edges[inode]
			GG.add_edge!(g, inode, knode)
		end
	end

	ecc = GG.eccentricity(g)
	ecc[ecc .> nnodes] .= 0
	return g, ecc

end


#-----------------------------
#  Clouds internal functions
#-----------------------------

function _xcoors(cells  ::Vector{<:CartesianIndex},
	 			 edges  ::Tuple{N{VF}})
	ndim    = length(edges)
	centers = [(ex[2:end] + ex[1:end-1])/2.0 for ex in edges]
	indices = [getindex.(cells, i) for i in 1:ndim]
	coors   = Tuple(ci[ii] for (ci, ii) in zip(centers, indices))
	return coors
end

function _hcoors(ucoors  ::T where T <: Array{<:Number},
			     move    ::VN)
	ndim = length(move)
	z    = ucoors .- move'
	zt   = Tuple(z[:, i] for i in 1:ndim)
	return zt
end

function _assert(coors   ::Tuple{N{VN}},
				 energy  ::VN,
				 steps   ::TNN)
	# asset inputs of clouds

	ndim  = length(coors)
	nsize = length(coors[1])

	# assert inputs
	if (ndim < 2) || (ndim > 3)
		throw(ArgumentError("clouds, only 2 or 3 coordinates allowed"))
	end

	for i in 2:ndim
		if (length(coors[i]) != nsize)
			throw(ArgumentError("clouds, length of all coordinates must be equal"))
		end
	end

	if (length(energy) != nsize)
		throw(ArgumentError("clouds, length of the energy must be equal to the coordinates"))
	end

	if (length(steps) != ndim)
		throw(ArgumentError("clouds, dimension of steps and coordinates must be equal"))
	end

	if (all([step > 0 for step in steps]) != true)
		throw(ArgumentError("clouds, steps in each coordinate must be greather than 0"))
	end

end

function _deltas(ucoors     ::T where T <: Array{<:Number},
	             energy     ::VN,
				 edges      ::Tuple{N{VN}},
				 steps      ::TNN,
				 m          ::Moves,
				 threshold  ::Float64)
	# Compute deltas in each of the moves

    his = [SB.fit(SB.Histogram, _hcoors(ucoors, steps .* move),
		SB.weights(energy), edges) for move in m.moves]
    contents = deepcopy(his[m.i0].weights)
	asize = Base.size(contents)
	mask = contents .<= threshold
	deltas   = [h.weights .- contents for h in his]
	for (delta, h) in zip(deltas, his)
		delta[h.weights .<= threshold] .= 0.0
		delta[mask] .= 0.0
	end

    dsteps   = [LA.norm(steps .* move) for move in m.moves]
    dsteps[m.i0] = 1.
    deltas  = [delta ./dstep for (delta, dstep) in zip(deltas, dsteps)]

    return deltas
end

function _gradient(deltas ::VAN,
	 			   m      ::Moves)
	# compute the gradiend from the deltas in each move
    dims   = Base.size(deltas[m.i0])
    d0     = deltas[m.i0]
    grad   = deepcopy(d0)
    igrad  = m.i0 .* ones(Int, dims...)
    for (i, di) in enumerate(deltas)
        imask         = di .> grad
        grad[imask]  .= di[imask]
        igrad[imask] .= i
    end
    return grad, igrad
end


function _curvatures(deltas ::VAN,
	                 m      ::Moves)
	# compute the curvatures in each direction from the deltas

    curvs = Dict{Int64, AN}()
    for smove in m.isym
       #curvs[smove[1]] = reduce(.+, [deltas[kmove] for kmove in m.iortho[smove]])
	   curvs[smove[1]] = reduce(.+, [deltas[kmove] for kmove in smove])
    end
    return curvs
end

function _maxmin_curvatures(nsize   ::TNI,
	 		                curves  ::Dict{Int64, AN},
							m       ::Moves)
	# compute the min, max, curvature from the curvatures grouped in a dictionary
    curmin  =  1e6 .* ones(nsize...)
    icurmin = m.i0 .* ones(Int, nsize...)
    curmax  = -1e6 .* ones(nsize...)
    icurmax = m.i0 .* ones(Int, nsize...)
    for smove in m.isym
		imove = smove[1]
        dd = curves[imove]
        mask1 = dd .> curmax
        curmax[mask1] .= dd[mask1]
        icurmax[mask1] .= imove
        mask2 = dd .< curmin
        curmin[mask2] .= dd[mask2]
        icurmin[mask2] .= imove
	end
    return curmax, icurmax, curmin, icurmin
end

function _node(cell   ::CartesianIndex,
	           igrad  ::AI,
			   m      ::Moves)
	# assing a node to a cell following the gradient
    imove = igrad[cell]
    if (imove == m.i0)
        return cell
    else
        #cindex_ = tuple(cindex) + m.moves[imove]
		nextcell = Tuple(Tuple(cell) .+ m.moves[imove])
        return _node(CartesianIndex(nextcell), igrad, m)
    end
end

function _nodes(igrad   ::AI,
	 	        cells   ::Vector{<:CartesianIndex},
				m       ::Moves)
	# compute the nodes to which the cells belong to

	cnodes  = [_node(cell, igrad, m) for cell in cells]
	ucnodes = unique(cnodes)
	dicnodes = Dict()
	for (i, cnode) in enumerate(ucnodes)
    	dicnodes[Tuple(cnode)] = i
	end
	nodes = [dicnodes[Tuple(cnode)] for cnode in cnodes]
	return nodes
end
#
# function _nodes_equal_values(cells, deltas, m)
# 	ncells    = length(cells)
# 	xnode     = zeros(Int64, ncells)
# 	currentid = 0
# 	for (i, cell) in enumerate(cells)
# 		if xnode[i] == 0
# 			currentid += 1
# 			xnode[i]   = currentid
# 		nodeid = xnode[i]
# 		mis    = [i for (i, delta) in enumerate(deltas) if delta[cell] == 0.0]
# 		if (length(mis) <= 1) continue
# 		for imove in mis
# 			if (imove == m.i0)
# 				continue
# 			end
# 			nextcell = Tuple(Tuple(cell) .+ m.moves[imove])
# 			jindices   = findall(x -> x .== nextcell, cells)
# 			@assert(length(jindices) == 1, "nodes equal values : only one cell should be accesible")
# 			jindex = jindices[1]
# 			if xnode[jindex] == 0
# 				xnode[jindex] = nodeid
# 			@asset(xnode[jindex] == nodeid), "nodes equal values : cells should have the same node id")
# 		end
# 	end
# 	return xnode
# end

function _neighbour_node(ucoors ::T where T <: Array{<:Number},
	                     nodes  ::VI,
						 edges  ::Tuple{N{VN}},
						 steps  ::TNN,
						 m      ::Moves)
	# compute the borders and neighbour cells

    his = [SB.fit(SB.Histogram, _hcoors(ucoors, steps .* move),
		SB.weights(nodes), edges) for move in m.moves]

    contents = deepcopy(his[m.i0].weights)
    mask     = contents .> 0
    borders  = [(h.weights .>0) .* (h.weights .!= contents) for h in his]
    nborders = reduce(.+, borders) .* mask

    return nborders, [h.weights for h in his]
end


function _nodes_edges(nodes  ::VI,
	 				  neighs ::VAN)
	# compute the dictionary that associates to a node the list of nodes
	# which are connected via adjacent cells

	#TODO Is this correct?
	imove0 = length(neighs) > 9 ? 14 : 5

	dus     = Dict{Int64, Array{Int64, 1}}()
	for inode in nodes
		imask = neighs[imove0] .== inode
		us = []
		for neigh in neighs
			kmask = imask .* (neigh .> 0) .* (neigh .!= inode)
			ius   = unique(vec(neigh[kmask]))
			for k in ius
				if !(k in us)
					append!(us, k)
				end
			end
		end
		dus[inode] = us
	end
	return dus
end

function _cloud_nodes(nodes_edges ::Dict{Int64, VI})
	# associate nodes to cloud indices
	# TODO: clean this ugly function!
	used   = []
	function _add(graph, inode)
		append!(graph, inode)
		append!(used , inode)
		knodes = nodes_edges[inode]
		for knode in knodes
			if !(knode in graph)
				_add(graph, knode)
			end
		end
	end
	graphs = Dict{Int64, VI}()
	i = 0
	for inode in keys(nodes_edges)
		if inode in used
			continue
		else
			igraph = [] # Array{Int64, 1}()
			_add(igraph, inode)
			i += 1
			graphs[i] = igraph
		end
	end
	return graphs
end

function _cloudid(xnodes        ::VI,
	              xclouds_nodes ::Dict{Int, VI})
	# associated to the nodes a cloud id

	graphid = zeros(Int64, length(xnodes))
	for i in keys(xclouds_nodes)
		for inode in xclouds_nodes[i]
			graphid[xnodes .== inode] .= i
		end
	end
	return graphid
end
