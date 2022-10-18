
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
VI  = Vector{Int64}
VF  = Vector{Float64}
VN  = Vector{<:Number}


# Helper Data Type to hold one step movements in a 2D or 3D grid
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
end



#----------
# Moves
#----------

"""
Function that returns a Moves Data Type. Moves holds the one step moves in a 2D or 3D grid.
It serves as a help data type for Clouds

Parameters:
	ndim: int, 2 or 3
Returns:
	Moves: DataType with the movemeents

"""
function moves(ndim::Int64)
	moves = ndim == 2 ? Tuple([i, j] for i in -1:1:1 for j in -1:1:1) :
	     Tuple([i, j, k] for i in -1:1:1 for j in -1:1:1 for k in -1:1:1)
	move0 = ndim == 2 ? [0, 0] : [0, 0, 0]
	kmove0 = [i for (i, move) in enumerate(moves) if (move == move0)][1]
	smoves = Tuple((i, j) for (i, movei) in enumerate(moves)
	          for (j, movej) in enumerate(moves)
		      if (movei == -1 .* movej ) &  (i > j))
	omoves = Dict{T2I, VI}()
	for ii in smoves
		movei = moves[ii[1]]
		omoves[ii] = [j for (j, movej) in enumerate(moves) if ((sum(movej .* movei) == 0.0) & (sum(movej .* movej) != 0.0))]
	end
	#return (moves = moves, i0 = kmove0, isym = smoves, iortho = omoves)
	return Moves(moves, kmove0, smoves, omoves)
end


#-----
# Clouds
#-----

"""

Clouds Main Function

From a collection of points (coors) with energy. The Clouds algorithm creates a
pixel or voxelixed space where it compute the local gradient, laplacian of each
pixel or voxel (called cells).
Nodes based on each cell or the cell connected by the gradient
are used to define a graph. Extra information of the nodes is also provided.

Parameters:
	coors : Tuple with the vector of the coordinates of the points
	energy: Vector with the values of the contents or energy in each point
	steps : Tuple with the size of the pizel/voxels in each coordinate
	threshold: minimum value of the contents to be create a cell, default = 0.
	cellnode : create a node for each cell (true), otherwise create a node from
		the cells which are connected via the gradient

Returns:
	DataCells: DataType with information of the cells, include the coordinates,
				contents, the gradient, laplacian, etc...
	DataNodes: DataType with the information of the nodes, include coordinates,
				content, max and min laplacian, etc..
	Graph    : GraphType from Graphs containing the graph of the cloud
	edges    : Tuple with the ranges in the different coordinates

"""
function clouds(coors     ::Tuple{N{VN}},  # Tuple with the point coordinates
	 	        energy    ::VN,            # energy or content of the point coordinates
				steps     ::TNV;           # Tuple with the step size in each coordinate
				threshold ::Float64 = 0.0, # energy threshold, default = 0.
				cellnode  ::Bool = false)  # if cellnode == true create a node for each cell
						                   # otherwise the node is composed by the cells whose gradient ends
										   # into a cell with null gradient! (that is, the higuest energy cell)

    ndim  = length(coors)
    nsize = length(coors[1])

	if (ndim < 2) || (ndim > 3)
		throw(ArgumentError("only 2 or 3 coordinates allowed"))
	end

    # assert dimensions
    for i in 2:ndim
		if (length(coors[i]) != nsize)
			throw(ArgumentError("length of all coordinates must be equal"))
		end
    end

    if (length(energy) != nsize)
		throw(ArgumentError("length of the energy must be equal to the coordinates"))
	end

	if (length(steps) != ndim)
		throw(ArgumentError("dimension of steps and coordinates must be equal"))
	end

    # define the extended edges
    edges = Tuple(minimum(x) - 1.5*step : step : maximum(x) + 1.5*step for (x, step) in zip(coors, steps))

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
    xnodes  = cellnode ?  (1:length(cells)) : _nodes(igrad, cells, mm)
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
	graph, ecc = _graph(xnodes_edges)
	dfn     = merge(dfn, (ecc = ecc,))

	# create DataNodes
	dfnodes = DataNodes(dfn.size, dfn.contents, dfn.maxcontent,
			dfn.maxgrad, dfn.maxlap, dfn.minlap, dfn.curmax, dfn.curmin,
			dfn.nedges, dfn.cloud, dfn.coors, dfn.coors_std, dfn.coors_cell,
			dfn.ecc)

	return dfcells, dfnodes, graph, edges

end


#-----------------------------
# Nodes
#-----------------------------

function _dfnodes(xcl, nodes_edges; cellnode = false)
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

	#weights  = _w(xcl)

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

function _graph(nodes_edges)
	# grates the graph
	nnodes = length(keys(nodes_edges))
	graphs = []
	g = GG.Graph(nnodes)
	for inode in keys(nodes_edges)
		for knode in nodes_edges[inode]
			GG.add_edge!(g, inode, knode)
		end
	end

	ecc = GG.eccentricity(g)
	return g, ecc

end


#-----------------------------
#  Clouds internal functions
#-----------------------------


function _xcoors(cells, edges)
	ndim    = length(edges)
	centers = [(ex[2:end] + ex[1:end-1])/2.0 for ex in edges]
	indices = [getindex.(cells, i) for i in 1:ndim]
	coors   = Tuple(ci[ii] for (ci, ii) in zip(centers, indices))
	return coors
end

function _hcoors(ucoors, move)
	ndim = length(move)
	z    = ucoors .- move'
	zt   = Tuple(z[:, i] for i in 1:ndim)
	return zt
end


function _deltas(ucoors, energy, edges, steps, m, threshold)

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

function _gradient(deltas, m)
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


function _curvatures(deltas, m)

    curvs = Dict()
    for smove in m.isym
       #curvs[smove[1]] = reduce(.+, [deltas[kmove] for kmove in m.iortho[smove]])
	   curvs[smove[1]] = reduce(.+, [deltas[kmove] for kmove in smove])
    end
    return curvs
end

function _maxmin_curvatures(nsize, curves, m)
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

function _node(cell, igrad, m)
    imove = igrad[cell]
    if (imove == m.i0)
        return cell
    else
        #cindex_ = tuple(cindex) + m.moves[imove]
		nextcell = Tuple(Tuple(cell) .+ m.moves[imove])
        return _node(CartesianIndex(nextcell), igrad, m)
    end
end

function _nodes(igrad, cells, m)

	cnodes  = [_node(cell, igrad, m) for cell in cells]
	ucnodes = unique(cnodes)
	dicnodes = Dict()
	for (i, cnode) in enumerate(ucnodes)
    	dicnodes[Tuple(cnode)] = i
	end
	nodes = [dicnodes[Tuple(cnode)] for cnode in cnodes]
	return nodes
end


function _neighbour_node(ucoors, nodes, edges, steps, m)

    his = [SB.fit(SB.Histogram, _hcoors(ucoors, steps .* move),
		SB.weights(nodes), edges) for move in m.moves]

    contents = deepcopy(his[m.i0].weights)
    mask     = contents .> 0
    borders  = [(h.weights .>0) .* (h.weights .!= contents) for h in his]
    nborders = reduce(.+, borders) .* mask

    return nborders, [h.weights for h in his]
end


function _nodes_edges(nodes, neighs)

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

function _cloud_nodes(nodes_edges)
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
	graphs = Dict()
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

function _cloudid(xnodes, xclouds_nodes)

	graphid = zeros(Int64, length(xnodes))
	for i in keys(xclouds_nodes)
		for inode in xclouds_nodes[i]
			graphid[xnodes .== inode] .= i
		end
	end
	return graphid
end
