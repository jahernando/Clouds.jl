
#using Base
import StatsBase     as SB
import LinearAlgebra as LA
import Graphs        as GG
#import DataFrames as DF


#export mesh, quiver
export moves, clouds, nodes

#----------
# Moves
#----------

"""
Function that returns a NamedTuple with the one step movement vectors
in 2D and 3D.

It also returns the null move, the symetric pair of movees,
and the ortogonal move associated to any symmetric movement.
The null (*i0*), symmetric (*isym*), and ortogonal movs (*iorto*)
are indexed respect the main vector of moves (*moves*)

i.e 2D: null [0., 0.], symmetric pair ([1, 0], [-1, 0]),
ortogonal directions of the previous pair ([0, 1], [0, -1])

"""
function moves(ndim)
	moves = ndim == 2 ? [[i, j] for i in -1:1:1 for j in -1:1:1] : [[i, j, k] for i in -1:1:1 for j in -1:1:1 for k in -1:1:1]
	move0 = ndim == 2 ? [0, 0] : [0, 0, 0]
	kmove0 = [i for (i, move) in enumerate(moves) if (move == move0)][1]
	smoves = [(i, j) for (i, movei) in enumerate(moves) for (j, movej) in enumerate(moves) if (movei == -1 .* movej ) &  (i > j)]
	omoves = Dict()
	for ii in smoves
		movei = moves[ii[1]]
		omoves[ii] = [j for (j, movej) in enumerate(moves) if ((sum(movej .* movei) == 0.0) & (sum(movej .* movej) != 0.0))]
	end
	return (moves = moves, i0 = kmove0, isym = smoves, iortho = omoves)
end


#-----
# Clouds
#-----

"""

Create Clouds

From a (x, y) or (x, y, z) tuples of points associated to contents.

The space is binned with *step* dimensiones and a histogram is created with the contents.

Only bins with contents > threshold are considered valid cells.

For each cell several information is provided:

contents, gradient, laplacian, minimum and maximum curvature, node number,
and number of neighbour cells belonging to another node, and graph id

Other information:

  * A list, graphs, graph_nodes, with the list of nodes in each graph
  * A dictionary, nodes_edges, with the edes between the nodes.

"""
function clouds(coors, energy, steps; threshold = 0., cellnode = false)

    ndim  = length(coors)
    nsize = length(coors[1])

    # assert dimensions
    for i in 2:ndim
        @assert(length(coors[i]) == nsize)
    end
    @assert(length(energy) == nsize)
    @assert(length(steps)  == ndim)

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

    # deltas
    deltas = _deltas(ucoors, energy, edges, steps, mm, threshold)

    # gradient
    grad, igrad = _gradient(deltas, mm)

    # curvatures
    curves = _curvatures(deltas, mm)
    lap    = reduce(.+, deltas)

    # maximum and monimum curvatures
    curmax, icurmax, curmin, icurmin = _maxmin_curvatures(Base.size(lap), curves, mm)

    #xs, ys = mesh(edges)

    # nodes
    xnodes  = cellnode ?  (1:length(cells)) : _nodes(igrad, cells, mm)
    xborders, xneigh = _neighbour_node(ucoors, xnodes, edges, steps, mm)

	xlinks  = _links(xnodes, xneigh)

	xgraphs  = _graph_nodes(xlinks)
	xgraphid = _graphid(xnodes, xgraphs)

    # output
    xcl = (edges = edges,  coors = coors, contents = contents[cells],
			cells = cells,
            grad = grad[cells], igrad = igrad[cells],
            lap = lap[cells], #curves = curves,
            curmax = curmax[cells], icurmax = icurmax[cells],
            curmin = curmin[cells], icurmin = icurmin[cells],
            node   = xnodes, graph = xgraphid,
            nborders = xborders[cells],
			graphs_nodes = xgraphs,
			nodes_edges = xlinks)
	return xcl

end


#-----------------------------
# Nodes table
#-----------------------------

function nodes(xcl)
	nnodes   = maximum(xcl.node)
	nsize    = [sum(xcl.node .== inode) for inode in 1:nnodes]
	contents = [sum(xcl.contents[xcl.node .== inode]) for inode in 1:nnodes]
	maxcontent  = [maximum(xcl.contents[xcl.node .== inode]) for inode in 1:nnodes]
	maxgrad  = [maximum(xcl.grad[xcl.node .== inode]) for inode in 1:nnodes]
	maxlap   = [maximum(xcl.lap[xcl.node .== inode]) for inode in 1:nnodes]
	minlap   = [minimum(xcl.lap[xcl.node .== inode]) for inode in 1:nnodes]
	maxcur   = [maximum(xcl.curmax[xcl.node .== inode]) for inode in 1:nnodes]
	mincur   = [minimum(xcl.curmin[xcl.node .== inode]) for inode in 1:nnodes]
	nlinks   = [length(xcl.nodes_edges[inode]) for inode in 1:nnodes]
	data = (size = nsize, contents = contents, maxcontent = maxcontent,
	        maxgrad = maxgrad, maxlap = maxlap, minlap = minlap,
			maxcur = maxcur, mincur = mincur,
			nedges = nlinks)
	return data
end

#---------------
#  Graphs
#----------------


function create_graphs(nodes_edges)

	nnodes = length(keys(nodes_edges))
	graphs = []
	g = GG.Graph(nnodes)
	for inode in keys(nodes_edges)
		for knode in nodes_edges[inode]
			GG.add_edge!(g, inode, knode)
		end
	end

	ecc = GG.eccentricity(g)
	_mst = GG.prim_mst(g)
	mst = [mi.src for mi in _mst]
	append!(mst, _mst[length(_mst)].dst)

	return g, ecc, _mst

end


#-----------------------------
#  Clouds internal functions
#-----------------------------

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


function _links(nodes, neighs)

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


function _graph_nodes(nodes_edges)
	function _add(graph, inode)
		if (inode in graph)
			return graph
		end
		append!(graph, inode)
		knodes = nodes_edges[inode]
		for knode in knodes
			graph = _add(graph, knode)
		end
		return graph
	end
	graphs = []
	for inode in keys(nodes_edges)
		if any([inode in graph for graph in graphs])
			continue
		else
			igraph = Array{Int64, 1}()
			igraph = _add(igraph, inode)
			append!(graphs, igraph)
		end
	end
	return graphs
end

function _graphid(xnodes, xgraphs)

	graphid = zeros(Int64, length(xnodes))
	for (i, graph) in enumerate(xgraphs)
		for inode in graph
			graphid[xnodes .== inode] .= i
		end
	end
	return graphid
end
