#--------------------
#
# Main methods to discrtize, compute the gradient, laplacian, etc
#
#--------------------

import StatsBase     as SB
import LinearAlgebra as LA

export discretize, discrete_gradient

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


#----------------------------
# Main functions
#----------------------------

function discretize(coors     ::Tuple{N{VN}},  # Tuple with the point coordinates
                    contents  ::VN,            # energy or content of the point coordinates
                    steps     ::TNN;           # Tuple with the step size in each coordinate
                    threshold ::Float64 = 0.0) # energy threshold, default = 0.

    ndim  = length(coors)
              
    # assert inputs
    _assert(coors, contents, steps)
              
    # define the extended edges
    edges = (minimum(x) - 1.5*step : step : maximum(x) + 1.5*step for (x, step) in zip(coors, steps))
    edges = Tuple(Vector(ei) for ei in edges)
              
    # alias
    mm  = moves(ndim)
    m0  = mm.moves[mm.i0]
    ucoors = reduce(hcat, coors)
              
    # main histogram
    histo    = SB.fit(SB.Histogram, _hcoors(ucoors, m0), SB.weights(contents), edges)
    xcontents = deepcopy(histo.weights)
    cells    = findall(x -> x .> threshold, xcontents)
    coors_cells  = _xcoors(cells, edges)
    #ucoors_cells = reduce(hcat, coors_cells)
              
    return (coors = coors_cells, contents = xcontents, cells = cells, edges = edges)          

end

function discrete_gradient(coors, contents, steps; threshold = 0.0, full_output = false)

    ndim  = length(coors)
    mm    = moves(ndim)

    #digitalize
    xcoors, xcontents, cells, edges = discretize(coors, contents, steps; threshold = threshold)

    # deltas
    ucoors = reduce(hcat, coors)
    deltas = _deltas(ucoors, contents, edges, steps, mm, threshold)

    # gradient
    grad, igrad = _gradient(deltas, mm)

    # nodes
    xnodes = _nodes(cells, igrad)

    if (full_output == false)
        res = (coors = xcoors, contents = xcontents, cells = cells,
               grad = grad[cells], graddir = igrad[cells], node = xnodes)
        return res, edges
    end

    # curvatures
    curves = _curvatures(deltas, mm)
    lap    = reduce(.+, deltas)
    
    # maximum and monimum curvatures
    curmax, icurmax, curmin, icurmin = _maxmin_curvatures(Base.size(lap), curves, mm)
    
    res = (coors = xcoors, contents = xcontents, cells = cells,
           grad   = grad[cells], graddir = igrad[cells], node = xnodes,
           lap = lap[cells],
           curmax = curmax[cells], curmaxdir = icurmax[cells],
           curmin = curmin[cells], curmindir = icurmin[cells])
    return res, edges
end

#-----------------------------
#  internal functions
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
    his = [SB.fit(SB.Histogram, _hcoors(ucoors, steps .* move), SB.weights(energy), edges) for move in m.moves]
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


function _nodes(cells   ::Vector{<:CartesianIndex},
                igrad   ::AI)
    # compute the nodes to which the cells belong to
    ndim = length(cells[1])
    m    = moves(ndim)

    cnodes  = [_node(cell, igrad, m) for cell in cells]
    ucnodes = unique(cnodes)
    dicnodes = Dict()
    for (i, cnode) in enumerate(ucnodes)
        dicnodes[Tuple(cnode)] = i
    end
    nodes = [dicnodes[Tuple(cnode)] for cnode in cnodes]
    return nodes
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


# #
# # function _nodes_equal_values(cells, deltas, m)
# # 	ncells    = length(cells)
# # 	xnode     = zeros(Int64, ncells)
# # 	currentid = 0
# # 	for (i, cell) in enumerate(cells)
# # 		if xnode[i] == 0
# # 			currentid += 1
# # 			xnode[i]   = currentid
# # 		nodeid = xnode[i]
# # 		mis    = [i for (i, delta) in enumerate(deltas) if delta[cell] == 0.0]
# # 		if (length(mis) <= 1) continue
# # 		for imove in mis
# # 			if (imove == m.i0)
# # 				continue
# # 			end
# # 			nextcell = Tuple(Tuple(cell) .+ m.moves[imove])
# # 			jindices   = findall(x -> x .== nextcell, cells)
# # 			@assert(length(jindices) == 1, "nodes equal values : only one cell should be accesible")
# # 			jindex = jindices[1]
# # 			if xnode[jindex] == 0
# # 				xnode[jindex] = nodeid
# # 			@asset(xnode[jindex] == nodeid), "nodes equal values : cells should have the same node id")
# # 		end
# # 	end
# # 	return xnode
# # end

# function _neighbour_node(ucoors ::T where T <: Array{<:Number},
#            nodes  ::VI,
#            edges  ::Tuple{N{VN}},
#            steps  ::TNN,
#            m      ::Moves)
# # compute the borders and neighbour cells

# his = [SB.fit(SB.Histogram, _hcoors(ucoors, steps .* move),
# SB.weights(nodes), edges) for move in m.moves]

# contents = deepcopy(his[m.i0].weights)
# mask     = contents .> 0
# borders  = [(h.weights .>0) .* (h.weights .!= contents) for h in his]
# nborders = reduce(.+, borders) .* mask

# return nborders, [h.weights for h in his]
# end


# function _nodes_edges(nodes  ::VI,
#          neighs ::VAN)
# # compute the dictionary that associates to a node the list of nodes
# # which are connected via adjacent cells

# #TODO Is this correct?
# imove0 = length(neighs) > 9 ? 14 : 5

# dus     = Dict{Int64, Array{Int64, 1}}()
# for inode in nodes
# imask = neighs[imove0] .== inode
# us = []
# for neigh in neighs
# kmask = imask .* (neigh .> 0) .* (neigh .!= inode)
# ius   = unique(vec(neigh[kmask]))
# for k in ius
#   if !(k in us)
#       append!(us, k)
#   end
# end
# end
# dus[inode] = us
# end
# return dus
# end

# function _cloud_nodes(nodes_edges ::Dict{Int64, VI})
# # associate nodes to cloud indices
# # TODO: clean this ugly function!
# used   = []
# function _add(graph, inode)
# append!(graph, inode)
# append!(used , inode)
# knodes = nodes_edges[inode]
# for knode in knodes
# if !(knode in graph)
#   _add(graph, knode)
# end
# end
# end
# graphs = Dict{Int64, VI}()
# i = 0
# for inode in keys(nodes_edges)
# if inode in used
# continue
# else
# igraph = [] # Array{Int64, 1}()
# _add(igraph, inode)
# i += 1
# graphs[i] = igraph
# end
# end
# return graphs
# end

# function _cloudid(xnodes        ::VI,
#     xclouds_nodes ::Dict{Int, VI})
# # associated to the nodes a cloud id

# graphid = zeros(Int64, length(xnodes))
# for i in keys(xclouds_nodes)
# for inode in xclouds_nodes[i]
# graphid[xnodes .== inode] .= i
# end
# end
# return graphid
# end
