#import StatsBase as SB
#import LinearAlgebra as LA
import Graphs as GG

export Spine, spine

N   = Vararg
T2I = Tuple{Int64, Int64}
TNI = Tuple{N{Int64}}
TNV = Tuple{N{Float64}}
TNN = Tuple{N{<:Number}}
VI  = Vector{Int64}
VF  = Vector{Float64}
VN  = Vector{<:Number}
VAN = Vector{T} where T <: Array{<:Number}
AN  = T where T <: Array{<:Number}
AI  = T where T <: Array{Int64}

struct Spine

	spine        #:: Graph
	eccentricity #:: VI
	dists        #:: Matrix
	extremes     #:: Tuple of Int Pairs
	extremes_maxcontents #::Tuple of Int Pairs
end


function spine(nodes_edges :: Dict{Int64, VI},
		       contents    :: VN)
	# create the graph of the Cloud
	nnodes = length(keys(nodes_edges))
	g = GG.Graph(nnodes)
	for inode in keys(nodes_edges)
		for knode in nodes_edges[inode]
			GG.add_edge!(g, inode, knode)
		end
	end

	ecc = GG.eccentricity(g)
	ecc[ecc .> nnodes] .= 0

	extremes, dist = _extremes(g)

	extremes_maxcontents = _extremes_maxcontents(extremes, contents)

	sp = Spine(g, ecc, dist, extremes, extremes_maxcontents)

	return sp

end

#-----
# Internal functions

function _extremes(graph)
	nvertices = GG.nv(graph)
	dists     = zeros(Int64, nvertices, nvertices)
	for i in 1:nvertices
		ds = GG.dijkstra_shortest_paths(graph, i)
		vv = copy(ds.dists)
		vv[vv .> nvertices] .= 0
		dists[i, :] = vv
	end
	extremes = findall(x -> x .>= maximum(dists), dists)
	extremes = [(i, j) for (i, j) in Tuple.(extremes) if i < j]
	return extremes, dists
end

function _extremes_maxcontents(extremes, contents)
	contents = [contents[i] + contents[j] for (i, j) in extremes]
	ids      = findall(x -> x .== maximum(contents), contents)
	return extremes[ids]
end
