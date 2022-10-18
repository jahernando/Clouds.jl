using StatsBase
using LinearAlgebra
import Graphs as GG
using Clouds
using Test


function box_to_coors(aa)
    ndim = length(size(aa))
    steps = Tuple(ones(ndim))
    cells = CartesianIndices(aa)
    coors = Tuple(vec([c[i] for c in cells]) for i in 1:ndim)
    return coors, vec(aa), steps
end


function simple_clouds(bs; threshold = 0.0)

	# set the input for clouds
	ndim    = length(size(bs))
	indices = CartesianIndices(bs)
	steps   = ones(ndim)

	# prepare the local matrix to compute gradients, laplacian and curvatures
	nsize  = size(bs) .+ 2
	aa    = zeros(Float64, nsize...)
	for index in indices
		kindex = CartesianIndex(Tuple(index) .+ 1)
		aa[kindex] = bs[index]
	end

	# the mask
	mask = aa  .>  threshold
	nmask = aa .<= threshold
	# set the moves
	mm = moves(ndim)

	# compute the deltas
	function _delta(move)
		delta = zeros(Float64, nsize...)
		mod   = sqrt(sum(move .* move))
		mod   = mod == 0.0 ? 1 : mod
		for index in indices
			uindex = CartesianIndex(Tuple(index)  .+ 1)
			kindex = CartesianIndex(Tuple(Tuple(uindex) .+ move))
			dd = aa[kindex] > 0.0 ? (aa[kindex] - aa[uindex])/mod : 0.
			delta[uindex] = dd
		end
		return delta
	end

	# compute deltas
	deltas = [_delta(move) for move in mm.moves]

	# compute gradient
	ugrad = zeros(Float64, nsize...)
	igrad = mm.i0 .* ones(Int, nsize...)
	for (i, delta) in enumerate(deltas)
		imask = delta .> ugrad
		ugrad[imask] = delta[imask]
		igrad[imask] .= i
	end


	# compute and test laplacian
	lap = reduce(.+, deltas)

	# compute and test curves and max, min curve
	curves = [reduce(.+, [deltas[i] for i in s]) for s in mm.isym]
	#curves = [reduce(.+, [deltas[i] for i in mm.iortho[s]]) for s in mm.isym]
	curmax = -1e6*ones(nsize...)
	curmin = +1e6*ones(nsize...)
	icurmax = mm.i0 .* ones(Int, nsize...)
	icurmin = mm.i0 .* ones(Int, nsize...)
	for (i, curve) in enumerate(curves)
		imove = mm.isym[i][1]
		imask = curve .> curmax
		curmax[imask] .= curve[imask]
		icurmax[imask] .= imove
		imask = curve .< curmin
		curmin[imask] .= curve[imask]
		icurmin[imask] .= imove
	end
	icurmin[nmask] .= mm.i0
	icurmax[nmask] .= mm.i0

    data = (contents = aa[mask],
            grad = ugrad[mask], igrad = igrad[mask], lap = lap[mask],
            curmax = curmax[mask], icurmax = icurmax[mask],
            curmin = curmin[mask], icurmin = icurmin[mask])
    return data

end

function _stats(coors, contents)
	norma = sum(contents)
	means = [ sum(coor .* contents)/norma for coor in coors]
	stds  = [ sqrt(sum(contents .* (coor .- mean).^2)/norma)
			  for (mean, coor) in zip(means, coors)]
	return means, stds
end

#--------

# Test
#---------------------------

@testset "Clouds.jl" begin

    @testset "moves2d" begin
        movs = moves(2)
        ms = [[i, j] for i in -1:1:1 for j in -1:1:1]
        xdim = length(movs.moves[1])
        @test length(movs.moves) == 3^xdim
        @test findall(x -> all( x .== zeros(xdim)), movs.moves)[1] == movs.i0
        @test sum([all(v1 .== v2) for (v1, v2) in zip(movs.moves, ms)]) == length(ms)
        @test sum([movs.moves[i] == -movs.moves[j] for (i, j) in movs.isym]) == length(movs.isym)
        @test sum([all([dot(movs.moves[isym[1]], movs.moves[k]) == 0 for k in movs.iortho[isym]]) for isym in movs.isym]) == length(movs.isym)
    end

    @testset "moves3d" begin
        movs = moves(3)
        ms = [[i, j, k] for i in -1:1:1 for j in -1:1:1 for k in -1:1:1]
        xdim = length(movs.moves[1])
        @test length(movs.moves) == 3^xdim
        @test findall(x -> all( x .== zeros(xdim)), movs.moves)[1] == movs.i0
        @test sum([all(v1 .== v2) for (v1, v2) in zip(movs.moves, ms)]) == length(ms)
        @test sum([movs.moves[i] == -movs.moves[j] for (i, j) in movs.isym]) == length(movs.isym)
        @test sum([all([dot(movs.moves[isym[1]], movs.moves[k]) == 0 for k in movs.iortho[isym]]) for isym in movs.isym]) == length(movs.isym)
    end


    @testset "box2d" begin
        mat  = box2d()
        coors, contents, steps = box_to_coors(mat.contents)
        ndim = 2
        xcl, _, _, _  = clouds(coors, contents, steps)
        @test sum(contents) == sum(mat.contents)
        @test sum(isapprox.(xcl.contents, vec(mat.contents))) == 3^ndim
        @test sum(isapprox.(xcl.grad    , vec(mat.grad)))     == 3^ndim
        @test sum(isapprox.(xcl.lap     , vec(mat.lap)))      == 3^ndim
        @test sum(isapprox.(xcl.curmax  , vec(mat.curmax)))   == 3^ndim
        @test sum(isapprox.(xcl.curmin  , vec(mat.curmin)))   == 3^ndim
    end

    @testset "box3d" begin
        mat  = box3d()
        coors, contents, steps = box_to_coors(mat.contents)
        ndim = 3
        xcl, _, _, _  = clouds(coors, contents, steps)
        @test sum(contents) == sum(mat.contents)
        @test sum(isapprox.(xcl.contents, vec(mat.contents))) == 3^ndim
        @test sum(isapprox.(xcl.grad    , vec(mat.grad)))     == 3^ndim
        @test sum(isapprox.(xcl.lap     , vec(mat.lap)))      == 3^ndim
        @test sum(isapprox.(xcl.curmax  , vec(mat.curmax)))   == 3^ndim
        @test sum(isapprox.(xcl.curmin  , vec(mat.curmin)))   == 3^ndim
    end

    @testset "simple" begin
        mat  = box3d()
        ndim = 3
        scl  = simple_clouds(mat.contents)
        coors, contents, steps = box_to_coors(mat.contents)
        xcl, _, _, _  = clouds(coors, contents, steps)
        @test sum(isapprox.(xcl.contents, vec(scl.contents))) == 3^ndim
        @test sum(isapprox.(xcl.grad    , vec(scl.grad)))     == 3^ndim
        @test sum(isapprox.(xcl.lap     , vec(scl.lap)))      == 3^ndim
        @test sum(isapprox.(xcl.curmax  , vec(scl.curmax)))   == 3^ndim
        @test sum(isapprox.(xcl.curmin  , vec(scl.curmin)))   == 3^ndim
        @test sum(scl.igrad   .== scl.igrad)                  == 3^ndim
        @test sum(scl.icurmax .== scl.icurmax)                == 3^ndim
        @test sum(scl.icurmin .== scl.icurmin)                == 3^ndim
    end

	@testset "nodes" begin
		nn   = sample(1:10)
		b2   = box2d()
		b2n  = repeat(b2.contents, nn)
		coors, contents, steps = box_to_coors(b2n)
		xcl, xnd, xgraph, _ = clouds(coors, contents, steps)
		@test maximum(xcl.node)    == nn
		@test length(xnd.contents) == nn
		@test GG.nv(xgraph)        == nn
		@test maximum(xcl.cloud)   == 1
		xcl, xnd, xgraph, _ = clouds(coors, contents, steps, cellnode = true)
		@test maximum(xcl.node)    == length(xcl.node)
		@test length(xnd.contents) == length(xcl.node)
		@test GG.nv(xgraph)        == length(xcl.node)
		@test maximum(xcl.cloud)   == 1
		b3  = box3d()
		b3n = repeat(b3.contents, nn)
		coors, contents, steps = box_to_coors(b3n)
		xcl, xnd, xgraph, _ = clouds(coors, contents, steps)
		@test maximum(xcl.node)    == nn
		@test length(xnd.contents) == nn
		@test GG.nv(xgraph)        == nn
		@test maximum(xcl.cloud)   == 1
		xcl, xnd, xgraph, _ = clouds(coors, contents, steps, cellnode = true)
		@test maximum(xcl.node) == length(xcl.node)
		@test length(xnd.contents) == length(xcl.node)
		@test GG.nv(xgraph)        == length(xcl.node)
		@test maximum(xcl.cloud)   == 1
	end

	@testset "cloudid" begin
		nn   = sample(1:10)
		b2   = box2d()
		m2   = b2.contents
		b2n  = repeat(vcat(m2, 0 .* m2), nn)
		coors, contents, steps = box_to_coors(b2n)
		xcl, xnd, _, _ = clouds(coors, contents, steps)
		@test maximum(xcl.node)  == nn
		@test maximum(xcl.cloud) == nn
		b3  = box3d()
		m3  = b3.contents
		b3n = repeat(vcat(m3, 0 .* m3), nn)
		coors, contents, steps = box_to_coors(b3n)
		xcl, xnd, xgraph, _ = clouds(coors, contents, steps)
		@test maximum(xcl.node)  == nn
		@test maximum(xcl.cloud) == nn
	end

	@testset "nodedata" begin
		b2   = box2d()
		m2   = b2.contents
		coors, contents, steps = box_to_coors(b2.contents)
		_, xnd, _, _ = clouds(coors, contents, steps)
		@test xnd.size[1]       == length(contents)
		@test xnd.contents[1]   == sum(contents)
		@test xnd.maxcontent[1] == maximum(contents)
		@test xnd.maxgrad[1]    == maximum(b2.grad)
		@test xnd.maxlap[1]     == maximum(b2.lap)
		@test xnd.minlap[1]     == minimum(b2.lap)
		@test xnd.curmax[1]     == maximum(b2.curmax)
		@test xnd.curmin[1]     == minimum(b2.curmin)
		@test xnd.nedges[1]     == 0
	end

	@testset "nodescoors" begin
		b2 = box2d()
		coors, contents, steps = box_to_coors(b2.contents)
		xcl, xnd, graph, edges = clouds(coors, contents, steps)
		i0  = argmax(contents)
		means, stds = _stats(coors, contents)
		@test all(xnd.coors_cell .== Tuple([c[i0],] for c in coors))
		@test all([c[1] .== mean for (c, mean) in zip(xnd.coors, means)])
		@test all([c[1] .== std  for (c, std)  in zip(xnd.coors_std, stds)])
		b3 = box3d()
		coors, contents, steps = box_to_coors(b3.contents)
		xcl, xnd, graph, edges = clouds(coors, contents, steps)
		i0  = argmax(contents)
		means, stds = _stats(coors, contents)
		@test all(xnd.coors_cell .== Tuple([c[i0],] for c in coors))
		@test all([c[1] .== mean for (c, mean) in zip(xnd.coors, means)])
		@test all([c[1] .== std  for (c, std)  in zip(xnd.coors_std, stds)])
		b2 = box2d()
		coors, contents, steps = box_to_coors(b2.contents)
		xcl, xnd, graph, edges = clouds(coors, contents, steps, cellnode = true)
		@test all(xnd.coors_cell .== xnd.coors) == true
		@test sum(sum(xnd.coors_std)) == 0.0
	end

end
