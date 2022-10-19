import StatsBase as SB
import Images.ImageFiltering as IF

export box2d, box3d, line, box_to_coors

# sort-cuts type definitions
N   = Vararg
T2I = Tuple{Int64, Int64}
TNI = Tuple{N{Int64}}
TNV = Tuple{N{Float64}}
TNN = Tuple{N{<:Number}}
VI  = Vector{Int64}
VF  = Vector{Float64}
VN  = Vector{<:Number}
AN  = T where T <: Array{<:Number}
AI  = Array{<:Int64}

"""
returns a (3, 3) 2D matrix:

  A = [a b a; b c b; a b c]

requires  a < b < c
"""
function box2d(vals::VN = [1, 2, 4])

	@assert length(vals) ==3

	a, b, c = vals[1], vals[2], vals[3]
	@assert (a < b) & (b < c)

	cs = [a b a; b c b; a b a]

	mat = cs

	ba = b-a
	cb = c-b
	ca = (c-a)/sqrt(2.)

	c1, c2, c3 = max(ba, ca), cb, 0.
	grad = [c1 c2 c1; c2 c3 c2; c1 c2 c1]

	c1 = 2*ba + ca
	c2 = -2*ba + cb
	c3 = -4*cb - 4*ca
	lap = [c1 c2 c1; c2 c3 c2; c1 c2 c1]

	c1, c2, c3 = ca, cb, 2*max(-ca, -cb)
	maxc = [c1 c2 c1; c2 c3 c2; c1 c2 c1]

	# notice that a diagonal step from a corner out of the cloud can have curvature 0.
	c1, c2, c3 = min(ba, ca, 0.0), -2*ba, min(-2*cb, -2*ca)
	minc = [c1  c2 c1; c2 c3 c2; c1 c2 c1]

	return (contents = mat, grad = grad, lap = lap, curmax = maxc,
    curmin = minc)

end

"""
return a (3, 3, 3) 3D array A:
 	A[:, :, 1] = [a b a; b c b; a b c]
	A[:, :, 2] = [b c b; c d c; b c b]
	A[:, :, 3] = [a b a; b c b; a b a]

requires a < b < c < d
"""
function box3d(vals::VN = [1, 2, 3, 4])

	@assert (length(vals) == 4)

	a, b, c, d = vals[1], vals[2], vals[3], vals[4]
	@assert (a < b) & (b < c) & (c < d)

	bs = zeros(3, 3, 3)
	bs[:, :, 1] = [a b a; b c b; a b a]
	bs[:, :, 2] = [b c b; c d c; b c b]
	bs[:, :, 3] = [a b a; b c b; a b a]

	ba = b-a
	ca = c-a
	da = d-a
	cb = c-b
	db = d-b
	dc = d-c
	d2 = sqrt(2.)
	d3 = sqrt(3.)

	function _mat(c1, c2, c3, c4)
		m = zeros(3, 3, 3)
		m[:, :, 1] = [c1 c2 c1; c2 c3 c2; c1 c2 c1]
		m[:, :, 2] = [c2 c3 c2; c3 c4 c3; c2 c3 c2]
		m[:, :, 3] = [c1 c2 c1; c2 c3 c2; c1 c2 c1]
		return m
	end

	c1, c2, c3, c4 = max(ba, ca/d2, da/d3), max(cb, db/d2), dc, 0.
	grad = _mat(c1, c2, c3, c4)

	c1     =  3*ba + 3*ca/d2 + da/d3
	c2     = -2*ba + 2*cb + 2*cb/d3 + db/d2
	c3     = -4*ca/d2 - 4*cb - 4*cb/d3 + dc
	c4     = -8*da/d3 - 12*db/d2 - 6*dc
	lap    = _mat(c1, c2, c3, c4)

	c1     = max(ba, ca/d2, da/d3)
	c2     = max(cb, db/d2)
	c3     = dc
	c4     = 2*max(-da/d3, -db/d2, -dc)
	curmax = _mat(c1, c2, c3, c4)

	c1     = 0 # due to the effect of the diagonal movement outside the cloud!
	c2     = -2ba
	c3     = min(-2cb, -2ca/d2, -dc)
	c4     = min(-2dc, -2db/d2, -2da/d3)
	curmin = _mat(c1, c2, c3, c4)

	return (contents = bs, grad = grad, lap = lap,
		    curmax = curmax, curmin = curmin)


end


function box_to_coors(aa::AN)
    ndim = length(size(aa))
    steps = Tuple(ones(ndim))
    cells = CartesianIndices(aa)
    coors = Tuple(vec([c[i] for c in cells]) for i in 1:ndim)
    return coors, vec(aa), steps
end

"""
returns 3D or 2D points in a line smeared with a gaussian

"""

"""
Just a smeared line!
"""
function line(;ndim::Int64 = 3, threshold::Float64 = 0.)

	tstep = 0.1
	ts = 0:tstep:1.
	ax, bx, cx = 5., 0., 0.
	ay, by, cy = -5., -5., 0.
	az, bz, cz = 5., -5., 0.
	xx = cx .* ts .* ts + ax .* ts .+ bx
	yy = cy .* ts .* ts + ay .* ts .+ by
	zz = cz .* ts .* ts + az .* ts .+ bz

	zsig  = 5.
	sigma = 2 * tstep #* ma(ax, ay, az)
	xxbins = minimum(xx) - zsig .* sigma : sigma : maximum(xx) + zsig .*sigma
	yybins = minimum(yy) - zsig .* sigma : sigma : maximum(yy) + zsig .*sigma
	zzbins = minimum(zz) - zsig .* sigma : sigma : maximum(zz) + zsig .*sigma

	heigth  = 1000
	xxcontents = heigth * ones(Base.size(xx))

	coors = ndim == 2 ? (xx, yy) : (xx, yy, zz)
	edges = ndim == 2 ? (xxbins, yybins) : (xxbins, yybins, zzbins)
	hh  = SB.fit(SB.Histogram, coors, SB.weights(xxcontents), edges)

	centers = [(edge[2:end] + edge[1:end-1])./2 for edge in edges]

	factor = 10.
	sigma_kernel = (sigma * factor) .* ones(ndim)
	weights_ = IF.imfilter(hh.weights, IF.Kernel.gaussian(sigma_kernel))

	cells    = findall(x -> x .> threshold, weights_)
	coors    = Tuple([centers[i][cell[i]] for cell in cells] for i in 1:1:ndim)
	contents = [hh.weights[cell] for cell in cells]

	contents = [weights_[cell] for cell in cells]

	return (coors = coors, cells = cells, contents = contents, edges = edges)

end
