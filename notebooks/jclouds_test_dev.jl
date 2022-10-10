### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5dcb2929-115e-459c-b98d-43ae7bcabd3a
using Pkg; Pkg.activate("/Users/hernando/work/investigacion/NEXT/software/julias/Clouds")

# ╔═╡ a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
begin
using Markdown
#using HDF5
#using DataFrames
import StatsBase as SB
using Plots
import LinearAlgebra as LA
#using Statistics
using PlutoUI
using Random
using Base
import Images.ImageFiltering as IF
#using TestImages
#using Distributions
import Graphs  as GG
import GraphPlot as GP
end

# ╔═╡ a57cdb41-c388-4976-bec8-ec0650fb139c
import Clouds as jc

# ╔═╡ cdc50171-b288-40b6-9d0d-9511901218e0
md"""

## Description


Test clouds 3D in Julia


J.A. Hernado,

Santiago, October 2022

---
"""

# ╔═╡ 3922eba2-f322-4b06-b9e0-83bc723d7930
PlutoUI.TableOfContents(title = "Clouds in Julia (dev)", indent = true, aside = true)

# ╔═╡ 3aedeb39-f255-4fd5-9ac3-29888a129e90
plotly();

# ╔═╡ 7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
md"""

## Generate Matrix

"""

# ╔═╡ e8848fd9-205e-4b56-b192-62f1acda8d7e
begin

bndim = @bind ndim Select([2, 3])

a  = 1
bb = @bind b Slider(2:10, default = a+1)
bc = @bind c Slider(b:10, default = a+2)
bd = @bind d Slider(c:10, default = a+3)

#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))

md"""

Select dimensions of the line $(bndim)

Select the rage of the box entries (a = 1, b, c, d)

b = $(bb)

"""
end

# ╔═╡ 5f0109a3-0e2d-4528-af51-652790a3c23e
md"""
c = $(bc)

"""

# ╔═╡ 836c7884-5a2c-49a1-aab0-8a124219e39e
if ndim == 3
	md"""
	d = $(bd)
	"""
end

# ╔═╡ c27a1991-24c6-4ead-9807-b2e2b37dfcdb
begin
if ndim == 2
	md"""
	Selected values (a = 1, b = $(b), c = $(c))
	"""
else
	md"""
	Selected values (a = 1, b = $(b), c = $(c), d = $d)
	"""
end
end

# ╔═╡ 416cb989-f4f9-4acd-b9de-8b754a15f79e
begin

function box2d(vals = [1, 2, 4])

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

	return (contents = mat, grad = grad, lap = lap, curmax = maxc, curmin = minc)

end



end #begin

# ╔═╡ 34820285-e673-4f45-9593-fc5cb409d3d1
begin


function box3d(vals = [1, 2, 3, 4])

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

end

# ╔═╡ 6a0e9248-957b-4d06-90a3-66cf4c6b54fe
begin
mat = ndim == 2 ? jc.box2d([1, b, c]) : jc.box3d([1, b, c, d])
end;

# ╔═╡ 2e39663f-be86-40e7-ab62-7fca98a4489f
mat.contents

# ╔═╡ 233ce695-1cc2-41ab-b798-4326fba77966
begin
aa   = mat.contents
xndim = length(size(aa))
steps = ones(xndim)
cells = CartesianIndices(aa)
coors = [vec([c[i] for c in cells]) for i in 1:xndim]
xcl = jc.clouds(coors, vec(aa), steps)
end;

# ╔═╡ 40f1d9a6-3f07-464d-8a97-bf627f6028e2
begin

function test(xcl, mat)

	@assert sum(isapprox.(xcl.contents, vec(mat.contents))) == 3^ndim
	@assert sum(isapprox.(xcl.grad, vec(mat.grad)))         == 3^ndim
	@assert sum(isapprox.(xcl.lap, vec(mat.lap)))           == 3^ndim
	@assert sum(isapprox.(xcl.curmax, vec(mat.curmax)))     == 3^ndim
	@assert sum(isapprox.(xcl.curmin, vec(mat.curmin)))     == 3^ndim

	return true
end

ok = test(xcl, mat)
md"""

**Test** of gradient, laplacian, curvatures ? : **$(ok)**

"""
end

# ╔═╡ Cell order:
# ╠═5dcb2929-115e-459c-b98d-43ae7bcabd3a
# ╟─a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
# ╟─a57cdb41-c388-4976-bec8-ec0650fb139c
# ╟─cdc50171-b288-40b6-9d0d-9511901218e0
# ╟─3922eba2-f322-4b06-b9e0-83bc723d7930
# ╟─3aedeb39-f255-4fd5-9ac3-29888a129e90
# ╟─7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
# ╟─e8848fd9-205e-4b56-b192-62f1acda8d7e
# ╟─5f0109a3-0e2d-4528-af51-652790a3c23e
# ╟─836c7884-5a2c-49a1-aab0-8a124219e39e
# ╟─c27a1991-24c6-4ead-9807-b2e2b37dfcdb
# ╟─416cb989-f4f9-4acd-b9de-8b754a15f79e
# ╟─34820285-e673-4f45-9593-fc5cb409d3d1
# ╟─6a0e9248-957b-4d06-90a3-66cf4c6b54fe
# ╠═2e39663f-be86-40e7-ab62-7fca98a4489f
# ╟─233ce695-1cc2-41ab-b798-4326fba77966
# ╟─40f1d9a6-3f07-464d-8a97-bf627f6028e2
