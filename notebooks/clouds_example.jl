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
import DataFrames as DF
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
import Graphs     as GG
import MetaGraphs as MG
import GraphPlot  as GP
end

# ╔═╡ a57cdb41-c388-4976-bec8-ec0650fb139c
import Clouds as jc

# ╔═╡ cdc50171-b288-40b6-9d0d-9511901218e0
md"""
## Description
Clouds example - a line in 2D/3D 
J.A. Hernado,
Santiago, September 2022
---
"""

# ╔═╡ 3922eba2-f322-4b06-b9e0-83bc723d7930
PlutoUI.TableOfContents(title = "Clouds in Julia (dev)", indent = true, aside = true)

# ╔═╡ 3aedeb39-f255-4fd5-9ac3-29888a129e90
plotly();

# ╔═╡ 7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
md"""
## Generate Image
Produces a smeared curve in 2D or ·D
"""

# ╔═╡ e8848fd9-205e-4b56-b192-62f1acda8d7e
begin
bndim = @bind nndim Select([2, 3])

#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))

md"""
Select dimensions of the line $(bndim)
"""
end

# ╔═╡ 6c8bf138-8fec-4c69-b4dd-4284faddeed0
begin
img = jc.line(ndim = nndim, threshold = 6.)
end;

# ╔═╡ 1a8e9aa9-a47d-40fd-84c6-cfa49f9b1cc4
begin

md"""
**Image**
"""
end

# ╔═╡ 5a1832c1-33ff-45dc-8f47-212179dbe862
md"""
## Clouds
"""

# ╔═╡ 7cc12053-2bb8-4f24-b1a3-81c7fb679e19
begin
bcellnode = @bind xcellnode Select([:bygradient, :bycell])

#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))

md"""
Select the mode nodes are build: using the gradient or create a node for each cell $(bcellnode)
"""
end

# ╔═╡ f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
begin
cellnode = xcellnode == :bycell
steps = Tuple(edge[2]-edge[1] for edge in img.edges)
xcl, xnd, graph, edges  = jc.clouds(img.coors, img.contents, steps; cellnode = cellnode)
end;

# ╔═╡ 4e43c8e3-89e2-44ca-a6ed-48a364d90486
begin
md"""
steps of the voxels: $(steps[1])
"""
end

# ╔═╡ 13ac9fdf-46d0-4940-80e3-8619f0609108
md"""
## Plots
"""

# ╔═╡ a689debb-8763-45c4-a03d-94c8e970b243
begin

blabel = @bind label Select([:contents, :grad, :lap, :curmax, :curmin, :node, :nborders, :cloud])

#blabel = @bind typeevt Select(coll(:contents, :grad, :lap, :curmin, :curmax, :nodes, :nbordes))

md"""
Select label to plot $(blabel)
"""
end

# ╔═╡ 8dca9736-1140-495c-98a3-4cb5acc8ffc1
begin
vals = getfield(xcl, label)
minv, maxv = minimum(vals), maximum(vals)
end;

# ╔═╡ f17d0274-4a61-423c-a76f-870dcef41a60
begin
brange0 = @bind v0 Slider(minv:maxv, default = minv)
brange1 = @bind v1 Slider(minv:maxv, default = maxv)
md"""
Selec range for variable $(label):
minimum $(brange0)
maximum  $(brange1)
"""
end

# ╔═╡ e7544908-23e0-4e3a-ad93-2af5e0dc11f1
md"""
Selected range : [ $(v0), $(v1) ]
"""

# ╔═╡ 2814ba8e-58fa-4b68-af7b-b9e6656dcc19
md"""
## Nodes

Number of nodes            : $(length(xnd.contents))
"""

# ╔═╡ 7b7981ca-1540-48a1-88e1-4f27e7787b70
md"""
## Graph
"""

# ╔═╡ 1c402508-afd3-46a1-8dbc-a23fd9bd63e1
begin
GP.gplot(graph, nodelabel=1:GG.nv(graph), nodesize = xnd.contents)
end

# ╔═╡ 53ed70f6-8ed8-410b-8448-e6f3373240a7
md"""
#### Minimum Spanning Tree
"""

# ╔═╡ 98481aff-2884-4774-ba7e-f0df54c1c8ca
begin
mst = GG.Graph(GG.prim_mst(graph))
GP.gplot(mst, nodelabel = 1:GG.nv(mst), nodesize = xnd.contents)
end

# ╔═╡ a779ac6e-5bac-46f1-b8ef-1e3d5b111f4d
md"""
## Code
"""

# ╔═╡ dfa64554-5fb1-4d63-80d3-19aee7a476b8
begin
function cplot(cl, label, title, vrange, edges)
	ndim = length(cl.coors)
	vals = getfield(cl, label)
	mask = (vals .>= vrange[1]) .* (vals .<= vrange[2])
	coors = [c[mask] for c in cl.coors]
	vvals = vals[mask]
	theme(:dark)
	p1 = ndim == 2 ? histogram2d(coors..., weights = vvals, nbins = edges) : p1 = scatter(coors..., zcolor = vvals, alpha = 0.1)
	p2 = histogram(vvals, nbins = 100)
	plot(p1, p2, title = title)
end
end

# ╔═╡ d26c89ae-1629-4e98-8bde-3e8abe8bfd8d
cplot(img, :contents, :contents, [minimum(img.contents), maximum(img.contents)], img.edges)

# ╔═╡ 1fab453f-5dab-48bb-87d2-1c92b3f6d7cc
cplot(xcl, label, label, [v0, v1], edges)

# ╔═╡ 97fd042c-39f4-403e-a09c-6765bc94a1fb
md"""
## Dev area
"""

# ╔═╡ Cell order:
# ╟─5dcb2929-115e-459c-b98d-43ae7bcabd3a
# ╟─a9d9186f-19aa-41d7-8ec6-ad5197a74b8b
# ╠═a57cdb41-c388-4976-bec8-ec0650fb139c
# ╟─cdc50171-b288-40b6-9d0d-9511901218e0
# ╟─3922eba2-f322-4b06-b9e0-83bc723d7930
# ╟─3aedeb39-f255-4fd5-9ac3-29888a129e90
# ╟─7408b8f3-31f8-4764-bbf1-2bf9b245a8bb
# ╟─e8848fd9-205e-4b56-b192-62f1acda8d7e
# ╠═6c8bf138-8fec-4c69-b4dd-4284faddeed0
# ╟─1a8e9aa9-a47d-40fd-84c6-cfa49f9b1cc4
# ╠═d26c89ae-1629-4e98-8bde-3e8abe8bfd8d
# ╟─5a1832c1-33ff-45dc-8f47-212179dbe862
# ╟─7cc12053-2bb8-4f24-b1a3-81c7fb679e19
# ╟─f5dbdc6d-6676-4e0c-a70e-a5daafbbd9db
# ╟─4e43c8e3-89e2-44ca-a6ed-48a364d90486
# ╟─13ac9fdf-46d0-4940-80e3-8619f0609108
# ╟─a689debb-8763-45c4-a03d-94c8e970b243
# ╟─8dca9736-1140-495c-98a3-4cb5acc8ffc1
# ╟─f17d0274-4a61-423c-a76f-870dcef41a60
# ╟─e7544908-23e0-4e3a-ad93-2af5e0dc11f1
# ╠═1fab453f-5dab-48bb-87d2-1c92b3f6d7cc
# ╟─2814ba8e-58fa-4b68-af7b-b9e6656dcc19
# ╟─7b7981ca-1540-48a1-88e1-4f27e7787b70
# ╟─1c402508-afd3-46a1-8dbc-a23fd9bd63e1
# ╟─53ed70f6-8ed8-410b-8448-e6f3373240a7
# ╟─98481aff-2884-4774-ba7e-f0df54c1c8ca
# ╟─a779ac6e-5bac-46f1-b8ef-1e3d5b111f4d
# ╠═dfa64554-5fb1-4d63-80d3-19aee7a476b8
# ╠═97fd042c-39f4-403e-a09c-6765bc94a1fb
