module Clouds

#import movements
include("moves.jl")
export moves, Moves

#import discretize and gradient functions
include("gradient.jl")
export discretize, discrete_gradient

# import and export the relevant functions and strucs
include("core.jl")
export clouds, nclouds

include("clustering.jl")
export cluster_nodes, clustering

include("spine.jl")
export spine

include("sources.jl")
export box2d, box3d, line


end
