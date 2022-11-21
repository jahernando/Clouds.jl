#-----
# Holds function and Struct for one-step movements in 2D and 3D
#-----

#using Base
#import StatsBase     as SB

export moves, Moves

#-----------------
# Data Types
#-----------------

# sort-cuts type definitions
N   = Vararg
T2I = Tuple{Int64, Int64}
TNI = Tuple{N{Int64}}
#TNV = Tuple{N{Float64}}
#TNN = Tuple{N{<:Number}}
VI  = Vector{Int64}
#VB  = Vector{Bool}
#VF  = Vector{Float64}
#VN  = Vector{<:Number}
#VAN = Vector{T} where T <: Array{<:Number}
#AN  = T where T <: Array{<:Number}
#AI  = T where T <: Array{Int64}


# Helper Data Type to hold one step movements in a 2D or 3D grid
# one steps 1D movements are : [1, 0], [0, 1], [1, 1], [-1, 0], ...
struct Moves
	moves  ::Tuple{N{VI}}   # list with the vector of unitary movements in 2D or 3D
	i0     ::Int64          # index of the null movement (0., 0) or (0, 0, 0)
	isym   ::Tuple{N{T2I}}  # list of the pair of indices of symmetric movements i.e (1, 0), (-1, 0)
	iortho ::Dict{T2I, VI}  # dictionary with the indices of the ortogonal movements
	                        # i.e movements (0, 1) (0, -1) are orthogornal to (1, 0), (-1, 0)
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

