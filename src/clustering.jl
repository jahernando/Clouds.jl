
import StatsBase     as SB
import LinearAlgebra as LA
import Graphs        as GG

export cluster_nodes, clustering

VTI = Vector{Tuple{Vararg{Int64}}}


"""
return a connection table:
    the table are integer, the rows are the cells index and the columns the n-dim movements,
    the vlues of the tables at the i-cell and j-movements are the k-cell 
       that connect the i-cell and the j-movement, that is
       cell-k = cell-i + movement-j
    i.e cell--i, in the i-th position of cells, is (1, 1) and the movement-j is (1, 2)
    if the cell-k, in the k-th position of cells is (2, 2) therefore the value of
    the connetion matrix A[i, j] = k, if there is no cell (2, 2) in the cells list then A[i, j] = 0
        cells[i]    = (1, 1)
        movement[j] = (1, 1)
        A[i, j]     = k  
        if k > 0
            cells[k]    = (2, 2)
        if k = 0
            there is no (2, 2) in the list of *cells*
"""
function _connection_table(cells  ::VTI)
    ndim   = length(cells[1])
    nsize  = length(cells)
    m      = moves(ndim)
    nmoves = length(m.moves)

    connection_table = zeros(Int64, nsize, nmoves)

    function _ipos(icell)
        ids = findall(x -> x == icell, cells)
        return length(ids) != 1 ? 0 : ids[1]
    end

    for (i, icell) in enumerate(cells)
        for (j, imove) in enumerate(m.moves)
            kcell = Tuple(icell .+ imove)
            connection_table[i, j] = _ipos(kcell)
        end
    end
    return connection_table
end


function _filter_same_label(label ::Vector{Int64},
                            table ::Matrix{Int64})
    function _filter(id)
        i, j = id[1], table[id]
        if (i <= 0) | (j <= 0)
            return 0
       end
       return label[i] == label[j] ? j : 0
    end
    return _filter
end

function _filter_diff_label(label ::Vector{Int64},
                            table ::Matrix{Int64})
    function _filter(id)
        i, j = id[1], table[id]
        if (i <= 0) | (j <= 0)
            return 0
       end
       label[i] != label[j] ? j : 0
    end
    return _filter
end


function _filter_connection_table(table  ::Matrix{Int64},
                                 _filter ::Function)
    indices     = findall(table .> 0)
    xtable = zeros(Int64, size(table)...)
    xtable[indices] .= _filter.(indices)
    return xtable 
end


function _set_nodes(table ::Matrix{Int64})
    nsize = size(table)[1]
    ids   = 1:nsize
    nodes = zeros(Int64, nsize)
    visit = zeros(Bool , nsize)
    node     = 1
    nodes[1] = node
    function _next_node()
        sel = (nodes .> 0) .& (visit .== false)
        if (sum(sel) <= 0)
            sel = visit .== false
        end
        return (sum(sel) <= 0) ? 0 : ids[sel][1]
    end
    id = _next_node()
    while (id in ids)
        node = nodes[id] > 0 ? nodes[id] : node + 1
        nodes[filter(x -> x > 0, table[id, :])] .= node
        visit[id] = true
        id        = _next_node()
    end
    return nodes
end

function _get_links(nodes ::Vector{Int64},
                    table ::Matrix{Int64})
    nsize = size(table)[1]
    nnodes = maximum(nodes) 
    links = Dict{Int64, Vector{Int64}}()
    for i in 1:nnodes
        links[i] = []
    end
    for i in 1:nsize
        node = nodes[i]
        ids = table[i, :]
        ids = ids[ids.>0]
        if (length(ids) == 0)
            continue
        end 
        inodes = Set(nodes[ids])
        inodes = filter(x -> (x in links[node]) == false, inodes)
        for inode in inodes
            push!(links[node], inode)
        end
    end
    return links
end

"""
fron a vector *cells* (a vector of cartexian index in 2d or 3d) each one with a Int label, 
given by the *label*, a vector of Int, return a list of the same dimension of the cells,
with the association of each cell to a node.

A node is a collection of cells that are connected in the n-dimension by the a face, side or vertex,
and that share the same label.

i.e.
    inputs: 
        cells = [(1, 1), (1, 2), (3, 3), (3, 4)], label = [1, 1, 1, 2]
    returns:
        nodes = [1, 1, 2, 3]

"""
function cluster_nodes(cells ::VTI,
                       label ::Vector{Int64})
    ctable = _connection_table(cells)
    ntable = _filter_connection_table(ctable, _filter_same_label(label, ctable))
    nodes  = _set_nodes(ntable)
    return nodes
end    

"""
fron a vector *cells* (a vector of cartexian index in 2d or 3d) each one with a Int label, 
given by the *label*, a vector of Int, return a list of the same dimension of the cells,
with the association of each cell to a node, and a dictionary which keys are nodes and the valures
are the nodes associated to it.

A node is a collection of cells that are connected in the n-dimension by the a face, side or vertex,
and that share the same label.
Two nodes are linked in they have connecting cells (with a face, sice or vertex) each cell belonging
to different nodes.

i.e.
    inputs: 
        cells = [(1, 1), (1, 2), (3, 3), (3, 4)], label = [1, 1, 1, 2]
    returns:
        nodes = [1, 1, 2, 3]
        links = Dict(1 => [], 2 => [3], 3 => [2])

"""

function clustering(cells  ::VTI,
                    label  ::Vector{Int})
    ctable = _connection_table(cells)
    ntable = _filter_connection_table(ctable, _filter_same_label(label, ctable))
    ltable = _filter_connection_table(ctable, _filter_diff_label(label, ctable))
    nodes  = _set_nodes(ntable)
    links  = _get_links(nodes, ltable)
    return (nodes = nodes, edges = links)
end