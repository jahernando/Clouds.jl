
import StatsBase     as SB
import LinearAlgebra as LA
import Graphs        as GG

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
function connection_table(cells)
    ndim   = length(cells[1])
    nsize  = length(cells)
    m      = moves(ndim)
    nmoves = length(m.moves)

    connection_table = zeros(Int64, nsize, nmoves)
    #aa = [[Tuple(icell .+ mi) for mi in m.moves] for icell in cells]

    function _ipos(icell)
        ids = findall(x -> x == icell, cells)
        return length(ids) != 1 ? 0 : ids[1]
    end

    for (i, icell) in enumerate(cells)
        for (j, imove) in enumerate(m.moves)
            kcell = Tuple(icell .+ imove)
            print(kcell)
            connection_table[i, j] = _ipos(kcell)
        end
    end
    return connection_table
end
