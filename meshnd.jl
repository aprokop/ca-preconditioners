# Copyright 2007, Timothy A. Davis, Univ. of Florida
# Copyright 2017, Andrey Prokopenko, ORNL

# Adapted from Tim Davis' meshnd code
function meshnd(m, n, k = 1)
    """MESHND creation and nested dissection of a regular 2D or 3D mesh.

    G = meshnd(m,n)   constructs a m-by-n 2D mesh G
    G = meshnd(m,n,k) constructs a m-by-n-by-k 3D mesh G

    Example:
    G = meshnd(4,5)

    returns
    G =
    1     2     3     4     5
    6     7     8     9    10
    11    12    13    14    15
    16    17    18    19    20

    With no inputs, a few example meshes are generated and plotted.

    See also nested, numgrid.
    """

    if k == 1
        # create the m-by-n-by-k mesh in "natural" (row-major) order.  This is how
        # a typical 2D mesh is ordered.  A column-major order would be better, since
        # in that case G(:) would equal 1:(m*n) ... but let's stick with tradition.

        # FIXME: this does not work as expected
        # G = reshape(1:(m*n*k), n, m, k)'
        G = zeros(Int64, m, n, k)
        for ii = 1:m
            for jj = 1:n
                G[ii, jj, 1] = (ii-1)*m + jj
            end
        end
    else
        # create the m-by-n-by-k mesh in column-major order. The first m-by-n-by-1
        # slice is in column-major order, followed by all the other slices 2 to k.
        G = reshape(1:(m*n*k), m, n, k)
    end

    return G

end
