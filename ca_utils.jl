# Copyright 2017, Andrey Prokopenko, ORNL

function partition_manual(A, nxyz, mxyz)
    dim = size(nxyz, 2)
    @assert(dim <= 3)

    if dim == 1
        nxyz = [nxyz[1] 0 0]
        mxyz = [mxyz[1] 1 1]
    elseif dim == 2
        nxyz = [nxyz[1] nxyz[2] 0]
        mxyz = [mxyz[1] mxyz[2] 1]
    end
    m = prod(mxyz)

    row_gids = Array{Any}(1,m)
    for K = 1:mxyz[3]
        for J = 1:mxyz[2]
            for I = 1:mxyz[1]
                IJK = (K-1)*mxyz[2]*mxyz[1] + (J-1)*mxyz[1] + I

                (min_x, max_x) = sub1D_bnd(nxyz[1], mxyz[1], I)

                if nxyz[2] != 0
                    (min_y, max_y) = sub1D_bnd(nxyz[2], mxyz[2], J)
                else
                    min_y = 1
                    max_y = 1
                end

                if nxyz[3] != 0
                    (min_z, max_z) = sub1D_bnd(nxyz[3], mxyz[3], K)
                else
                    min_z = 1
                    max_z = 1
                end

                # Calculate inner part othe subdomain
                rows = zeros(Int64, 1, (max_x-min_x+1)*(max_y-min_y+1)*(max_z-min_z+1))
                ind  = 1
                for k = min_z:max_z
                    for j = min_y:max_y
                        for i = min_x:max_x
                            rows(ind) = (k-1)*nxyz[2]*nxyz[1] +
                                        (j-1)*nxyz[1] + i
                            ind = ind+1
                        end
                    end
                end

                row_gids[IJK] = rows
            end
        end
    end

    return row_gids
end

# function [row_gids] = partition_metis(A, m)
# % Requires metismex MATLAB interface to METIS
# % Download here: https://github.com/dgleich/metismex
# map = metismex('PartGraphRecursive', A, m)
#
# row_gids = cell(1,m)
# for I = 1:m
    # row_gids{I} = find(map == (I-1))
# end
# end
#
# function [row_gids] = partition_schwarz(A, row_gids0, schwarz)
#
# % Add overlapped parts
# % The nodes are sorted with respect to the overlap level:
# % first go inner subdomain nodes, then overlap one, then overlap
# % two, and so on
# rows = row_gids0
# for s = 1:schwarz
    # % Find nodes connected to the previous overlap level,
    # % starting from 0
    # [~,cols] = find(A(rows,:))
    # if size(cols,1) == 1
        # cols = cols'
    # end
    # % Update nodes to include those.
    # rows = cat(2, rows, setdiff(cols', rows))
# end
#
# row_gids = rows
# end

function construct_cols(A, row_gids)
    # col_gids is similar to column map
    # It contains a list of columns with nonzero elements in row_gids rows, and
    # is ordered so that the part corresponding to row_gids goes first.
    (_, cols) = findn(A[row_gids,:])
    if size(cols,1) == 1
        cols = cols'
    end
    col_gids = cat(2, row_gids, setdiff(unique(cols'), row_gids))

    return col_gids
end

# function [neigh] = get_neigh_list(sub_index)
# global col_gids0
# global gid2proc
# procs = gid2proc(col_gids0{sub_index})
# neigh = setdiff(unique(procs), sub_index)
# end

function sub1D_bnd(n, m, i)
    @assert(i > 0 && i <= m, "Incorrect i value")

    min_width = trunc(Int64, n/m)
    num_wides = n - m*min_width

    if i <= num_wides
        min = (i-1) *   (min_width+1)
        max = min +      min_width+1
    else
        min = num_wides*(min_width+1) + (i-1 - num_wides)*min_width
        max = min +      min_width
    end
    min = min+1

    return min, max
end

# function [] = plot_stencil(A, coords, color)
# n = size(A,1)
#
# % 1. Separate the matrix into two parts: one with symmetric pattern, and a remainder
# % Get the square part of A (A may be short and wide rectangular matrix)
# S = A(1:n,1:n)
# Z = S
# % Create a pattern matrix
# Z(Z ~= 0) = 1
# % Create the symmetric subpattern
# [r,c] = find((Z+Z') == 2)
# % Create final pattern
# Z = sparse(r,c,ones(1,size(r,1)),n,n)
# S = S .* Z
#
# line_width = 3
# marker = 'o'
#
# % 2. Plot the symmetric part
# [rows,cols] = find(S)
# for k = 1:size(rows)
    # i = rows(k)
    # j = cols(k)
    # if i <= j % need to plot only once
        # plot3([coords(i,1), coords(j,1)], ...
              # [coords(i,2), coords(j,2)], ...
              # [coords(i,3), coords(j,3)], ...
              # 'LineWidth',          line_width, ...
              # 'Marker',             marker, ...
              # 'MarkerFaceColor',    color, ...
              # 'Color',              color)
    # end
# end
#
# % 3. Plot the nonsymmetric part
# [rows,cols] = find(A-S)
# scaling = 0.73
# u = coords(cols,1)-coords(rows,1)
# v = coords(cols,2)-coords(rows,2)
# w = coords(cols,3)-coords(rows,3)
# quiver3(coords(rows,1), coords(rows,2), coords(rows,3), ...
        # scaling*u, scaling*v, scaling*w, ...
        # 'LineWidth',        line_width, ...
        # 'Marker',           marker, ...
        # 'MarkerFaceColor',  color, ...
        # 'AutoScale',        'off', ...
        # 'Color',            color)
#
# axis square
# % axis off
# view(0,90)
#
# set(gcf, 'PaperPositionMode',   'auto')
# set(gcf, 'PaperPosition',       [0 0 100 100])
#
# end
#
# function [alessb, blessa] = compare(ind1, ind2)
# cmp = [ind1 <= ind2]
# alessb = ~sum(find(cmp == 0))
#
# cmp = [ind2 <= ind1]
# blessa = ~sum(find(cmp == 0))
#
# end
