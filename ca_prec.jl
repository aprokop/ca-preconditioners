#!/usr/bin/env julia
# Copyright 2017, Andrey Prokopenko, ORNL

include("meshnd.jl")
include("meshsparse.jl")
include("ca_utils.jl")

## Input data
# Note: the matrix is assumed to be symmetric. If that is not true, one would
# need to edit the underlap() function to properly compute the first underlap
# set (other sets are based on this one)
global m
if true
    ## Matrix parameters
    matrix_source = "generator"
    # matrix_source = 'file'

    if matrix_source == "generator"
        # Domain and subdomains parameters
        #  nxyz - number of elements in array sets the problem dimension.
        #  mxyz - number of subdomains in each direction
        #         Note that this array can be larger than the nxyz array, in
        #         which case we reduce its size.
        nxyz = [4 4]
        mxyz = [1 1]

    elseif matrix_source == "file"
        # Number of subdomains
        m = 16
        # matrix_file = 'G3_circuit/G3_circuit.mtx'
        # matrix_file = 'Boeing/msc04515.mtx'
        matrix_file = "A.mm"
    end
end

global ilu_setup
if true
    ## Solver parameters
    # Outer solver parameters
    # NOTE: RAS preconditioners cannot be used with CG (due to nonsymmetry)

    # solver = 'CG'
    # solver = 'GMRES-left'
    solver = "GMRES-right"

    max_its = 1000
    tol     = 1e-12
    restart = 20

    # Subdomain solver parameters
    sub_solver = "direct"

    # sub_solver = 'ilu'
    if sub_solver == "ilu"
        # ilu_setup.type    = "crout"
        ilu_setup.type    = "nofill"
        ilu_setup.milu    = "row"
        ilu_setup.droptol = 0.1
    end
end

## Check input parameters consistency
if matrix_source == "generator"
    # Check that the dimension of subdomain array is at least that of nxyz
    # array (if smaller, last dimension is not used)
    @assert(size(mxyz,2) >= size(nxyz,2))
    mxyz = mxyz[1:size(nxyz,2)]

    # Check that there is at least one point per dimension per subdomain
    # FIXME Julia complains about the assert
    # @assert(size(find(mxyz > nxyz),2) == 0, "Elements of mxyz must be <= than corresponding elements of nxyz")

    m = prod(mxyz)
end

## Construct matrix and rhs
global A
global coords
if matrix_source == "generator"
    coords = zeros(prod(nxyz), 3)

    if size(nxyz,2) == 1
        nx = nxyz[1]
        A = meshsparse(meshnd(nx))
        coords[:,1] = [1:nx]

    elseif size(nxyz,2) == 2
        nx = nxyz[1]
        ny = nxyz[2]

        # Construct the matrix
        # NOTE: we swap x and y dimensions as meshnd uses YX ordering
        A = meshsparse(meshnd(ny, nx), 5);  # 5-point stencil
        # A = meshsparse(meshnd(ny, nx), 9);  # 9-point stencil

        # Construct coordinates
        for j = 1:ny
            for i = 1:nx
                coords[(j-1)*nx + i,:] = [i j 0]
            end
        end

    else
        nx = nxyz[1]
        ny = nxyz[2]
        nz = nxyz[3]

        # Construct the matrix
        # NOTE: we swap x and y dimensions as meshnd uses YXZ ordering
        A = meshsparse(meshnd(ny, nx, nz),  7); #  7-point stencil
        # A = meshsparse(meshnd(ny, nx, nz), 27); # 27-point stencil

        # Construct coordinates
        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    coords[(k-1)*ny*nx + (j-1)*nx + i,:] = [i j k]
                end
            end
        end
    end

elseif matrix_source == "file"
    # Read in the matrix
    A = mmread(matrix_file)

    # No coordinates are available
end
# Make sure the matrix is square
@assert(size(A,1) == size(A,2))

n = size(A,1)

# Set seed for reproducibility
# rng(1426)

b = A*rand(n,1)

## Construct partitioning
# row_gids0 and col_gids0 are both cell arrays, and follow Petra notion, i.e.
#   row_gids{I} contains global DOFs of a row    map of subdomain I
#   col_gids{I} contains global DOFs of a column map of subdomain I
# gid2proc is the global ID -> subdomain ID mapping
global row_gids0
global col_gids0
global gid2proc
if matrix_source == "generator"
    row_gids0 = partition_manual(A, nxyz, mxyz)
elseif matrix_source == "file"
    row_gids0 = partition_metis(A, m)
end
col_gids0 = cell(1,m)
gid2proc  = zeros(1,n)
for I = 1:m
    col_gids0{I} = construct_cols(A, row_gids0{I})
    gid2proc(row_gids0{I}) = I
end


# ================ NOT CONVERTED STARTING THIS POINT ================
# % This is a debugging option
# % 0 : plot the residual history for selected preconditioners [default]
# % 1 : plot the connectivity structure of constructed preconditioner
# % 2 : plot the dropped connectivity structure
# global do_stencil_plot
#
# do_stencil_plot = 1
#
# %% Select preconditioners
# % Select a list of precondioners to run
# % Each run can be commented/uncommented independetly
# % WARNING: do not use split_type = 'full' with Schwarz overlap != 0
# precs = {}
# %% None
# % PROROTYPE: (A,name,{schwarz_level,schwarz_method},sub_solver,underlap_level)
# % precs = [precs {get_prec(A, 'none',                         {0},        sub_solver, 0)}]
# % precs = [precs {get_prec(A, 'subdomain Jacobi',             {0},        sub_solver, 0)}]
# % precs = [precs {get_prec(A, 'underlap + diagonal',          {0},        sub_solver, 1)}]
# % precs = [precs {get_prec(A, 'underlap + block diagonal',    {0},        sub_solver, 1)}]
# % precs = [precs {get_prec(A, 'underlap + block GS',          {0},        sub_solver, 1)}]
# % precs = [precs {get_prec(A, 'underlap + diagonal',          {0},        sub_solver, 2)}]
# % precs = [precs {get_prec(A, 'underlap + block diagonal',    {0},        sub_solver, 2)}]
# precs = [precs {get_prec(A, 'underlap + block GS',          {0},        sub_solver, 2)}]
# % precs = [precs {get_prec(A, 'l1 Gauss-Seidel',      {0},        sub_solver, {})}]
#
# %% Solve
# nprecs = size(precs,2)
# labels = cell(nprecs,1)
#
# % Reuse same colors for < 10 plots
# cc = linspecer(max(9,nprecs))
#
# for i = 1:nprecs
    # % Setup preconditioner
    # % Setting up a preconditioner is done by calling the preconditioner
    # % function with an empty array. This constructs all the necessary data for
    # % solve, and returns the preconditioner label which is going to be used in
    # % a legend
    # tic
    # label = precs{i}([])
    # setup_time = toc
#
    # if ~do_stencil_plot
        # disp(['Running: ', label])
        # fprintf('  setup time = %f\n', setup_time)
#
        # tic
        # if strcmp(solver, 'CG')
            # [~,flag,~,~,res_hist] = pcg(A, b, tol, max_its, precs{i})
        # else
            # rst    = min(restart,                size(A,1))
            # mx_its = min(fix(max_its/restart)+1, size(A,1))
            # if strcmp(solver, 'GMRES-left')
                # [~,flag,~,~,res_hist] = gmres(A, b, rst, tol, mx_its, precs{i})
            # elseif strcmp(solver, 'GMRES-right')
                # params = [tol, max_its, 1]
                # x0 = zeros(n, 1)
                # [~,flag,~,res_hist] = gmres_r(A, precs{i}, b, params, x0)
            # end
        # end
        # solve_time = toc
#
        # % Check result for convergence
        # if flag
            # fprintf('%s has not converged: flag = %d\n', solver, flag)
            # continue
        # else
            # fprintf('  solve time = %f\n', solve_time)
        # end
#
        # labels{i} = label
#
        # % Plot residual history
        # semilogy(1:size(res_hist), res_hist/norm(b), '-o', ...
            # 'Color',     cc(i,:), ...
            # 'LineWidth', 3)
        # hold on
    # end
# end
#
# if ~do_stencil_plot
    # % Construct title
    # if strcmp(matrix_source, 'generator')
        # if     size(nxyz,2) == 1, titlestr{1} = 'Laplace1D: '
        # elseif size(nxyz,2) == 2, titlestr{1} = 'Laplace2D: '
        # else                      titlestr{1} = 'Laplace3D: '
        # end
        # titlestr{1} = [titlestr{1}, 'size=[', int2str(nxyz), ']', ...
            # ', subdomains=[', int2str(mxyz), ']']
    # elseif strcmp(matrix_source, 'file')
        # titlestr{1} = [matrix_file, ': size=', int2str(size(A,1)), ...
            # ', subdomains=', int2str(m)]
    # end
    # titlestr{2} = solver
    # if strcmp(solver, 'GMRES-left')
        # titlestr{2} = [titlestr{2}, ': restart=', int2str(restart)]
    # elseif strcmp(solver, 'GMRES-right')
        # titlestr{2} = [titlestr{2}, ': no restart']
    # end
#
    # title(titlestr)
    # xlabel('Iteration number')
    # ylabel('Relative residual')
    # legend(labels)
    # % set(gca,'ygrid','on')
# end
# hold off
#
# end
#
# %% Helper functions
# % underlap() constructs level s underlap
# % It can operator in two modes:
# %   - 'full'  where level 1 underlap contains subdomain nodes within distance 1 from *any* neighbor
# %   - 'split' where 'full' underlap set is split into subsets (the article mode)
# function [sets] = underlap(A, i, s, stype)
# global m
# global gid2proc
# global row_gids0
# global col_gids0
#
# assert(s > 0)
# assert(strcmp(stype, 'full') || strcmp(stype, 'split'))
#
# % Get subdomain row and column maps
# row_gids = row_gids0{i}
# col_gids = col_gids0{i}
#
# % Construct local matrix (equivalent to the distributed variant in Trilinos)
# % It is a rectangular matrix (short and fat)
# lA = A(row_gids, col_gids)
# n  = size(lA,1)
#
# % excl set is used to exclude points which were already identified to be a
# % part of some s-set. We need it as, for instance, underlap level 2 set is
# % connected to both underlap level 1 and underlap level 3. At that point,
# % excl would contain underlap level 1, so we would get underlap level 3.
# excl = []
#
# if strcmp(stype, 'full')
    # [~,cols] = find(lA)
    # cols = setdiff(cols, 1:n)
    # % Find rows with nonzero elements in ghost columns
    # [rows,~] = find(lA(:,cols))
    # rows = unique(rows')
#
    # sets = cell(s,1)
#
    # sets{1} = rows
#
    # excl = sets{1}
    # for l = 2:s
        # if size(excl, 1) == 0
            # sets{l} = []
            # continue
        # end
#
        # % Find rows connected to the last s-set
        # [rows,~] = find(lA(:, sets{l-1}))
        # rows = unique(rows')
#
        # % Construct a new s-set
        # sets{l} = setdiff(rows, excl)
#
        # % Update excl set to include current s-set
        # excl = cat(2, excl, sets{l})
    # end
#
# else % stype = 'split'
    # % Find subdomains within distance s of the current one
    # overlap_gids = partition_schwarz(A, row_gids, s)
    # procs        = gid2proc(overlap_gids)
    # neighs       = setdiff(unique(procs), i)
#
    # %% Construct matrix D
    # % Matrix D is such that D(i,j) is the distance from node i to subdomain j
    # % 1. We consider only neighbor subdomains, the distance to non-neighbors is ignored
    # % 2. D(i,j) = 0 means infinity
    # D = zeros(n,m)
    # for k = 1:size(neighs,2)
        # neigh = neighs(k)
#
        # rows = row_gids0{neigh}
        # for l = 1:s
            # % Find "partition_schwarz(A, row_gids0{k}, l) \ partition_schwarz(A, row_gids0{k}, l-1)"
            # rows1 = partition_schwarz(A, rows, 1)
            # sub_k_over_l_gids = setdiff(rows1, rows)
            # rows = rows1
#
            # % Find local indices of the subdomain that correspond to the
            # % overlapped part of the "neighboring subdomain"
            # sub_i_under_gids  = setdiff(row_gids0{i}, setdiff(row_gids0{i}, sub_k_over_l_gids))
            # sub_i_under_lids  = []
            # for t = 1:size(sub_i_under_gids,2)
                # sub_i_under_lids = [sub_i_under_lids find(row_gids0{i} == sub_i_under_gids(t))]
            # end
#
            # % Fill in correspoding rows of D with distance l
            # D(sub_i_under_lids, k) = l
        # end
    # end
#
    # %% Construct sets
    # % sets{., 1} : contains a level subset (for varios underlap levels)
    # % sets{., 2} : contains distances from each node in the level subset to
    # %              each of the neighboring subdomains
    # sets = cell(0,2)
#
    # % S is a matrix with each row being a unique row of D
    # S = unique(D, 'rows')
    # pos = 1
    # for t = 1:size(S,1)
        # distances = S(t,:)
        # if size(find(distances),2) == 0
            # % Row with all zeros. Nodes corresponding to such rows are futher
            # % away from any neighbor than distance s, therefore we ignore them
            # continue
        # end
#
        # % Level subset = all rows of D with the same distances as the row of S
        # sets{pos,1} = find(ismember(D, distances, 'rows'))'
#
        # % We replace 0 with inf so that we have a partial sorting later
        # distances(distances == 0) = inf
        # sets{pos,2} = distances
#
        # pos = pos+1
    # end
# end
# end
#
# % Standard block Jacobi method
# function [M, row_gids] = prec_block_jacobi(A, sub_index, schwarz_level)
# global row_gids0
#
# if size(A,1) == 0
    # M = 'Block Jacobi'
    # return
# end
#
# row_gids = partition_schwarz(A, row_gids0{sub_index}, schwarz_level)
# M = A(row_gids, row_gids)
# end
#
# % Standard unpreconditioned method (do nothing)
# function [M, row_gids] = prec_none(A, sub_index, schwarz_level)
# global row_gids0
#
# if size(A,1) == 0
    # M = 'None'
    # return
# end
# assert(schwarz_level == 0)
#
# row_gids = partition_schwarz(A, row_gids0{sub_index}, schwarz_level)
# M = eye(size(row_gids,2))
# end
#
# function [M, row_gids] = prec_l1_gs(A, sub_index, schwarz_level)
# global row_gids0
#
# if size(A,1) == 0
    # M = 'l1 Gauss-Seidel'
    # return
# end
# assert(schwarz_level == 0)
#
# row_gids = partition_schwarz(A, row_gids0{sub_index}, schwarz_level)
# M = tril(A(row_gids, row_gids))
#
# end
#
# % SC14 paper method: interior + diagonalized underlap
# function [M, row_gids] = prec_underlap_diag(A, sub_index, schwarz_level, s)
# global m
# global row_gids0
#
# assert(schwarz_level <= 0, ...
    # 'schwarz_level > 0 was most likely incorrect for underlap and was removed')
#
# if size(A,1) == 0
    # M = ['Underlap (', int2str(s), ') + diagonal']
    # return
# end
# assert(s > 0)
#
# % Construct overlapped matrix
# row_gids = row_gids0{sub_index}
# M = A(row_gids, row_gids)
#
# % Construct underlap sets
# uset = underlap(A, sub_index, s, 'full')
#
# % Fix interior part
# n = size(row_gids0{sub_index},2)
# inner = 1:n
# for i = 1:s
    # inner = setdiff(inner, uset{i})
# end
# outer = setdiff(1:n, inner)
#
# M(inner, outer) = 0
# M(outer, inner) = 0
# M(outer, outer) = diag(diag(M(outer,outer)))
#
# end
#
# % New method 1: inner + block diagonal underlap
# function [M, row_gids] = prec_underlap_blockdiag(A, sub_index, schwarz_level, s)
# global m
# global row_gids0
#
# assert(schwarz_level <= 0, ...
    # 'schwarz_level > 0 was most likely incorrect for underlap and was removed')
#
# if size(A,1) == 0
    # M = ['Underlap (', int2str(s), ') + block diagonal']
    # return
# end
# assert(s > 0)
#
# % Construct overlapped matrix
# row_gids = row_gids0{sub_index}
# M = A(row_gids, row_gids)
#
# % Construct underlap sets
# uset       = underlap(A, sub_index, s, 'full')
# uset_split = underlap(A, sub_index, s, 'split')
#
# % Fix interior part
# n = size(row_gids0{sub_index},2)
# inner = 1:n
# for i = 1:s
    # inner = setdiff(inner, uset{i})
# end
#
# for i = 1:s
    # % Zero out bi-directional connections between inner and underlap
    # M(inner,   uset{i}) = 0
    # M(uset{i}, inner) = 0
    # % Zero out bi-directional connections between different underlap levels
    # for j = i+1:s
        # M(uset{i}, uset{j}) = 0
        # M(uset{j}, uset{i}) = 0
    # end
#
    # % Zero out bi-directional connections between subsets within the same underlap level
    # nsets = size(uset_split, 1)
    # for j = 1:nsets
        # set_j = uset_split{j,1}
        # for k = 1:nsets
            # if k == j
                # continue
            # end
            # set_k = uset_split{k,1}
#
            # M(set_j,set_k) = 0
            # M(set_k,set_j) = 0
        # end
    # end
# end
#
# end
#
# % New method 1: inner + block GS underlap
# function [M, row_gids] = prec_underlap_blockgs(A, sub_index, schwarz_level, s)
# global m
# global row_gids0
#
# assert(schwarz_level <= 0, ...
    # 'schwarz_level > 0 was most likely incorrect for underlap and was removed')
#
#
# if size(A,1) == 0
    # M = ['Underlap (', int2str(s), ') + block GS']
    # return
# end
# assert(s > 0)
#
# % Construct overlapped matrix
# row_gids = row_gids0{sub_index}
# M = A(row_gids, row_gids)
#
# % Construct underlap sets
# uset       = underlap(A, sub_index, s, 'full')
# uset_split = underlap(A, sub_index, s, 'split')
#
# % Fix interior part
# n = size(row_gids0{sub_index},2)
# inner = 1:n
# for i = 1:s
    # inner = setdiff(inner, uset{i})
# end
# % M(inner,inner) = tril(M(inner,inner))
#
# for i = 1:s
    # % Zero out uni-directional connections from underlap to inner
    # M(uset{i}, inner) = 0
#
    # % Zero out uni-directional connections between all underlap subsets, and between underlap sets and inner
    # nsets = size(uset_split, 1)
    # for j = 1:nsets
        # set_j = uset_split{j,1}
        # ind_j = uset_split{j,2}
        # for k = 1:nsets
            # if k == j
                # continue
            # end
            # set_k = uset_split{k,1}
            # ind_k = uset_split{k,2}
#
            # % This is a very tricky part
            # % We zero out connections based on a partial sorting algorithms based on array operator <
            # % If we can compare two distance arrays, and one is smaller than
            # % the other, then we zero out connections. An example:
            # % (1,1) is     smaller than (1,2)   => we do     zero out
            # % (1,2) is NOT smaller than (1,1)   => we do not zero out
            # % (1,2) is NOT smaller than (2,1)   => we do not zero out
            # % (1,2) is     smaller than (2,2)   => we do     zero out
            # % (1,2) is     smaller than (2,inf) => we do     zero out
            # % etc
            # [cmp1, cmp2] = compare(ind_j, ind_k)
            # if cmp1 || ...          % ind_j <= ind_k
               # (~cmp1 && ~cmp2)     % ind_j and ind_k are incomparable (in the partial sorting sense)
                # M(set_j,set_k) = 0
            # end
        # end
    # end
# end
#
# end
#
# function [x] = gen_prec(A, schwarz_level, schwarz_method, solver, P, b)
# global ilu_setup
# global m
# global coords
# global do_stencil_plot
# global row_gids0
#
# assert(strcmp(schwarz_method, 'RAS') || strcmp(schwarz_method, 'AS') || strcmp(schwarz_method, 'AAS'))
# assert(strcmp(solver, 'direct') || strcmp(solver, 'ilu'))
#
# persistent M
# persistent L
# persistent U
# persistent row_gids
#
# if size(b,1) == 0
    # % Setup stage
    # row_gids = cell(1,m)
#
    # colors  = ['k', 'b', 'g', 'r', 'm', 'c']
    # ncolors = size(colors,2)
#
    # for I = 1:m
        # tic
        # [M{I}, row_gids{I}] = P(A, I, schwarz_level)
        # if do_stencil_plot
            # hold on
            # lcoords = coords(row_gids{I},:)
            # if do_stencil_plot == 1
                # % Show remaining connections
                # X = M{I}
            # else
                # % Show dropped connections
                # X = A(row_gids{I},row_gids{I}) - M{I}
            # end
            # color = colors(rem(I,ncolors)+1)
            # plot_stencil(X, lcoords, color)
            # % To look at the individual subdomains, uncomment the following lines
            # % Then, step through subdomains in command line with dbcont
            # if 0
                # keyboard
                # if I ~= m
                    # clf
                # end
            # end
        # end
#
        # if strcmp(solver, 'direct')
            # [L{I}, U{I}] = lu(M{I})
        # else
            # % Current approach is not very good, as it does
            # % ILU on the overlapped matrix
            # % What we would like to do is to solve "tails"
            # % exactly and then perform ILU only on the proper
            # % inner part (underlap s)
            # [L{I}, U{I}] = ilu(M{I}, ilu_setup)
        # end
    # end
#
    # % Construct label by attaching Schwarz to the label of P
    # x = P([], -1, -1)
    # if schwarz_level
        # x = [x, ', schwarz=', schwarz_method, ...
            # '(', int2str(schwarz_level), ')']
    # end
    # x = [x, ', solver=', solver]
#
# else
    # % Solve stage
    # n = size(b, 1)
    # x         = zeros(n, 1)
    # indicator = zeros(n, 1)
#
    # for I = 1:m
        # lb = b(row_gids{I})
#
        # lx = U{I} \ (L{I} \ lb)
#
        # if strcmp(schwarz_method, 'RAS')
            # % RAS (Restrictive Additive Schwarz)
            # % Assign only owned DOFs
            # x(row_gids0{I}) = lx(1:size(row_gids0{I},2))
#
        # elseif strcmp(schwarz_method, 'AS') || strcmp(schwarz_method, 'AAS')
            # % AS (Additive Schwarz)
            # % Contributions from multiple subdomains are summed
            # x(row_gids{I}) = x(row_gids{I}) + lx
            # indicator(row_gids{I}) = indicator(row_gids{I})+1
        # end
    # end
    # if strcmp(schwarz_method, 'AAS')
        # % AAS (Averaged Additive Schwarz)
        # % Contributions from multiple subdomains to the same dof are averaged
        # % by the number of contributing subdomain
        # x = x ./ indicator
    # end
# end
# end
#
# function [prec] = get_prec(A, name, schwarz_cell, sub_solver, underlap_level)
# if size(schwarz_cell,2) == 1
    # if schwarz_cell{1} == 0
        # schwarz_cell = [schwarz_cell{1} {'AS'}]
    # else
        # error('Please specify RAS/AS for overlap level > 0')
    # end
# end
#
# if     strcmp(name, 'none')
    # prec = @(b)(gen_prec(A, ...
        # schwarz_cell{1}, schwarz_cell{2}, sub_solver, ...
        # @(X,i,l)(prec_none(X,i,l)), b))
#
# elseif strcmp(name, 'subdomain Jacobi')
    # prec = @(b)(gen_prec(A, ...
        # schwarz_cell{1}, schwarz_cell{2}, sub_solver, ...
        # @(X,i,l)(prec_block_jacobi(X,i,l)), b))
#
# elseif strcmp(name, 'l1 Gauss-Seidel')
    # prec = @(b)(gen_prec(A, ...
        # schwarz_cell{1}, schwarz_cell{2}, sub_solver, ...
        # @(X,i,l)(prec_l1_gs(X,i,l)), b))
#
# elseif strcmp(name, 'underlap + diagonal')
    # prec = @(b)(gen_prec(A, ...
        # schwarz_cell{1}, schwarz_cell{2}, sub_solver, ...
        # @(X,i,l)(prec_underlap_diag(X, i, l, underlap_level)), b))
#
# elseif strcmp(name, 'underlap + block diagonal')
    # prec = @(b)(gen_prec(A, ...
        # schwarz_cell{1}, schwarz_cell{2}, sub_solver, ...
        # @(X,i,l)(prec_underlap_blockdiag(X, i, l, ...
        # underlap_level)), b))
#
# elseif strcmp(name, 'underlap + block GS')
    # prec = @(b)(gen_prec(A, ...
        # schwarz_cell{1}, schwarz_cell{2}, sub_solver, ...
        # @(X,i,l)(prec_underlap_blockgs(X, i, l, ...
        # underlap_level)), b))
#
# else
    # error(['Unknown preconditioner: ', name])
# end
# end
