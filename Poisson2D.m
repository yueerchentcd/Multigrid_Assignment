
function A = Poisson2D(n)
%POISSON2D Constructs the sparse 2D Poisson matrix.
%
%   A = Poisson2D(n) returns the sparse matrix corresponding
%   to the five-point finite difference discretisation of
%   -Delta u = f on (0,1)^2 with homogeneous Dirichlet
%   boundary conditions.
%
%   The grid has n interior points in each spatial direction,
%   so A has size n^2 by n^2.

    h = 1 / (n + 1);

    e = ones(n, 1);

    % 1D second-order finite difference matrix
    T = spdiags([-e, 2*e, -e], -1:1, n, n);

    % Identity matrix
    I = speye(n);

    % 2D five-point stencil matrix
    A = (kron(I, T) + kron(T, I)) / h^2;

end



