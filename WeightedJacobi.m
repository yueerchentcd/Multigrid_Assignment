

function x = WeightedJacobi(A, b, x, omega, num_steps)
%WEIGHTEDJACOBI Applies weighted Jacobi iterations.
%
%   x = WeightedJacobi(A, b, x, omega, num_steps)
%
%   Solves approximately Ax = b using num_steps iterations of
%   weighted Jacobi:
%
%       x_new = x + omega * D^{-1} * (b - A*x)
%
%   where D is the diagonal of A.

    D = diag(A);

    for k = 1:num_steps
        r = b - A*x;
        x = x + omega * (r ./ D);
    end

end

