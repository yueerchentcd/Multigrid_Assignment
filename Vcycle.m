

function x = Vcycle(level, A_list, R_list, b, x0, direct_n, weight, pre_smooth, pos_smooth)

% level       Current level, initial = 1
% A_list      Cell array of coefficient matrices on each level
% R_list      Cell array of restriction operators on each level
% b           Right-hand side
% x0          Initial guess
% direct_n    Threshold for directly solving Ax = b
% weight      Scaling coefficient between restriction and prolongation operators
% pre_smooth  Number of iterations of pre-smoothing
% pos_smooth  Number of iterations of post-smoothing

    % Load coefficient matrix on current level
    A = A_list{level};

    % Solve directly if problem is small enough
    n = size(b, 1);
    if n <= direct_n
        x = A \ b;
        return;
    end

    % Weighted Jacobi relaxation parameter
    omega = 2/3;

    % Pre-smoothing
    x = WeightedJacobi(A, b, x0, omega, pre_smooth);

    % Load restriction operator and construct interpolation operator
    R = R_list{level};
    P = weight * R';

    coarse_n = size(R, 1);

    % Compute residual and restrict it to the coarse grid
    r = b - A*x;
    r_H = R * r;

    % Recursively solve the coarse-grid error equation
    e0 = zeros(coarse_n, 1);
    e = Vcycle(level + 1, A_list, R_list, r_H, e0, direct_n, weight, pre_smooth, pos_smooth);

    % Interpolate coarse-grid error back and correct fine-grid solution
    x = x + P * e;

    % Post-smoothing
    x = WeightedJacobi(A, b, x, omega, pos_smooth);

end