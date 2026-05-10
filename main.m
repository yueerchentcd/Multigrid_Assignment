clear; clc; close all;

% ============================================================
% Question 1: Construct the 2D Poisson matrix
% ============================================================

p = 5;
n = 2^p - 1;

A = Poisson2D(n);
N = n^2;

fprintf('n = %d\n', n);
fprintf('Size of A = %d by %d\n', size(A,1), size(A,2));
fprintf('Number of nonzeros = %d\n', nnz(A));

% Plot sparsity pattern of the 2D Poisson matrix
figure;
spy(A);
title('Sparsity pattern of 2D Poisson matrix');
xlabel('Column index');
ylabel('Row index');
grid on;

% Save the figure for the report
saveas(gcf, 'sparsity_pattern.png');


% ============================================================
% Question 3: Test multigrid V-cycle convergence
% ============================================================

rng(1);
f = randn(N, 1);

u = zeros(N, 1);

direct_N = 9;
weight = 4;
pre_smooth = 3;
pos_smooth = 3;

[A_list, R_list, max_level] = Restrict(A, direct_N);

fprintf('Number of multigrid levels = %d\n', max_level);

num_cycles = 12;
res_norm = zeros(num_cycles + 1, 1);

% Residual norm before any V-cycle
res_norm(1) = norm(f - A*u, 2);

% Apply V-cycles
for k = 1:num_cycles
    u = Vcycle(1, A_list, R_list, f, u, direct_N, weight, pre_smooth, pos_smooth);
    res_norm(k+1) = norm(f - A*u, 2);
end

fprintf('\nResidual norms:\n');
for k = 0:num_cycles
    fprintf('V-cycle %2d: %.6e\n', k, res_norm(k+1));
end

% Plot residual convergence
figure;
semilogy(0:num_cycles, res_norm, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Number of V-cycles');
ylabel('Residual norm ||f - Au||_2');
title('Multigrid V-cycle convergence');
grid on;

% Save the figure for the report
saveas(gcf, 'residual_convergence.png');
