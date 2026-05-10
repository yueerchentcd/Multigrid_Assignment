# Multigrid Assignment

This repository contains the MATLAB implementation for the Case Study 4 Multigrid Assignment.

## Files

- `main.m`  
  Main script for Questions 1--3. It constructs the 2D Poisson matrix, runs the multigrid V-cycle solver, records the residual norms, and generates the plots.

- `Poisson2D.m`  
  Constructs the sparse 2D Poisson matrix using the five-point finite difference stencil.

- `Restrict.m`  
  Builds the multigrid hierarchy and restriction operators.

- `Vcycle.m`  
  Implements the recursive multigrid V-cycle solver.

- `WeightedJacobi.m`  
  Implements the weighted Jacobi smoother with relaxation parameter omega = 2/3.

- `output.txt`  
  Contains the numerical output from running `main.m`.

- `sparsity_pattern.png`  
  Sparsity pattern of the 2D Poisson matrix.

- `residual_convergence.png`  
  Residual convergence plot for the multigrid V-cycle solver.

## How to run

Open MATLAB in this folder and run:

```matlab
main
