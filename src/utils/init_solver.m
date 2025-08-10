% configure_solver.m
function solver_opt = init_solver(max_iter, residual_tol)
%INIT_SOLVER Configure iterative linear solver parameters.
%   Sets convergence criteria for GMRES-based linear system solvers
%   used in the PISO algorithm.
%
%   Inputs:
%       max_iter     - Maximum solver iterations before failure
%       residual_tol - Relative residual tolerance for convergence
%
%   Output:
%       solver_opt - Solver configuration structure
%
%   Usage:
%       solver_opt = init_solver(1000, 1e-6);
%
%   See also SOLVE_LINEAR_SYSTEM, GMRES.

solver_opt.max_iter = max_iter;                  % Iteration limit
solver_opt.tol = residual_tol;                   % Convergence tolerance
end
