function [u_star_star, v_star_star, dp, precond_dp] = ...
         piso_correct1(u_star, v_star, grid_u, grid_v, grid_p, A_dp, dp_old, precond_dp, solver_opt)
%PISO_CORRECT1 First correction step of the PISO algorithm.
%   Computes first pressure correction and corrects intermediate velocity
%   field to better satisfy continuity equation.
%
%   Inputs:
%       u_star, v_star - Intermediate velocity fields from predictor step [m/s]
%       grid_u/v/p     - Staggered grid structures
%       A_dp           - Pressure Poisson matrix (constant)
%       dp_old         - Previous pressure correction (initial guess) [Pa]
%       precond_dp     - Pressure equation preconditioner
%       solver_opt     - Linear solver configuration
%
%   Outputs:
%       u_star_star, v_star_star - First-corrected velocity fields [m/s]
%       dp             - First pressure correction [Pa]
%       precond_dp     - Updated preconditioner (if applicable)
%
%   Algorithm Steps:
%   1. Assemble RHS: b = (ρ/Δt) * ∇·u*
%   2. Solve: ∇²p' = b
%   3. Correct: u** = u* - (Δt/ρ) * ∇p'
%
%   See also PISO_PREDICT, PISO_CORRECT2, GET_B_DP1, PRESSURE_CORRECTION.

% Assemble RHS of pressure Poisson equation from velocity divergence
b_dp = get_b_dp1(u_star, v_star, grid_p);

% Solve for first pressure correction
[dp, precond_dp] = solve_leq(A_dp, b_dp, dp_old, precond_dp, solver_opt);

% Apply pressure correction to velocity field
[u_star_star, v_star_star] = pressure_correction(u_star, v_star, dp, grid_u, grid_v);

end
