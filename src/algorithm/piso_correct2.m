function [u_new, v_new, dp_corr, precond_dp] =...
         piso_correct2(u_star, v_star, u_star_star, v_star_star, grid_u, grid_v, grid_p, A_dp, dp_corr_old, precond_dp, solver_opt)
%PISO_CORRECT2 Second correction step of the PISO algorithm.
%   Computes second pressure correction and produces final velocity field
%   for the current time step.
%
%   Inputs:
%       u_star, v_star         - Intermediate velocity from predictor [m/s]
%       u_star_star, v_star_star - First-corrected velocity fields [m/s]
%       grid_u/v/p             - Staggered grid structures
%       A_dp                   - Pressure Poisson matrix (constant)
%       dp_corr_old           - Previous second correction (initial guess) [Pa]
%       precond_dp            - Pressure equation preconditioner
%       solver_opt            - Linear solver configuration
%
%   Outputs:
%       u_new, v_new  - Final corrected velocity fields for time step [m/s]
%       dp_corr       - Second pressure correction [Pa]
%       precond_dp    - Updated preconditioner (if applicable)
%
%   Algorithm Steps:
%   1. Assemble RHS: b ≈ ρ * ∇·(L(u**) - L(u*))
%   2. Solve: ∇²p'' = b  
%   3. Correct: u_new = u** - (Δt/ρ) * ∇p''
%
%   See also PISO_CORRECT1, GET_B_DP2, PRESSURE_CORRECTION.

% Assemble RHS based on difference between corrected and intermediate velocities
b_dp_corr = get_b_dp2(u_star, v_star, u_star_star, v_star_star, grid_p);

% Solve for second pressure correction
[dp_corr, precond_dp] = solve_leq(A_dp, b_dp_corr, dp_corr_old, precond_dp, solver_opt);

% Apply second correction to obtain final velocity field
[u_new, v_new] = pressure_correction(u_star_star, v_star_star, dp_corr, grid_u, grid_v);

end
