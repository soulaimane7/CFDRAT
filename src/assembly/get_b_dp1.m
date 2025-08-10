function b = get_b_dp1(u_star, v_star, grid_p)
%GET_B_DP1 Assemble RHS for first pressure correction equation.
%   Computes right-hand side vector for the first pressure correction
%   Poisson equation in PISO algorithm by enforcing continuity constraint
%   on intermediate velocity field.
%
%   Inputs:
%       u_star, v_star - Intermediate velocity fields from predictor step [m/s]
%       grid_p         - Pressure grid structure
%
%   Output:
%       b - RHS vector for pressure correction: ∇²p' = b
%
%   Mathematical Background:
%   The intermediate velocity u* from momentum predictor generally violates
%   continuity (∇·u* ≠ 0). The pressure correction p' is computed to project
%   u* onto the divergence-free space:
%   
%   ∇²p' = (ρ/Δt)∇·u*
%   
%   Discretization on staggered grid:
%   - U-velocities on vertical faces, V-velocities on horizontal faces
%   - Divergence computed at cell centers (pressure nodes)
%   - 2nd-order central differences for spatial derivatives
%
%   See also GET_B_DP2, PISO_CORRECT1, PRESSURE_CORRECTION.

global rho dt h

%% Compute velocity divergence at pressure nodes
% Divergence calculation: ∇·u* = ∂u*/∂x + ∂v*/∂y
% On staggered grid: div = (u_east - u_west)/h + (v_north - v_south)/h

% X-direction velocity gradient (∂u*/∂x)
grad_u_star_x = (u_star(2:end-1, 2:end) - u_star(2:end-1, 1:end-1)) / h;

% Y-direction velocity gradient (∂v*/∂y)  
grad_v_star_y = (v_star(2:end, 2:end-1) - v_star(1:end-1, 2:end-1)) / h;

% Total divergence at each pressure node
divergence = grad_u_star_x + grad_v_star_y;

%% Assemble RHS vector
% Scale divergence by (ρ/Δt) for pressure correction equation
b_matrix = (rho / dt) * divergence;

% Zero RHS for solid nodes (no pressure correction needed)
b_matrix(grid_p.is_solid) = 0;

% Sign correction for negative Laplacian matrix convention
b_matrix = -b_matrix;

% Convert to column vector for linear solver
b_transposed = b_matrix';
b = b_transposed(:);

end
