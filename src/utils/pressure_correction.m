function [u_new, v_new] = pressure_correction(u_old, v_old, dp, grid_u, grid_v)
%PRESSURE_CORRECTION Apply pressure gradient correction to velocity field.
%   Implements the core step of projection method by correcting velocity
%   field with pressure gradient to enforce mass conservation.
%
%   Inputs:
%       u_old, v_old - Uncorrected velocity fields 
%       dp           - Pressure correction field 
%       grid_u/v     - Velocity grid structures
%
%   Outputs:
%       u_new, v_new - Pressure-corrected velocity fields
%
%   Mathematical Operation:
%   u_new = u_old - (Δt/ρ)∇dp
%   
%   The correction projects the velocity field onto the divergence-free
%   space while maintaining boundary condition compliance.
%
%   Algorithm:
%   1. Compute pressure gradients on staggered grid locations
%   2. Apply correction to interior fluid nodes only
%   3. Re-enforce all physical boundary conditions
%
%   See also PISO_CORRECT1, PISO_CORRECT2, GET_B_DP1, GET_B_DP2.

global rho dt h slip_opt

%% Compute pressure gradients at velocity locations
% X-gradient at U-velocity nodes (vertical faces)
grad_dp_x = (dp(:, 2:end) - dp(:, 1:end-1)) / h;
% Y-gradient at V-velocity nodes (horizontal faces) 
grad_dp_y = (dp(2:end, :) - dp(1:end-1, :)) / h;

%% Apply pressure correction to interior nodes
% Initialize corrected fields
u_new = u_old;
v_new = v_old;

% Correction formula: u_new = u_old - (dt/rho)*grad_p
u_new(2:end-1, 2:end-1) = u_old(2:end-1, 2:end-1) - (dt / rho) * grad_dp_x;
v_new(2:end-1, 2:end-1) = v_old(2:end-1, 2:end-1) - (dt / rho) * grad_dp_y;

%% Re-enforce boundary conditions
% Pressure correction only updates interior nodes; boundary conditions
% must be re-applied to maintain physical consistency

% Solid boundaries: zero velocity
u_new(grid_u.is_solid) = 0;
v_new(grid_v.is_solid) = 0;
u_new(grid_u.is_solid_boundary) = 0;
v_new(grid_v.is_solid_boundary) = 0;

% Wall boundary conditions for U-velocity (top/bottom)
if strcmp(slip_opt, 'no-slip')
   % No-slip: u_ghost = -u_interior (enforces u=0 at wall)
    u_new([1, end], :) = -u_new([2, end-1], :);
elseif strcmp(slip_opt, 'slip')
    % Slip: u_ghost = u_interior (enforces du/dy=0 at wall)
    u_new([1, end], :) = u_new([2, end-1], :);
else
    error('invalid slip_opt!');
end

% V-velocity boundary conditions
% Left boundary (inlet): no-penetration v=0
v_new(2:end-1, 1) = -v_new(2:end-1, 2);

% Right boundary (outlet): zero normal gradient
v_new(2:end-1, end) = v_new(2:end-1, end-1);

end
