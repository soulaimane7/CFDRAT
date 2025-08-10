function b = get_b_dp2(u_star, v_star, u_star_star, v_star_star, grid_p)
%GET_B_DP2 Assemble RHS for second pressure correction equation (PISO).
%   Computes right-hand side for the second pressure correction in PISO
%   algorithm to better satisfy momentum equation residuals.
%
%   Inputs:
%       u_star, v_star           - Intermediate velocities from predictor [m/s]
%       u_star_star, v_star_star - First-corrected velocities [m/s]  
%       grid_p                   - Pressure grid structure
%
%   Output:
%       b - RHS vector for second correction: ∇²p'' = b
%
%   Mathematical Background:
%   The second correction addresses momentum equation residuals by computing:
%   ∇²p'' ≈ ρ ∇·[L(u**) - L(u*)]
%   where L(u) = -convection + diffusion is the spatial momentum operator.
%   This correction improves pressure-velocity coupling beyond first correction.
%
%   See also GET_B_DP1, PISO_CORRECT2, GET_CONV, GET_DIFF, GET_DIV.

global rho  

%% Compute spatial momentum operators for both velocity fields
% L(u) = -convection + diffusion terms from momentum equations
[conv_u_s, conv_v_s] = get_conv(u_star, v_star);
[diff_u_s, diff_v_s] = get_diff(u_star, v_star);
L_u_star = -conv_u_s + diff_u_s;
L_v_star = -conv_v_s + diff_v_s;

[conv_u_ss, conv_v_ss] = get_conv(u_star_star, v_star_star);
[diff_u_ss, diff_v_ss] = get_diff(u_star_star, v_star_star);
L_u_star_star = -conv_u_ss + diff_u_ss;
L_v_star_star = -conv_v_ss + diff_v_ss;

%% Calculate divergence of momentum operator difference
% ∇·[L(u**) - L(u*)] represents momentum residual spatial distribution
delta_L_u = L_u_star_star - L_u_star;
delta_L_v = L_v_star_star - L_v_star;

div_delta_L = get_div(delta_L_u, delta_L_v);

%% Assemble RHS vector
b_matrix = rho * div_delta_L;

% Zero RHS for solid nodes
b_matrix(grid_p.is_solid) = 0;

% Sign correction for negative Laplacian matrix convention
b_matrix = -b_matrix;

% Convert to column vector for linear solver
b_transposed = b_matrix';
b = b_transposed(:);

end

function [conv_u, conv_v] = get_conv(u, v)
%GET_CONV Compute convection terms at pressure nodes.
%   Calculates -(u·∇)u and -(u·∇)v using central differences.

global  h
    % Interpolate velocities to pressure node locations (cell centers)
    u_on_p = (u(2:end-1, 2:end) + u(2:end-1, 1:end-1)) / 2;
    v_on_p = (v(2:end, 2:end-1) + v(1:end-1, 2:end-1)) / 2;

    % Velocity gradients at pressure nodes (2nd-order central differences)
    grad_u_x_on_p = (u(2:end-1, 2:end) - u(2:end-1, 1:end-1)) / h;
    grad_u_y = (u(3:end, :) - u(1:end-2, :)) / (2*h);
    grad_u_y_on_p = (grad_u_y(:, 1:end-1) + grad_u_y(:, 2:end)) / 2;

    grad_v_y_on_p = (v(2:end, 2:end-1) - v(1:end-1, 2:end-1)) / h;
    grad_v_x = (v(:, 3:end) - v(:, 1:end-2)) / (2*h);
    grad_v_x_on_p = (grad_v_x(1:end-1, :) + grad_v_x(2:end, :)) / 2;

    % Convection terms: u·∇u and u·∇v
    conv_u = u_on_p .* grad_u_x_on_p + v_on_p .* grad_u_y_on_p;
    conv_v = u_on_p .* grad_v_x_on_p + v_on_p .* grad_v_y_on_p;
end


function [diff_u, diff_v] = get_diff(u, v)
%GET_DIFF Compute viscous diffusion terms at pressure nodes.
%   Calculates μ∇²u and μ∇²v using 5-point Laplacian stencil.

global mu h slip_opt
    % Interpolate velocities to pressure node locations
    u_on_p = (u(2:end-1, 2:end) + u(2:end-1, 1:end-1)) / 2;
    v_on_p = (v(2:end, 2:end-1) + v(1:end-1, 2:end-1)) / 2;

    % Build padded U-velocity matrix with appropriate boundary conditions
    u_inlet_values = u(2:end-1, 1);
    u_p_ghost_left = 2 * u_inlet_values - u_on_p(:, 1);  % Inlet (Dirichlet)
    u_p_ghost_right = u_on_p(:, end);                    % Outlet (Neumann)
    u_on_p_padded_lr = [u_p_ghost_left, u_on_p, u_p_ghost_right];

    % Wall boundary conditions for U-velocity
    if strcmp(slip_opt, 'no-slip')
        u_p_ghost_bottom = -u_on_p_padded_lr(1, :);      % No-slip: u_ghost = -u_interior
        u_p_ghost_top = -u_on_p_padded_lr(end, :);
    elseif strcmp(slip_opt, 'slip')
        u_p_ghost_bottom = u_on_p_padded_lr(1, :);       % Slip: u_ghost = u_interior
        u_p_ghost_top = u_on_p_padded_lr(end, :);
    else
        error('invalid slip option! ');
    end

    u_on_p_padded = [u_p_ghost_bottom; u_on_p_padded_lr; u_p_ghost_top];

    % Build padded V-velocity matrix with boundary conditions
    % Walls: no-penetration condition (v=0 → v_ghost = -v_interior)
    v_p_ghost_bottom = -v_on_p(1, :);
    v_p_ghost_top = -v_on_p(end, :);
    v_on_p_padded_tb = [v_p_ghost_bottom; v_on_p; v_p_ghost_top];

    % Inlet/outlet conditions for V-velocity
    v_p_ghost_left = -v_on_p_padded_tb(:, 1);           % Inlet (Dirichlet v=0)
    v_p_ghost_right = v_on_p_padded_tb(:, end);         % Outlet (Neumann)
    
    v_on_p_padded = [v_p_ghost_left, v_on_p_padded_tb, v_p_ghost_right];

    % Apply 5-point Laplacian stencil for diffusion terms
    laplacian_u = (u_on_p_padded(2:end-1, 1:end-2) + u_on_p_padded(2:end-1, 3:end) + ...
                   u_on_p_padded(1:end-2, 2:end-1) + u_on_p_padded(3:end, 2:end-1) - ...
                   4 * u_on_p_padded(2:end-1, 2:end-1)) / h^2;
               
    laplacian_v = (v_on_p_padded(2:end-1, 1:end-2) + v_on_p_padded(2:end-1, 3:end) + ...
                   v_on_p_padded(1:end-2, 2:end-1) + v_on_p_padded(3:end, 2:end-1) - ...
                   4 * v_on_p_padded(2:end-1, 2:end-1)) / h^2;

    diff_u = mu * laplacian_u;
    diff_v = mu * laplacian_v;
end


function div_L = get_div(L_u, L_v)
%GET_DIV Compute divergence of vector field at pressure nodes.
%   Calculates ∇·L using 2nd-order differences with proper boundary treatment.

    global h 
    [Ny, Nx] = size(L_u);
    grad_L_u_x = zeros(Ny, Nx);
    grad_L_v_y = zeros(Ny, Nx);

    % Interior nodes: 2nd-order central differences
    grad_L_u_x(:, 2:end-1) = (L_u(:, 3:end) - L_u(:, 1:end-2)) / (2 * h);
    grad_L_v_y(2:end-1, :) = (L_v(3:end, :) - L_v(1:end-2, :)) / (2 * h);
 
    % Boundary nodes: 2nd-order one-sided differences
    % Assumes zero normal gradient of divergence at boundaries
    
    % Left/right boundaries (x-direction)
    grad_L_u_x(:, 1)   = (-3*L_u(:, 1) + 4*L_u(:, 2) - L_u(:, 3)) / (2 * h);
    grad_L_u_x(:, end) = (3*L_u(:, end) - 4*L_u(:, end-1) + L_u(:, end-2)) / (2 * h);
    
    % Bottom/top boundaries (y-direction)
    grad_L_v_y(1, :)   = (-3*L_v(1, :) + 4*L_v(2, :) - L_v(3, :)) / (2 * h);
    grad_L_v_y(end, :) = (3*L_v(end, :) - 4*L_v(end-1, :) + L_v(end-2, :)) / (2 * h);

    div_L = grad_L_u_x + grad_L_v_y;
end
