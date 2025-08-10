function grid_u = init_grid_u(block_info)
%INIT_GRID_U Initialize staggered grid structure for U-velocity field.
%   Sets up the staggered grid for the U-momentum equation, including node
%   classifications, boundary condition masks, and pre-assembled constant
%   matrix coefficients. U-velocities are located on vertical cell faces.
%
%   Input:
%       block_info - Cell array of obstacle boundary information
%
%   Output:
%       grid_u - Complete grid structure containing:
%                .is_1st_order/.is_2nd_order - Discretization scheme masks
%                .I_const/.J_const/.V_const  - Pre-assembled BC coefficients
%                .p_idx, .s_idx, etc.        - Node index mappings
%                (plus geometric parameters and boundary masks)
%
%   Grid Layout:
%   - Staggered grid: U-velocities on vertical faces (Ny+2 × Nx+1)
%   - Includes ghost cells for top/bottom boundary conditions
%   - Inlet/outlet boundaries handled via Dirichlet/Neumann conditions
%
%   See also INIT_GRID_V, INIT_GRID_P, GET_AB_U.

global Nx Ny h H L slip_opt

%% Staggered grid geometry for U-velocity
% U-velocities located on vertical faces of pressure cells
Ny_u = Ny + 2;                       % Include top/bottom ghost rows
Nx_u = Nx + 1;                       % Vertical faces per row

% Physical coordinates of U-velocity nodes
xx = h * (0:Nx);                     % Face centers: x = 0, h, 2h, ...
yy = h * (-0.5 : 1 : Ny + 0.5);      % Includes ghost cells: y = -h/2, h/2, ...
[XX, YY] = meshgrid(xx, yy);

%% Node classification based on geometry
% Primary solid/fluid classification from obstacle geometry
is_solid = get_blocked_mask(block_info, XX, YY);
is_fluid = ~is_solid;

% Immersed boundary detection (fluid nodes adjacent to solid)
is_solid_s = false(Ny_u, Nx_u); is_solid_s(2:end, :) = is_solid(1:end-1, :);
is_solid_n = false(Ny_u, Nx_u); is_solid_n(1:end-1, :) = is_solid(2:end, :);
is_solid_w = false(Ny_u, Nx_u); is_solid_w(:, 2:end) = is_solid(:, 1:end-1);
is_solid_e = false(Ny_u, Nx_u); is_solid_e(:, 1:end-1) = is_solid(:, 2:end);
is_solid_boundary = is_fluid & (is_solid_s | is_solid_n | is_solid_w | is_solid_e);

% Domain boundary identification
is_inlet        = false(Ny_u, Nx_u); is_inlet(2:end-1, 1) = true;      % Left boundary
is_outlet       = false(Ny_u, Nx_u); is_outlet(2:end-1, end) = true;   % Right boundary
is_ghost_bottom = false(Ny_u, Nx_u); is_ghost_bottom(1, :) = true;     % Bottom ghost
is_ghost_top    = false(Ny_u, Nx_u); is_ghost_top(end, :) = true;      % Top ghost

%% Equation type assignment for each node
% Dirichlet nodes: prescribed velocity values
is_dirichlet = is_solid | is_solid_boundary | is_inlet;

% Ghost nodes: wall boundary conditions (slip/no-slip)
is_ghost = is_ghost_bottom | is_ghost_top;

% Neumann nodes: zero-gradient outflow condition
is_neumann = is_outlet;

% PDE nodes: full momentum equation solved
is_pde_solve = is_fluid & ~is_dirichlet & ~is_ghost & ~is_neumann;

%% Discretization scheme selection
% Use 1st-order near boundaries, 2nd-order in interior
is_not_pde_solve = ~is_pde_solve;
is_neighbor_s = false(Ny_u, Nx_u); is_neighbor_s(2:end, :) = is_not_pde_solve(1:end-1, :);
is_neighbor_n = false(Ny_u, Nx_u); is_neighbor_n(1:end-1, :) = is_not_pde_solve(2:end, :);
is_neighbor_w = false(Ny_u, Nx_u); is_neighbor_w(:, 2:end) = is_not_pde_solve(:, 1:end-1);
is_neighbor_e = false(Ny_u, Nx_u); is_neighbor_e(:, 1:end-1) = is_not_pde_solve(:, 2:end);
is_near_boundary = is_neighbor_s | is_neighbor_n | is_neighbor_w | is_neighbor_e;

is_1st_order = is_pde_solve & is_near_boundary;   % Conservative near boundaries
is_2nd_order = is_pde_solve & ~is_near_boundary;  % Accurate in interior

%% Pre-assembly of constant matrix coefficients
% Create global node indexing for sparse matrix assembly
[MM, NN] = meshgrid(1:Nx_u, 1:Ny_u);
p_idx = (NN-1)*Nx_u + MM;                        % Current node index
s_idx = circshift(p_idx, [1, 0]);                % South neighbor
n_idx = circshift(p_idx, [-1, 0]);               % North neighbor
w_idx = circshift(p_idx, [0, 1]);                % West neighbor

% Pre-allocate triplet arrays for sparse matrix assembly
nnz_estimate = nnz(is_dirichlet) + 2*nnz(is_ghost) + 2*nnz(is_neumann);
I_const = zeros(nnz_estimate, 1);
J_const = zeros(nnz_estimate, 1);
V_const = zeros(nnz_estimate, 1);
current_pos = 0;

% Dirichlet boundary conditions: u = prescribed_value
[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_dirichlet, ...
    {p_idx}, 1, ...
    I_const, J_const, V_const, current_pos, p_idx);

% Wall boundary conditions (ghost cells)
if strcmp(slip_opt, 'no-slip')
    % No-slip: (u_ghost + u_interior)/2 = 0
    [I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_bottom, ...
        {p_idx, n_idx}, [0.5; 0.5], ...
        I_const, J_const, V_const, current_pos, p_idx);
    [I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_top, ...
        {p_idx, s_idx}, [0.5; 0.5], ...
        I_const, J_const, V_const, current_pos, p_idx);
elseif strcmp(slip_opt, 'slip')
    % Slip: u_ghost = u_interior
    [I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_bottom, ...
        {p_idx, n_idx}, [1; -1], ...
        I_const, J_const, V_const, current_pos, p_idx);
    [I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_top, ...
        {p_idx, s_idx}, [1; -1], ...
        I_const, J_const, V_const, current_pos, p_idx);
else
    error('Invalid slip_opt. Choose "slip" or "no-slip".');
end

% Neumann outflow: ∂u/∂x = 0 → u_outlet = u_neighbor
[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_outlet, ...
    {p_idx, w_idx}, [1; -1], ...
    I_const, J_const, V_const, current_pos, p_idx);

% Remove unused pre-allocated entries
I_const = I_const(1:current_pos);
J_const = J_const(1:current_pos);
V_const = V_const(1:current_pos);

%% Package complete grid structure
grid_u.h = h;
grid_u.L = L;
grid_u.H = H;
grid_u.Nx = Nx_u;
grid_u.Ny = Ny_u;
grid_u.xx = xx;
grid_u.yy = yy;

% Node classification masks
grid_u.is_solid = is_solid;
grid_u.is_solid_boundary = is_solid_boundary;
grid_u.is_inlet = is_inlet;
grid_u.is_outlet = is_outlet;
grid_u.is_ghost_bottom = is_ghost_bottom;
grid_u.is_ghost_top = is_ghost_top;

% Discretization scheme classification
grid_u.is_1st_order = is_1st_order;
grid_u.is_2nd_order = is_2nd_order;

% Pre-assembled constant matrix coefficients
grid_u.I_const = I_const;
grid_u.J_const = J_const;
grid_u.V_const = V_const;

% Node index mappings for stencil assembly
grid_u.p_idx = p_idx;
grid_u.s_idx = s_idx;
grid_u.n_idx = n_idx;
grid_u.w_idx = w_idx;
grid_u.e_idx = circshift(p_idx, [0, -1]);
grid_u.ss_idx = circshift(p_idx, [2, 0]);
grid_u.nn_idx = circshift(p_idx, [-2, 0]);
grid_u.ww_idx = circshift(p_idx, [0, 2]);
grid_u.ee_idx = circshift(p_idx, [0, -2]);

end
