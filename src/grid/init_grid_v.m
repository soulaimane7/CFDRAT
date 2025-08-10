function grid_v = init_grid_v(block_info)
%INIT_GRID_V Initialize staggered grid structure for V-velocity field.
%   Sets up the staggered grid for the V-momentum equation, including node
%   classifications, boundary condition masks, and pre-assembled constant
%   matrix coefficients. V-velocities are located on horizontal cell faces.
%
%   Input:
%       block_info - Cell array of obstacle boundary information
%
%   Output:
%       grid_v - Complete grid structure containing:
%                .is_1st_order/.is_2nd_order - Discretization scheme masks
%                .I_const/.J_const/.V_const  - Pre-assembled BC coefficients
%                .p_idx, .s_idx, etc.        - Node index mappings
%                (plus geometric parameters and boundary masks)
%
%   Grid Layout:
%   - Staggered grid: V-velocities on horizontal faces (Ny+1 × Nx+2)
%   - Includes ghost cells for left/right boundary conditions  
%   - Top/bottom walls enforce v=0 (no penetration condition)
%
%   See also INIT_GRID_U, INIT_GRID_P, GET_AB_V.

global Nx Ny h H L

%% Staggered grid geometry for V-velocity
% V-velocities located on horizontal faces of pressure cells
Ny_v = Ny + 1;                       % Horizontal faces per column
Nx_v = Nx + 2;                       % Include left/right ghost columns

% Physical coordinates of V-velocity nodes
xx = h * (-0.5 : 1 : Nx + 0.5);      % Includes ghost cells: x = -h/2, h/2, ...
yy = h * (0:Ny);                     % Face centers: y = 0, h, 2h, ...
[XX, YY] = meshgrid(xx, yy);

%% Node classification based on geometry
% Primary solid/fluid classification from obstacle geometry
is_solid = get_blocked_mask(block_info, XX, YY);
is_fluid = ~is_solid;

% Immersed boundary detection (fluid nodes adjacent to solid)
is_solid_s = false(Ny_v, Nx_v); is_solid_s(2:end, :) = is_solid(1:end-1, :);
is_solid_n = false(Ny_v, Nx_v); is_solid_n(1:end-1, :) = is_solid(2:end, :);
is_solid_w = false(Ny_v, Nx_v); is_solid_w(:, 2:end) = is_solid(:, 1:end-1);
is_solid_e = false(Ny_v, Nx_v); is_solid_e(:, 1:end-1) = is_solid(:, 2:end);
is_solid_boundary = is_fluid & (is_solid_s | is_solid_n | is_solid_w | is_solid_e);

% Domain boundary identification
is_wall_bottom  = false(Ny_v, Nx_v); is_wall_bottom(1, :) = true;       % Bottom wall (v=0)
is_wall_top     = false(Ny_v, Nx_v); is_wall_top(end, :) = true;        % Top wall (v=0)
is_ghost_left   = false(Ny_v, Nx_v); is_ghost_left(2:end-1, 1) = true;  % Left ghost
is_ghost_right  = false(Ny_v, Nx_v); is_ghost_right(2:end-1, end) = true; % Right ghost

%% Equation type assignment for each node
% Dirichlet nodes: v=0 (no penetration at walls and immersed boundaries)
is_dirichlet = is_solid | is_solid_boundary | is_wall_bottom | is_wall_top;

% Neumann ghost cells: zero-gradient at inlet/outlet
is_neumann = is_ghost_left | is_ghost_right;

% PDE nodes: full momentum equation solved
is_pde_solve = is_fluid & ~is_dirichlet & ~is_neumann;

%% Discretization scheme selection
% Use robust 1st-order near boundaries, 2nd-order in interior
is_not_pde_solve = ~is_pde_solve;
is_neighbor_s = false(Ny_v, Nx_v); is_neighbor_s(2:end, :) = is_not_pde_solve(1:end-1, :);
is_neighbor_n = false(Ny_v, Nx_v); is_neighbor_n(1:end-1, :) = is_not_pde_solve(2:end, :);
is_neighbor_w = false(Ny_v, Nx_v); is_neighbor_w(:, 2:end) = is_not_pde_solve(:, 1:end-1);
is_neighbor_e = false(Ny_v, Nx_v); is_neighbor_e(:, 1:end-1) = is_not_pde_solve(:, 2:end);
is_near_boundary = is_neighbor_s | is_neighbor_n | is_neighbor_w | is_neighbor_e;

is_1st_order = is_pde_solve & is_near_boundary;   % Conservative near boundaries
is_2nd_order = is_pde_solve & ~is_near_boundary;  % Accurate in interior

%% Pre-assembly of constant matrix coefficients
% Create global node indexing for sparse matrix assembly
[MM, NN] = meshgrid(1:Nx_v, 1:Ny_v);
p_idx = (NN-1)*Nx_v + MM;                        % Current node index
w_idx = circshift(p_idx, [0, 1]);                % West neighbor
e_idx = circshift(p_idx, [0, -1]);               % East neighbor

% Pre-allocate triplet arrays for sparse matrix assembly
nnz_estimate = nnz(is_dirichlet) + 2*nnz(is_neumann);
I_const = zeros(nnz_estimate, 1);
J_const = zeros(nnz_estimate, 1);
V_const = zeros(nnz_estimate, 1);
current_pos = 0;

% Dirichlet boundary conditions: v = 0 (no penetration)
[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_dirichlet, ...
    {p_idx}, 1, ...
    I_const, J_const, V_const, current_pos, p_idx);

% Neumann ghost cells: ∂v/∂x = 0 → v_ghost = v_interior
[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_left, ...
    {p_idx, e_idx}, [1; -1], ...
    I_const, J_const, V_const, current_pos, p_idx);

[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_right, ...
    {p_idx, w_idx}, [1; -1], ...
    I_const, J_const, V_const, current_pos, p_idx);

% Remove unused pre-allocated entries
I_const = I_const(1:current_pos);
J_const = J_const(1:current_pos);
V_const = V_const(1:current_pos);

%% Package complete grid structure
grid_v.h = h;
grid_v.L = L;
grid_v.H = H;
grid_v.Nx = Nx_v;
grid_v.Ny = Ny_v;
grid_v.xx = xx;
grid_v.yy = yy;

% Node classification masks
grid_v.is_solid = is_solid;
grid_v.is_solid_boundary = is_solid_boundary;
grid_v.is_wall_bottom = is_wall_bottom;
grid_v.is_wall_top = is_wall_top;
grid_v.is_ghost_left = is_ghost_left;
grid_v.is_ghost_right = is_ghost_right;

% Discretization scheme classification
grid_v.is_1st_order = is_1st_order;
grid_v.is_2nd_order = is_2nd_order;

% Pre-assembled constant matrix coefficients
grid_v.I_const = I_const;
grid_v.J_const = J_const;
grid_v.V_const = V_const;

% Node index mappings for stencil assembly
grid_v.p_idx = p_idx;
grid_v.s_idx = circshift(p_idx, [1, 0]);
grid_v.n_idx = circshift(p_idx, [-1, 0]);
grid_v.w_idx = w_idx;
grid_v.e_idx = e_idx;
grid_v.ss_idx = circshift(p_idx, [2, 0]);
grid_v.nn_idx = circshift(p_idx, [-2, 0]);
grid_v.ww_idx = circshift(p_idx, [0, 2]);
grid_v.ee_idx = circshift(p_idx, [0, -2]);

end
