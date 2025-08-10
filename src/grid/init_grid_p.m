function [A_dp, grid_p] = init_grid_p(block_info)
%INIT_GRID_P Initialize pressure grid and assemble Poisson matrix.
%   Creates the pressure field grid structure and assembles the constant
%   coefficient matrix for the pressure Poisson equation ∇²p = RHS.
%   Handles complex geometries with comprehensive boundary condition
%   classification.
%
%   Input:
%       block_info - Cell array of obstacle boundary information
%
%   Outputs:
%       A_dp    - Sparse pressure Poisson matrix (positive-definite)
%       grid_p  - Complete pressure grid structure with node classifications
%
%   Key Features:
%   - Mutually exclusive node classification (interior, boundary, corner, etc.)
%   - 5-point finite difference stencil with boundary modifications
%   - Efficient IJV triplet format for sparse matrix assembly
%   - Handles Neumann (wall) and Dirichlet (outlet) boundary conditions
%
%   Node Types:
%   - Interior: standard 5-point Laplacian
%   - Single boundary: 4-point stencil  
%   - Corner: 3-point stencil
%   - Concave corner: 2-point stencil
%
%   See also INIT_GRID_U, INIT_GRID_V, CONFIRM_GRID.

global Nx Ny h H L

%% Pressure grid geometry (cell centers)
Ny_p = Ny;
Nx_p = Nx;
NyNx_p = Ny_p * Nx_p;                            % Total pressure nodes

% Physical coordinates at cell centers
xx = h * (0.5 : 1 : Nx-0.5);                    % x = h/2, 3h/2, ...
yy = h * (0.5 : 1 : Ny-0.5);                    % y = h/2, 3h/2, ...
[XX, YY] = meshgrid(xx, yy);

%% Primary geometric classification
% Determine solid/fluid status from obstacle geometry
is_solid = get_blocked_mask(block_info, XX, YY);
is_fluid = ~is_solid;

% Detect immersed boundaries (fluid nodes adjacent to solid)
is_solid_s = false(Ny_p, Nx_p); is_solid_s(2:end, :) = is_solid(1:end-1, :); % South neighbor is solid
is_solid_n = false(Ny_p, Nx_p); is_solid_n(1:end-1, :) = is_solid(2:end, :); % North neighbor is solid
is_solid_w = false(Ny_p, Nx_p); is_solid_w(:, 2:end) = is_solid(:, 1:end-1); % West neighbor is solid
is_solid_e = false(Ny_p, Nx_p); is_solid_e(:, 1:end-1) = is_solid(:, 2:end); % East neighbor is solid
is_solid_boundary = is_fluid & (is_solid_s | is_solid_w | is_solid_e | is_solid_n);

% Domain boundaries
is_domain_b = false(Ny_p, Nx_p); is_domain_b(1, :) = true;     % Bottom wall
is_domain_t = false(Ny_p, Nx_p); is_domain_t(end, :) = true;   % Top wall
is_domain_l = false(Ny_p, Nx_p); is_domain_l(:, 1) = true;     % Left wall (inlet)
is_domain_r = false(Ny_p, Nx_p); is_domain_r(:, end) = true;   % Right wall (outlet)

%% Boundary condition type identification
% Neumann boundaries: ∂p/∂n = 0 (walls and immersed boundaries)
is_neumann_s = (is_domain_b | (is_fluid & is_solid_s)) & ~is_solid;
is_neumann_n = (is_domain_t | (is_fluid & is_solid_n)) & ~is_solid;
is_neumann_w = (is_domain_l | (is_fluid & is_solid_w)) & ~is_solid;
is_neumann_e = (is_fluid & is_solid_e) & ~is_solid;
% Dirichlet boundary: p = 0 (outlet reference pressure)
is_dirichlet_e = is_domain_r & ~is_neumann_e & ~is_solid;

%% Mutually exclusive node classification
% Count boundary conditions per node for systematic classification
num_neumann_boundaries = double(is_neumann_s) + double(is_neumann_n) + double(is_neumann_w) + double(is_neumann_e);
num_dirichlet_boundaries = double(is_dirichlet_e);

% Interior fluid nodes (standard 5-point stencil)
is_fluid_interior = is_fluid & (num_neumann_boundaries == 0) & (num_dirichlet_boundaries == 0);

% Single boundary nodes (4-point stencil)
is_neumann_s_only = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 0) & is_neumann_s;
is_neumann_n_only = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 0) & is_neumann_n;
is_neumann_w_only = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 0) & is_neumann_w;
is_neumann_e_only = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 0) & is_neumann_e;
is_dirichlet_e_only = is_fluid & (num_dirichlet_boundaries == 1) & (num_neumann_boundaries == 0) & is_dirichlet_e;

% Corner nodes with two boundaries (3-point stencil)
is_corner_S_neumann_W_neumann = is_fluid & (num_neumann_boundaries == 2) & (num_dirichlet_boundaries == 0) & is_neumann_s & is_neumann_w;
is_corner_N_neumann_W_neumann = is_fluid & (num_neumann_boundaries == 2) & (num_dirichlet_boundaries == 0) & is_neumann_n & is_neumann_w;
is_corner_S_neumann_E_neumann = is_fluid & (num_neumann_boundaries == 2) & (num_dirichlet_boundaries == 0) & is_neumann_s & is_neumann_e;
is_corner_N_neumann_E_neumann = is_fluid & (num_neumann_boundaries == 2) & (num_dirichlet_boundaries == 0) & is_neumann_n & is_neumann_e;
% Mixed corners (Neumann + Dirichlet)
is_corner_S_neumann_E_dirichlet = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 1) & is_neumann_s & is_dirichlet_e;
is_corner_N_neumann_E_dirichlet = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 1) & is_neumann_n & is_dirichlet_e;

% Concave corner nodes with three boundaries (2-point stencil)
is_concave_NWS = is_fluid & (num_neumann_boundaries == 3) & (num_dirichlet_boundaries == 0) & is_neumann_n & is_neumann_w & is_neumann_s;
is_concave_NES = is_fluid & (num_neumann_boundaries == 3) & (num_dirichlet_boundaries == 0) & is_neumann_n & is_neumann_e & is_neumann_s;
is_concave_WNE = is_fluid & (num_neumann_boundaries == 3) & (num_dirichlet_boundaries == 0) & is_neumann_w & is_neumann_n & is_neumann_e;
is_concave_WSE = is_fluid & (num_neumann_boundaries == 3) & (num_dirichlet_boundaries == 0) & is_neumann_w & is_neumann_s & is_neumann_e;

%% Sparse matrix assembly using IJV triplet format
% Create node index mapping for efficient stencil assembly
[II, JJ] = meshgrid(1:Nx_p, 1:Ny_p);
p = (JJ-1)*Nx_p + II;                            % Current node index
s = circshift(p, [1, 0]); n = circshift(p, [-1, 0]); % South, North neighbors
w = circshift(p, [0, 1]); e = circshift(p, [0, -1]); % West, East neighbors

% Pre-allocate triplet arrays
nnz_estimate = nnz(is_solid) + 5*nnz(is_fluid);  % Conservative estimate
I = zeros(nnz_estimate, 1); J = zeros(nnz_estimate, 1); V = zeros(nnz_estimate, 1);
current_pos = 0;

% Assemble stencil coefficients for each node classification
% Solid nodes: identity equation (ensures non-singularity)
[I, J, V, current_pos] = assemble_stencil_const(is_solid, {p}, -1, I, J, V, current_pos, p);
% Interior: standard 5-point Laplacian
[I, J, V, current_pos] = assemble_stencil_const(is_fluid_interior, {s,w,p,e,n}, [1;1;-4;1;1], I, J, V, current_pos, p);
% Single boundary nodes: modified 4-point stencils
[I, J, V, current_pos] = assemble_stencil_const(is_neumann_s_only, {w,p,e,n}, [1;-3;1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_neumann_n_only, {s,w,p,e}, [1;1;-3;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_neumann_w_only, {s,p,e,n}, [1;-3;1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_neumann_e_only, {s,w,p,n}, [1;1;-3;1], I, J, V, current_pos, p);
% Corner nodes: 3-point stencils
[I, J, V, current_pos] = assemble_stencil_const(is_corner_S_neumann_W_neumann, {p,e,n}, [-2;1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_corner_N_neumann_W_neumann, {s,p,e}, [1;-2;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_corner_S_neumann_E_neumann, {w,p,n}, [1;-2;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_corner_N_neumann_E_neumann, {s,w,p}, [1;1;-2], I, J, V, current_pos, p);
% Concave corners: 2-point stencils  
[I, J, V, current_pos] = assemble_stencil_const(is_concave_NWS, {p,e}, [-1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_concave_NES, {p,w}, [-1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_concave_WNE, {p,s}, [-1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_concave_WSE, {p,n}, [-1;1], I, J, V, current_pos, p);
% Dirichlet boundary: enhanced stencil for reference pressure
[I, J, V, current_pos] = assemble_stencil_const(is_dirichlet_e_only, {s,w,p,n}, [1;1;-5;1], I, J, V, current_pos, p);
% Mixed boundary corners
[I, J, V, current_pos] = assemble_stencil_const(is_corner_S_neumann_E_dirichlet, {w,p,n}, [1;-4;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_corner_N_neumann_E_dirichlet, {s,w,p}, [1;1;-4], I, J, V, current_pos, p);

% Finalize sparse matrix
I = I(1:current_pos); J = J(1:current_pos); V = V(1:current_pos);
V = -V / h^2;                                    % Sign flip for positive-definiteness + grid scaling
A_dp = sparse(I, J, V, NyNx_p, NyNx_p);

%% Package complete grid structure
grid_p.h = h; grid_p.L = L; grid_p.H = H;
grid_p.Nx_p = Nx_p; grid_p.Ny_p = Ny_p;
grid_p.xx = xx; grid_p.yy = yy;
grid_p.is_solid = is_solid; grid_p.is_fluid = is_fluid;
grid_p.is_solid_boundary = is_solid_boundary;
grid_p.is_fluid_interior = is_fluid_interior;
grid_p.is_neumann_s_only = is_neumann_s_only;
grid_p.is_neumann_n_only = is_neumann_n_only;
grid_p.is_neumann_w_only = is_neumann_w_only;
grid_p.is_neumann_e_only = is_neumann_e_only;
grid_p.is_dirichlet_e_only = is_dirichlet_e_only;
grid_p.is_corner_S_neumann_W_neumann = is_corner_S_neumann_W_neumann;
grid_p.is_corner_N_neumann_W_neumann = is_corner_N_neumann_W_neumann;
grid_p.is_corner_S_neumann_E_neumann = is_corner_S_neumann_E_neumann;
grid_p.is_corner_N_neumann_E_neumann = is_corner_N_neumann_E_neumann;
grid_p.is_corner_S_neumann_E_dirichlet = is_corner_S_neumann_E_dirichlet;
grid_p.is_corner_N_neumann_E_dirichlet = is_corner_N_neumann_E_dirichlet;
grid_p.is_concave_NWS = is_concave_NWS;
grid_p.is_concave_NES = is_concave_NES;
grid_p.is_concave_WNE = is_concave_WNE;
grid_p.is_concave_WSE = is_concave_WSE;

end

