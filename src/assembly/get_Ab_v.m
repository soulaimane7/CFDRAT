function [A, b] = get_Ab_v(u, v, p_pressure, dt, mu, rho, grid_v)
%GET_AB_V Assemble coefficient matrix and RHS for V-momentum equation.
%   Constructs the linear system Av = b for the V-velocity component using
%   hybrid upwind discretization and handles time-dependent coefficients.
%
%   Inputs:
%       u, v        - Velocity fields at current time step [m/s]
%       p_pressure  - Pressure field from previous time step [Pa]
%       dt, mu, rho - Time step [s], dynamic viscosity [Pa·s], density [kg/m³]
%       grid_v      - V-velocity grid structure from init_grid_v
%
%   Outputs:
%       A - Sparse coefficient matrix for momentum equation
%       b - Right-hand side vector including pressure gradient terms
%
%   Discretization:
%   - Time: Implicit Euler (1st order)
%   - Convection: 1st-order upwind near boundaries, 2nd-order in interior
%   - Diffusion: 2nd-order central difference
%   - Pressure gradient: 2nd-order central difference (y-direction)
%
%   Note: V-velocity has no-penetration BC (v=0) at top/bottom walls
%   and zero-gradient BC at inlet/outlet.
%
%   See also GET_AB_U, PISO_PREDICT, INTERP_U_TO_V.

%% Extract grid parameters and constants
h = grid_v.h;
Ny_v = grid_v.Ny;
Nx_v = grid_v.Nx;

% Node classification masks
is_1st_order = grid_v.is_1st_order;
is_2nd_order = grid_v.is_2nd_order;

% Node index mappings for stencil assembly
p_idx = grid_v.p_idx; s_idx = grid_v.s_idx; n_idx = grid_v.n_idx;
w_idx = grid_v.w_idx; e_idx = grid_v.e_idx;
ss_idx = grid_v.ss_idx; nn_idx = grid_v.nn_idx;
ww_idx = grid_v.ww_idx; ee_idx = grid_v.ee_idx;

% Pre-compute frequently used coefficients
one_over_dt = 1/dt;
mu_over_h2 = mu/h^2;
one_over_h = 1/h;
one_over_2h = 1/(2*h);

%% Velocity interpolation and upwind factors
% Interpolate U-velocity to V-node locations for cross-derivative terms
u_on_v = interp_u_to_v(u);

% Compute velocity components for upwind discretization
u_on_v_pos = max(u_on_v, 0); u_on_v_neg = min(u_on_v, 0); u_on_v_abs = abs(u_on_v);
v_pos = max(v, 0); v_neg = min(v, 0); v_abs = abs(v);

%% Finite difference coefficient matrices
% 1st-order upwind coefficients (robust near boundaries)
C1_s = -mu_over_h2 - v_pos * one_over_h;
C1_w = -mu_over_h2 - u_on_v_pos * one_over_h;
C1_p = one_over_dt + u_on_v_abs * one_over_h + v_abs * one_over_h + 4*mu_over_h2;
C1_e = -mu_over_h2 + u_on_v_neg * one_over_h;
C1_n = -mu_over_h2 + v_neg * one_over_h;

% 2nd-order upwind coefficients (accurate in interior)
C2_ss = v_pos * one_over_2h;
C2_s  = -mu_over_h2 - 2*v_pos * one_over_h;
C2_ww = u_on_v_pos * one_over_2h;
C2_w  = -mu_over_h2 - 2*u_on_v_pos * one_over_h;
C2_p  = one_over_dt + 3*u_on_v_abs * one_over_2h + 3*v_abs * one_over_2h + 4*mu_over_h2;
C2_e  = -mu_over_h2 + 2*u_on_v_neg * one_over_h;
C2_ee = -u_on_v_neg * one_over_2h;
C2_n  = -mu_over_h2 + 2*v_neg * one_over_h;
C2_nn = -v_neg * one_over_2h;

%% Sparse matrix assembly using IJV triplet format
% Estimate non-zero entries for efficient memory allocation
nnz_const = length(grid_v.I_const);
nnz_1st = 5 * nnz(is_1st_order);
nnz_2nd = 9 * nnz(is_2nd_order);
total_nnz = nnz_const + nnz_1st + nnz_2nd;

I = zeros(total_nnz, 1);
J = zeros(total_nnz, 1);
V = zeros(total_nnz, 1);

% Insert pre-computed constant boundary condition coefficients
I(1:nnz_const) = grid_v.I_const;
J(1:nnz_const) = grid_v.J_const;
V(1:nnz_const) = grid_v.V_const;
current_pos = nnz_const;

% Assemble 1st-order stencils for boundary-adjacent nodes
[I, J, V, current_pos] = assemble_stencil_dynamic(is_1st_order, ...
    {s_idx, w_idx, p_idx, e_idx, n_idx}, ...
    {C1_s, C1_w, C1_p, C1_e, C1_n}, ...
    I, J, V, current_pos, p_idx);

% Assemble 2nd-order stencils for interior nodes
[I, J, V, current_pos] = assemble_stencil_dynamic(is_2nd_order, ...
    {ss_idx, s_idx, ww_idx, w_idx, p_idx, e_idx, ee_idx, n_idx, nn_idx}, ...
    {C2_ss, C2_s, C2_ww, C2_w, C2_p, C2_e, C2_ee, C2_n, C2_nn}, ...
    I, J, V, current_pos, p_idx);

% Create final sparse matrix
NyNx = Ny_v * Nx_v;
A = sparse(I, J, V, NyNx, NyNx);

%% Right-hand side vector assembly
b_matrix = zeros(Ny_v, Nx_v);
% All boundary conditions for V-velocity default to zero

% Interior nodes: temporal term minus pressure gradient (y-direction)
grad_p_y = (p_pressure(2:end, :) - p_pressure(1:end-1, :)) / h;
v_internal = v(2:end-1, 2:end-1);
b_matrix(2:end-1, 2:end-1) = v_internal * one_over_dt - grad_p_y / rho;

% Enforce zero RHS for solid nodes
b_matrix(grid_v.is_solid) = 0;
b_matrix(grid_v.is_solid_boundary) = 0;

% Convert to column vector (MATLAB sparse convention)
b_transposed = b_matrix';
b = b_transposed(:);

end

function u_on_v = interp_u_to_v(u)
%INTERP_U_TO_V Interpolate U-velocity from U-nodes to V-node locations.
%   Uses 4-point averaging for interior nodes with appropriate boundary
%   treatments for inlet (Dirichlet) and outlet (Neumann) conditions.

    % 4-point average for interior V-nodes
    u_on_v_interior = (u(1:end-1, 1:end-1) + u(1:end-1, 2:end) + ...
                       u(2:end, 1:end-1)   + u(2:end, 2:end)) / 4;
    
    % Outlet: zero-gradient condition (Neumann)
    u_on_v_ghost_right = u_on_v_interior(:, end);

    % Inlet: extrapolation from Dirichlet condition
    u_left_bound = (u(2:end,1) + u(1:end-1,1)) / 2; 
    u_on_v_ghost_left  =  2 * u_left_bound - u_on_v_interior(:, 1);

    u_on_v = [u_on_v_ghost_left, u_on_v_interior, u_on_v_ghost_right];
end
