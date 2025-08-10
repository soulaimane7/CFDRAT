function [A, b] = get_Ab_u(t, u, v, p_pressure, inlet_func, dt, mu, rho, grid_u)
%GET_AB_U Assemble coefficient matrix and RHS for U-momentum equation.
%   Constructs the linear system Au = b for the U-velocity component using
%   hybrid upwind discretization and handles time-dependent coefficients.
%
%   Inputs:
%       t           - Current time [s]
%       u, v        - Velocity fields at current time step [m/s]
%       p_pressure  - Pressure field from previous time step [Pa]
%       inlet_func  - Inlet boundary condition function handle
%       dt, mu, rho - Time step [s], dynamic viscosity [Pa·s], density [kg/m³]
%       grid_u      - U-velocity grid structure from init_grid_u
%
%   Outputs:
%       A - Sparse coefficient matrix for momentum equation
%       b - Right-hand side vector including BC and pressure gradient terms
%
%   Discretization:
%   - Time: Implicit Euler (1st order)
%   - Convection: 1st-order upwind near boundaries, 2nd-order in interior
%   - Diffusion: 2nd-order central difference
%   - Pressure gradient: 2nd-order central difference
%
%   See also GET_AB_V, PISO_PREDICT, INTERP_V_TO_U.

%% Extract grid parameters and constants
h = grid_u.h;
Ny_u = grid_u.Ny;
Nx_u = grid_u.Nx;

% Node classification masks
is_1st_order = grid_u.is_1st_order;
is_2nd_order = grid_u.is_2nd_order;
is_inlet = grid_u.is_inlet;

% Node index mappings for stencil assembly
p_idx = grid_u.p_idx; s_idx = grid_u.s_idx; n_idx = grid_u.n_idx;
w_idx = grid_u.w_idx; e_idx = grid_u.e_idx;
ss_idx = grid_u.ss_idx; nn_idx = grid_u.nn_idx;
ww_idx = grid_u.ww_idx; ee_idx = grid_u.ee_idx;

% Pre-compute frequently used coefficients
one_over_dt = 1/dt;
mu_over_h2 = mu/h^2;
one_over_h = 1/h;
one_over_2h = 1/(2*h);

%% Velocity interpolation and upwind factors
% Interpolate V-velocity to U-node locations for cross-derivative terms
v_on_u = interp_v_to_u(v);

% Compute velocity components for upwind discretization
u_pos = max(u, 0); u_neg = min(u, 0); u_abs = abs(u);
v_on_u_pos = max(v_on_u, 0); v_on_u_neg = min(v_on_u, 0); v_on_u_abs = abs(v_on_u);

%% Finite difference coefficient matrices
% 1st-order upwind coefficients (robust near boundaries)
C1_s = -mu_over_h2 - v_on_u_pos * one_over_h;
C1_w = -mu_over_h2 - u_pos * one_over_h;
C1_p = one_over_dt + u_abs * one_over_h + v_on_u_abs * one_over_h + 4*mu_over_h2;
C1_e = -mu_over_h2 + u_neg * one_over_h;
C1_n = -mu_over_h2 + v_on_u_neg * one_over_h;

% 2nd-order upwind coefficients (accurate in interior)
C2_ss = v_on_u_pos * one_over_2h;
C2_s  = -mu_over_h2 - 2*v_on_u_pos * one_over_h;
C2_ww = u_pos * one_over_2h;
C2_w  = -mu_over_h2 - 2*u_pos * one_over_h;
C2_p  = one_over_dt + 3*u_abs * one_over_2h + 3*v_on_u_abs * one_over_2h + 4*mu_over_h2;
C2_e  = -mu_over_h2 + 2*u_neg * one_over_h;
C2_ee = -u_neg * one_over_2h;
C2_n  = -mu_over_h2 + 2*v_on_u_neg * one_over_h;
C2_nn = -v_on_u_neg * one_over_2h;

%% Sparse matrix assembly using IJV triplet format
% Estimate non-zero entries for efficient memory allocation
nnz_const = length(grid_u.I_const);
nnz_1st = 5 * nnz(is_1st_order);
nnz_2nd = 9 * nnz(is_2nd_order);
total_nnz = nnz_const + nnz_1st + nnz_2nd;

I = zeros(total_nnz, 1);
J = zeros(total_nnz, 1);
V = zeros(total_nnz, 1);

% Insert pre-computed constant boundary condition coefficients
I(1:nnz_const) = grid_u.I_const;
J(1:nnz_const) = grid_u.J_const;
V(1:nnz_const) = grid_u.V_const;
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
NyNx = Ny_u * Nx_u;
A = sparse(I, J, V, NyNx, NyNx);

%% Right-hand side vector assembly
b_matrix = zeros(Ny_u, Nx_u);

% Inlet boundary condition: time-dependent parabolic profile
yy_inlet = grid_u.yy(is_inlet);
u_inlet_values = inlet_func(t, yy_inlet);
b_matrix(is_inlet) = u_inlet_values;
% Other BC nodes default to zero

% Interior nodes: temporal term minus pressure gradient
grad_p_x = (p_pressure(:, 2:end) - p_pressure(:, 1:end-1)) / h;
u_internal = u(2:end-1, 2:end-1);
b_matrix(2:end-1, 2:end-1) = u_internal * one_over_dt - grad_p_x / rho;

% Enforce zero RHS for solid nodes
b_matrix(grid_u.is_solid) = 0;
b_matrix(grid_u.is_solid_boundary) = 0;

% Convert to column vector (MATLAB sparse convention)
b_transposed = b_matrix';
b = b_transposed(:);

end


function v_on_u = interp_v_to_u(v)
%INTERP_V_TO_U Interpolate V-velocity from V-nodes to U-node locations.
%   Uses 4-point averaging for interior nodes and extrapolates for ghost
%   cells assuming no-slip boundary conditions.

    % 4-point average for interior U-nodes
    v_on_u_interior = (v(1:end-1, 1:end-1) + v(1:end-1, 2:end) + ...
                       v(2:end, 1:end-1)   + v(2:end, 2:end)) / 4;
    
    % Ghost cell values from no-slip condition: (v_ghost + v_interior)/2 = 0
    v_on_u_ghost_bottom = -v_on_u_interior(1, :);
    v_on_u_ghost_top    = -v_on_u_interior(end, :);
    
    v_on_u = [v_on_u_ghost_bottom; v_on_u_interior; v_on_u_ghost_top];
end
