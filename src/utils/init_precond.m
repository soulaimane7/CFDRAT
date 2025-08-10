function precond = init_precond(threshold_increasment_ratio, threshold_maxmin_ratio, uplimit_recorder, A_dp)
%INIT_PRECOND Initialize preconditioner configurations for linear solvers.
%   Sets up preconditioner structures for momentum and pressure equations.
%   Uses adaptive ILU preconditioning for momentum equations and static
%   incomplete Cholesky for the pressure Poisson equation.
%
%   Inputs:
%       threshold_increasment_ratio - Solver iteration increase threshold for updates
%       threshold_maxmin_ratio      - Max/min iteration ratio threshold for updates  
%       uplimit_recorder           - Maximum steps before forced update
%       A_dp                       - Pressure Poisson matrix (symmetric positive-definite)
%
%   Output:
%       precond - Structure containing preconditioner configurations:
%                 .precond_u  - Adaptive ILU for U-momentum equation
%                 .precond_v  - Adaptive ILU for V-momentum equation  
%                 .precond_dp - Static ICT for pressure equation
%
%   Preconditioning Strategy:
%   - Momentum: Adaptive ILU with performance-based updates
%   - Pressure: Static incomplete Cholesky (computed once)
%   - Drop tolerances scaled with grid size for robustness
%
%   See also SOLVE_LINEAR_SYSTEM, GET_AB_U, GET_AB_V.

global Ny Nx 

precond = struct();

%% Adaptive ILU preconditioners for momentum equations
% ILU factors are computed and updated dynamically during simulation
% based on solver performance metrics

% Scale drop tolerance with grid size for optimal performance
grid_size = Ny * Nx;
if grid_size >= 1e5
    droptol = 5e-6;                              % Strict tolerance for large grids
elseif grid_size >= 1e4 
    droptol = 5e-5;                              % Medium tolerance
else
    droptol = 5e-4;                              % Relaxed tolerance for small grids
end

% Common ILU configuration for both velocity components
setup_uv.type = 'ilutp';                         % ILU with threshold and pivoting
setup_uv.droptol = droptol;                      % Drop tolerance for factorization
setup_uv.udiag = true;                           % Replace zero diagonal elements

% U-momentum preconditioner configuration
precond_u.type = 'u';
precond_u.setup = setup_uv;
precond_u.threshold_increasment_ratio = threshold_increasment_ratio;
precond_u.threshold_maxmin_ratio = threshold_maxmin_ratio;
precond_u.failure_trigger = 200;                 % Iteration limit before forced update
precond_u.uplimit_recorder = uplimit_recorder;
precond_u.L = [];                                % ILU factors computed during simulation
precond_u.U = [];
precond_u.recorder = [];                         % Performance history

% V-momentum preconditioner (identical configuration)
precond_v.type = 'v';
precond_v.setup = setup_uv;
precond_v.threshold_increasment_ratio = threshold_increasment_ratio;
precond_v.threshold_maxmin_ratio = threshold_maxmin_ratio;
precond_v.failure_trigger = 200;
precond_v.uplimit_recorder = uplimit_recorder;
precond_v.L = [];
precond_v.U = [];
precond_v.recorder = [];

%% Static incomplete Cholesky preconditioner for pressure equation
% Pre-compute factorization once (pressure matrix is constant)
setup_p.type = "ict";
setup_p.droptol = 1e-8;                          % High accuracy for pressure solve
setup_p.diagcomp = 1e-6;                         % Diagonal regularization for stability

% Pre-compute Cholesky factorization: A_dp ≈ L*L'
L = ichol(A_dp, setup_p);

% Pressure preconditioner structure
precond_dp.L = L;                                % Lower triangular factor
precond_dp.U = L';                               % Upper triangular (transpose)
precond_dp.type = 'dp';

precond.precond_u = precond_u;
precond.precond_v = precond_v;
precond.precond_dp = precond_dp;
end
