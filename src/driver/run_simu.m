function run_simu(varargin)
%RUN_SIMU Main driver for 2D incompressible flow simulation.
%   RUN_SIMU(PARAMS) executes a complete CFD simulation using the PISO 
%   algorithm on a staggered grid. Handles complex geometries defined
%   via binary images.
%
%   RUN_SIMU() runs with built-in default parameters for quick testing.
%   RUN_SIMU(PARAMS) runs with user-specified parameters from GUI or script.
%
%   Input:
%       PARAMS - (optional) Struct containing simulation parameters:
%                .H                 - Domain height [m]
%                .u_inlet_func_str  - Inlet velocity expression string
%                .t_simu           - Total simulation time [s]
%                .T_record         - Data output interval [s]
%                .rho              - Fluid density [kg/m³]
%                .mu               - Dynamic viscosity [Pa*s]
%                .node_scale       - Approximate total grid nodes
%                .speed_opt        - CFL-based time step ('fast'|'medium')
%                .slip_opt         - Wall condition ('no-slip'|'slip')
%                .ramp_opt         - Inlet ramp type ('linear'|'smoothstep'|'none')
%                .image_filename   - Path to obstacle image file
%                .save_filename    - Output .mat file path
%                .image_display_opt- Visualization mode
%
%   Example:
%       % Quick test with defaults
%       run_simu();
%       
%       % Custom simulation
%       params.H = 0.5;
%       params.t_simu = 2.0;
%       params.image_filename = 'cylinder.png';
%       run_simu(params);
%
%   Simulation pipeline:
%   1. Initialize geometry and grids
%   2. Check Reynolds number stability
%   3. Setup linear solvers and preconditioners  
%   4. Execute time-stepping loop (PISO algorithm)
%   5. Save results and launch visualization
%
%   See also CFDRAT, SOLVE_PDE, INIT_BLOCK, DRAW_UV.

% Parse input arguments - use defaults if none provided
if nargin == 0 || isempty(varargin{1})
    % Default test case: flow past obstacle at moderate Reynolds number
    params = struct();
    params.H = 0.3;                                 % Domain height 
    params.u_inlet_func_str = '0.15*4*y*(H-y)/H^2'; % Parabolic inlet profile
    params.t_ramp = 0.5;                            % Inlet ramp-up time 
    params.slip_opt = "no-slip";                    % Wall boundary condition
    params.ramp_opt = "linear";                     % Inlet ramp type
    params.rho = 1;                                 % Fluid density 
    params.mu = 1e-4;                              % Dynamic viscosity 
    params.t_simu = 1;                             % Total simulation time 
    params.T_record = 0.05;                        % Output time interval 
    params.node_scale = 5e4;                       % Target grid resolution
    params.speed_opt = "fast";                     % CFL-based time stepping
    params.image_filename = "test_fig.png";        % Obstacle geometry file
    params.save_filename = "sim_result.mat";       % Results output file
    params.image_display_opt = "show image";       % Visualization preference
    params.enable_parallel = "off";
else
    % Use parameters from GUI or external script
    params = varargin{1};
end

% Initialize global simulation parameters
global H L Nx Ny h dt rho mu slip_opt simulation_control

enable_parallel = params.enable_parallel;         % Parallel computing flag
mu = params.mu;                                   % Extract fluid properties
rho = params.rho;
t_simu = params.t_simu;                          % Extract time parameters      
T_record = params.T_record;
node_scale = params.node_scale;                  % Extract grid parameters
H = params.H;                                    % Domain height 
u_inlet_func_str = params.u_inlet_func_str;     % Extract boundary conditions        
t_ramp = params.t_ramp;                          % Inlet velocity ramp-up time 
slip_opt = params.slip_opt;                      % Wall boundary condition type
ramp_opt = params.ramp_opt;                      % Inlet ramp function type
speed_opt = params.speed_opt;                    % Extract solver options
image_filename = params.image_filename;          % Extract I/O parameters
image_display_opt = params.image_display_opt;   % Post-processing visualization mode
save_filename = params.save_filename;            % Results output file path

% Numerical solver configuration (fine-tuned for stability)
Re_limit = 300;                                  % Maximum Reynolds number
threshold_increasment_ratio = 1.6;               % Preconditioner update trigger
threshold_maxmin_ratio = 2.5;                    % Min/max iteration ratio trigger
uplimit_recorder = 20;                           % Max preconditioner reuse count
max_iter = 200;                                  % Linear solver iteration limit
residual_tol = 5e-4;                            % Convergence tolerance

% Execute complete simulation pipeline
u_inlet_func = get_inlet_func(u_inlet_func_str, ramp_opt, t_ramp, H, slip_opt);  % Create inlet velocity function handle
[u_inlet_max, u_all_max] = get_u_max(u_inlet_func, H, t_simu, t_ramp);         % Estimate maximum velocities for CFL
[block_info, L, h, Nx, Ny] = init_block(image_filename, H, node_scale);        % Generate geometry and grid parameters
check_re(rho, u_inlet_max, mu, H, block_info, Re_limit);                       % Validate Reynolds number stability
dt = init_dt(speed_opt, u_all_max);                                            % Calculate stable time step
grid_u = init_grid_u(block_info);                                              % Initialize velocity grids
grid_v = init_grid_v(block_info);                                              % Initialize v-velocity grid and BC matrices
[A_dp, grid_p] = init_grid_p(block_info);                                      % Initialize pressure grid and matrix
continue_flag = confirm_grid(grid_p, block_info);                              % Visual grid validation
if ~continue_flag; return; end
precond = init_precond(threshold_increasment_ratio, threshold_maxmin_ratio, uplimit_recorder, A_dp);  % Setup adaptive preconditioners
solver_opt = init_solver(max_iter, residual_tol);                              % Configure iterative solver parameters
init_parpool(enable_parallel);                                                 % Setup parallel computing
[u_all, v_all, p_all, t_solve] = solve_pde(t_simu, T_record, u_inlet_func, grid_u, grid_v, grid_p, A_dp, precond, solver_opt);  % Execute PISO time-stepping loop
if strcmp(simulation_control, 'stopped')
    return;
else
    save_data(u_all, v_all, p_all, u_inlet_func, node_scale, dt, t_simu, T_record, t_solve, grid_u, grid_v, grid_p, block_info, save_filename);  % Package and save results
    draw_uv(save_filename, image_display_opt);                                     % Launch interactive visualization
end

end

