function [u_all, v_all, p_all, t_solve] = solve_pde(t_simu, T_record, u_inlet_func, grid_u, grid_v, grid_p, A_dp, precond, solver_opt)
%SOLVE_PDE Time-stepping loop for 2D incompressible Navier-Stokes equations.
%   Executes the main CFD simulation using the PISO algorithm on a staggered
%   grid. Handles momentum prediction and two-stage pressure correction for
%   accurate pressure-velocity coupling.
%
%   Inputs:
%       t_simu        - Total simulation time
%       T_record      - Data recording interval  
%       u_inlet_func  - Inlet velocity function handle @(t,y)
%       grid_u        - U-velocity grid structure with masks and indices
%       grid_v        - V-velocity grid structure with masks and indices
%       grid_p        - Pressure grid structure with boundary classifications
%       A_dp          - Pressure Poisson matrix (sparse, constant)
%       precond       - Preconditioner configuration struct
%       solver_opt    - Linear solver parameters (tolerance, max iterations)
%
%   Outputs:
%       u_all, v_all  - Velocity field time series [m/s]
%       p_all         - Pressure field time series [Pa]  
%       t_solve       - Total computation time [s]
%
%   Algorithm (PISO):
%   1. Momentum prediction: solve momentum equations explicitly
%   2. First pressure correction: enforce continuity equation  
%   3. Second pressure correction: improve momentum equation accuracy
%   4. Update pressure field and record data
%
%   Notes:
%   - Uses global variables (dt, Nx, Ny) for grid parameters
%   - Adaptive ILU preconditioning for momentum equations
%   - Progress updates printed to command window
%
%   See also PISO_PREDICT, PISO_CORRECT1, PISO_CORRECT2.

% Access global grid and time step parameters
global dt Nx Ny simulation_control

% Extract preconditioner configurations for each field
precond_u = precond.precond_u;                    % U-momentum preconditioner
precond_v = precond.precond_v;                    % V-momentum preconditioner  
precond_dp = precond.precond_dp;                  % Pressure correction preconditioner

% Calculate simulation timing parameters
total_steps = ceil(t_simu / dt);                  % Total time steps to execute
record_interval_steps = max(1, round(T_record / dt));  % Steps between data recording
num_records = floor(total_steps / record_interval_steps) + 1;  % Total data frames

% Initialize field variables on staggered grid
u = zeros(Ny+2, Nx+1);                           % U-velocity (vertical faces, includes ghost)
v = zeros(Ny+1, Nx+2);                           % V-velocity (horizontal faces, includes ghost)
p = zeros(Ny, Nx);                               % Pressure (cell centers)
dp1 = zeros(Ny, Nx);                             % First pressure correction
dp2 = zeros(Ny, Nx);                             % Second pressure correction

% Pre-allocate storage arrays for time series data
u_all = zeros(Ny+2, Nx+1, num_records);
v_all = zeros(Ny+1, Nx+2, num_records);
p_all = zeros(Ny, Nx, num_records);

% Store initial conditions (quiescent flow)
record_count = 1;
u_all(:,:,record_count) = u;
v_all(:,:,record_count) = v;
p_all(:,:,record_count) = p;
record_count = record_count + 1;

fprintf('--- Simulation Started ---\n');
fprintf('Total time: %.2fs | Time step (dt): %.4fs | Total steps: %d\n', t_simu, dt, total_steps);

% Main PISO time-stepping loop
tic; 
for step = 1:total_steps
    pause(0.001);
    if strcmp(simulation_control, 'stopped')
        break;  % simulation stopped, jump out 
    end

    t = step * dt;                               

    % PISO Step 1: Momentum prediction (solve momentum without pressure gradient)
    [u_star1, v_star1, precond_u, precond_v] = piso_predict(t, u, v, p, grid_u, grid_v, u_inlet_func, precond_u, precond_v, solver_opt);

    % PISO Step 2: First pressure correction (enforce continuity)
    [u_star2, v_star2, dp1, precond_dp] = piso_correct1(u_star1, v_star1, grid_u, grid_v, grid_p, A_dp, dp1, precond_dp, solver_opt);

    % PISO Step 3: Second pressure correction (improve momentum accuracy)
    [u, v, dp2, precond_dp] = piso_correct2(u_star1, v_star1, u_star2, v_star2, grid_u, grid_v, grid_p, A_dp, dp2, precond_dp, solver_opt);

    % Update pressure field with corrections
    p = p + dp1 + dp2;

    % Record simulation data at specified intervals
    if mod(step, record_interval_steps) == 0 && record_count <= num_records
        u_all(:,:, record_count) = u;
        v_all(:,:, record_count) = v;
        p_all(:,:, record_count) = p;
        record_count = record_count + 1;

         % Display progress to user
        fprintf('Time: %.3fs (%.1f%%), Step: %d/%d, Elapsed: %.2f s\n', ...
            t, 100*t/t_simu, step, total_steps, toc);
    end
end
t_solve = toc;                                   % Record total computation time

if strcmp(simulation_control, 'stopped') && step < total_steps
    fprintf('--- Simulation stopped by user, completed %d/%d steps in %.2f seconds ---\n', step, total_steps, t_solve);
else
    fprintf('--- Simulation finished in %.2f seconds ---\n', t_solve);
end

end

