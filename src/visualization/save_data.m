function save_filename = save_data(u_all, v_all, p_all, u_inlet_func, node_scale, dt, t_simu, T_record, t_solve, grid_u, grid_v, grid_p, block_info, save_filename)
%SAVE_DATA Process and save CFD simulation results to MAT file.
%   Post-processes simulation data by interpolating staggered velocity fields
%   to consistent grid locations, optimizes storage format, and archives
%   complete simulation state with metadata.
%
%   Inputs:
%       u_all, v_all, p_all - Time series of velocity and pressure fields
%       u_inlet_func        - Inlet boundary condition function handle
%       node_scale, dt      - Grid resolution and time step parameters
%       t_simu, T_record    - Simulation and recording time intervals
%       t_solve            - Total computation time [s]
%       grid_u/v/p         - Grid structures for staggered mesh
%       block_info         - Domain geometry information
%       save_filename      - Output file path
%
%   Output:
%       save_filename - Confirmed output file path
%
%   Processing:
%   - Interpolates U/V velocities from face centers to consistent locations
%   - Interpolates pressure from cell centers to grid nodes
%   - Converts to single precision for storage efficiency
%   - Bundles simulation parameters and metadata
%   - Saves in compressed MAT format for large dataset compatibility
%
%   See also LOAD, INTERP_STAGGERED_TO_NODES.

%% Ensure output directory exists
[save_path, ~, ~] = fileparts(save_filename);
if ~isempty(save_path) && ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('Created output directory: %s\n', save_filename);
end

%% Initialize data structure and allocate memory
info = struct();

[Ny_u, Nx_u, num_frames] = size(u_all);
[Ny_v, Nx_v, ~] = size(v_all);
[Ny_p, Nx_p, ~] = size(p_all);

% Pre-allocate with optimal data type for storage
info.u_all = zeros(Ny_u - 1, Nx_u, num_frames, 'single');
info.v_all = zeros(Ny_v, Nx_v - 1, num_frames, 'single');
info.p_all = zeros(Ny_p + 1, Nx_p + 1, num_frames, 'single');  % Pressure at nodes

%% Process velocity and pressure fields with interpolation
for k = 1:num_frames
    % U-velocity: interpolate from vertical faces to corner nodes
    u_slice = single(u_all(:,:,k));
    info.u_all(:,:,k) = (u_slice(1:end-1,:) + u_slice(2:end,:)) / 2;
    
    % V-velocity: interpolate from horizontal faces to corner nodes  
    v_slice = single(v_all(:,:,k));
    info.v_all(:,:,k) = (v_slice(:,1:end-1) + v_slice(:,2:end)) / 2;
    
    % Pressure: interpolate from cell centers to grid nodes
    p_slice = single(p_all(:,:,k));
    p_nodes = interp_p(p_slice);
    info.p_all(:,:,k) = p_nodes;
end

h = grid_u.h;
xx = h*(0:Nx_p);
yy = h*(0:Ny_p);
[XX, YY] = meshgrid(xx, yy);

%% Bundle simulation metadata and parameters
info.u_inlet_func = u_inlet_func;
info.dt = dt;
info.grid_u = grid_u;
info.grid_v = grid_v;
info.grid_p = grid_p;
info.h = h;
info.block_info = block_info;
info.T_record = T_record;
info.t_simu = t_simu;
info.node_scale = node_scale;
info.t_solve = t_solve;
info.XX = XX;
info.YY = YY;
info.save_time = datetime('now');
 
%% Write data to compressed MAT file
% Use v7.3 format for large file support and compression
save(save_filename, 'info', '-v7.3');
fprintf('Data successfully saved to: %s\n', save_filename);

end


function p_nodes_all = interp_p(p_centers_all)
%   Converts pressure field from (Ny x Nx x Nt) cell centers to 
%   (Ny+1 x Nx+1 x Nt) grid nodes using vectorized operations
%   Boundary conditions:
%   - Top/Left/Bottom: zero gradient (∂p/∂n = 0)  
%   - Right: zero value (p = 0)
%
%   Input:
%       p_centers_all - Pressure at cell centers [Ny x Nx]
%   Output:
%       p_nodes_all   - Pressure at grid nodes [Ny+1 x Nx+1]

[Ny, Nx] = size(p_centers_all);
p_nodes_all = zeros(Ny + 1, Nx + 1, 'single');

% Interior nodes (2:Ny, 2:Nx) use average of 4 surrounding cell centers
p_nodes_all(2:Ny, 2:Nx) = 0.25 * (...
    p_centers_all(1:Ny-1, 1:Nx-1) + ...  % Top-left
    p_centers_all(1:Ny-1, 2:Nx) + ...  % Top-right
    p_centers_all(2:Ny,   1:Nx-1) + ...  % Bottom-left
    p_centers_all(2:Ny,   2:Nx));      % Bottom-right

% Bottom boundary (i=1): zero gradient ∂p/∂y = 0
p_nodes_all(1, 2:Nx) = p_nodes_all(2, 2:Nx);

% Top boundary (i=Ny+1): zero gradient ∂p/∂y = 0
p_nodes_all(Ny+1, 2:Nx) = p_nodes_all(Ny, 2:Nx);

% Left boundary (j=1): zero gradient ∂p/∂x = 0  
p_nodes_all(2:Ny, 1) = p_nodes_all(2:Ny, 2);

% Right boundary (j=Nx+1): zero value p = 0
p_nodes_all(:, Nx+1) = 0;

% Bottom-left: zero gradient from both directions
p_nodes_all(1, 1) = p_nodes_all(2, 2);

% Top-left: zero gradient from both directions
p_nodes_all(Ny+1, 1) = p_nodes_all(Ny, 2);

% Bottom-right & Top-right: zero value (right boundary dominates) no need to allocate

end
