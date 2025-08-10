% check_reynolds_number.m
function check_re(rho, U_char, mu, H, block_info, Re_limit)
%CHECK_RE Validate Reynolds number within safe operating range.
%   Calculates Reynolds number based on flow properties and domain geometry.
%   Terminates simulation with detailed error if Re exceeds specified limit.
%
%   Inputs:
%       rho        - Fluid density [kg/m³]
%       U_char     - Characteristic velocity [m/s]
%       mu         - Dynamic viscosity [Pa·s]
%       H          - Domain height [m]
%       block_info - Obstacle geometry data from init_block
%       Re_limit   - Maximum allowable Reynolds number
%
%   Characteristic Length Selection:
%   - No obstacles: uses domain height H
%   - With obstacles: uses height of largest obstacle
%   - Ensures numerical stability for given discretization
%
%   See also INIT_BLOCK, SETUP_SIMULATION.

%% Determine characteristic length scale
if isempty(block_info)
    % Pure channel flow: characteristic length is domain height
    L_char = H;
else
    % Flow with obstacles: use largest obstacle height
    max_obstacle_height = 0;
    for i = 1:numel(block_info)
        y_min_obstacle = block_info{i}.y_coords(1);
        y_max_obstacle = block_info{i}.y_coords(2);
        current_obstacle_height = y_max_obstacle - y_min_obstacle;
        if current_obstacle_height > max_obstacle_height
            max_obstacle_height = current_obstacle_height;
        end
    end
    L_char = max_obstacle_height;
end

% Handle degenerate case (thin obstacles)
if L_char <= 1e-9
    L_char = H;  % Fallback to domain height
end

%% Calculate Reynolds number
if mu > 0
    Re = rho * U_char * L_char / mu;
else
    Re = inf;  % Inviscid flow limit
end

%% Validate against stability limit
if Re > Re_limit
    % Create detailed error message for exceeded Reynolds number
    error_message = sprintf([...
        'Reynolds number exceeds safe operating limit.\n\n' ...
        '  - Calculated Re: %.0f\n' ...
        '  - Maximum allowed: %.0f\n\n' ...
        'To maintain numerical stability, consider:\n' ...
        '  - Reduce inlet velocity\n' ...
        '  - Reduce fluid density\n' ...
        '  - Increase fluid viscosity'], ...
        Re, Re_limit);
    
    error('MyApp:ReynoldsNumberExceeded', error_message);
end

% Function completes silently if Reynolds number is acceptable

end
