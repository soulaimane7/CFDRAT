function inlet_func = get_inlet_func(expr_str, ramp_opt, t_ramp, H, slip_opt)
%GET_INLET_FUNC Create time-dependent inlet velocity function handle.
%   Converts a user-defined velocity expression into a robust, vectorized
%   function handle with optional temporal ramping and automatic boundary
%   condition enforcement for no-slip walls.
%
%   Inputs:
%       expr_str  - String expression for velocity profile using variables t, y, H
%                   Example: '0.15*4*y*(H-y)/H^2' (parabolic profile)
%       ramp_opt  - Temporal ramp type: 'linear'|'smoothstep'|'none'
%       t_ramp    - Ramp duration , set to 0 for no ramping
%       H         - Domain height
%       slip_opt  - Wall boundary condition: 'no-slip'|'slip'
%
%   Output:
%       inlet_func - Function handle @(t,y) returning velocity 
%                   where t is time and y is vertical coordinate
%
%   Features:
%   - Automatic vectorization of user expressions
%   - Smooth temporal startup to avoid impulsive flow initiation
%   - Tukey window enforcement for no-slip boundary consistency
%   - Robust error handling for invalid expressions
%
%   Example:
%       % Parabolic profile with linear ramp-up
%       inlet_func = get_inlet_func('0.1*4*y*(H-y)/H^2', 'linear', 0.5, 0.3, 'no-slip');
%       u_inlet = inlet_func(1.0, [0:0.01:0.3]'); % Evaluate at t=1s
%
%   See also STR2FUNC, SMOOTHSTEP.

% Tukey window parameter for no-slip boundary enforcement
tukey_alpha = 0.2;                               % Transition region fraction (20%)

% Convert user expression to element-wise operations for vector compatibility
safe_expr_str = strrep(expr_str, '*', '.*');    % Enable element-wise multiplication
safe_expr_str = strrep(safe_expr_str, '/', './'); % Enable element-wise division
safe_expr_str = strrep(safe_expr_str, '^', '.^'); % Enable element-wise power

% Parse user expression into executable function handle
try
    profile_handle = str2func(['@(t,y,H)' safe_expr_str]);
catch ME
    error('Invalid inlet velocity expression: "%s".\nError: %s\nEnsure it only uses variables t, y, and H.', safe_expr_str, ME.message);
end

% Configure temporal ramp-up function
if t_ramp <= 0
    ramp_factor_handle = @(t) 1.0;               % No ramping
else
    switch ramp_opt
        case 'linear' 
            ramp_factor_handle = @(t) min(1.0, t / t_ramp);        % Linear ramp
        case 'smoothstep' 
            ramp_factor_handle = @(t) smoothstep_internal(t / t_ramp); % Smooth ramp
        otherwise
            ramp_factor_handle = @(t) 1.0;       % Default: no ramping
    end
end

% Combine spatial profile with temporal ramping
% Element-wise multiplication handles scalar ramp with vector profile
base_inlet_func = @(t, y) profile_handle(t, y, H) .* ramp_factor_handle(t);

% Apply boundary condition enforcement for no-slip walls
if strcmp(slip_opt, 'no-slip')
    % Evaluate user function at wall boundaries over time range
    vel_y0 = max(abs(base_inlet_func(0:0.1:10, 0)));       % Bottom wall velocity
    vel_yH = max(abs(base_inlet_func(0:0.1:10, H)));       % Top wall velocity
    
    % Check boundary condition consistency
    if abs(vel_y0) > 1e-3 || abs(vel_yH) > 1e-3
        % User profile violates no-slip condition - apply corrective window
        fprintf('\n[WARNING] Inlet profile is non-zero at the walls for a "no-slip" case.\n');
        fprintf('          Solver will apply a flat-top window to enforce zero velocity at walls,\n');
        fprintf('          preserving the central 80%% of the profile.\n');
        fprintf('          User-defined velocity at y=0 was %.4f, at y=H was %.4f.\n\n', vel_y0, vel_yH);
        
        % Apply Tukey window to enforce zero velocity at walls
        inlet_func = @(t, y) base_inlet_func(t, y) .* tukey_window_vectorized(y, H, tukey_alpha);
    else
        % User profile already satisfies no-slip condition
        inlet_func = base_inlet_func;
    end
else
    % For slip condition, use user profile directly
    inlet_func = base_inlet_func;
end

end


function W = tukey_window_vectorized(y, H, alpha)
%TUKEY_WINDOW_VECTORIZED Flat-top cosine window for boundary enforcement.
%   Applies a Tukey window that smoothly transitions from 1 (interior)
%   to 0 (walls) to enforce no-slip boundary conditions.
    y_norm = y / H;                              % Normalize coordinate to [0,1]
    W = ones(size(y_norm));                      % Initialize window values
    transition_length = alpha / 2;               % Half-width of transition region
    
    % Apply cosine taper near bottom wall (y=0)
    left_mask = y_norm <= transition_length;
    if any(left_mask)
        W(left_mask) = 0.5 * (1 + cos(pi * (y_norm(left_mask) / transition_length - 1)));
    end
    
    % Apply cosine taper near top wall (y=H)
    right_mask = y_norm >= (1 - transition_length);
    if any(right_mask)
        W(right_mask) = 0.5 * (1 + cos(pi * (y_norm(right_mask) - (1 - transition_length)) / transition_length));
    end
end


function val = smoothstep_internal(x)
%SMOOTHSTEP_INTERNAL Smooth Hermite interpolation for ramping.
%   Provides C2-continuous transition from 0 to 1 using cubic Hermite polynomial.
    x_clamped = max(0, min(1, x));               % Clamp input to [0,1] range
    val = x_clamped.^2 .* (3.0 - 2.0 * x_clamped); % Cubic Hermite polynomial
end
